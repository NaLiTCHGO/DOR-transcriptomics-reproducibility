options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "../.."), winslash = "/")
} else {
  getwd()
}
root <- Sys.getenv("DOR_PROJECT_ROOT", unset = default_root)
m02 <- file.path(root, "05_analysis_steps/M02_COHORT_QC/runs/20260813_M02_B1_GSE274832")
run <- file.path(root, "05_analysis_steps/M03_WITHIN_COHORT_EFFECTS/runs/20260814_M03_B1_GSE274832")
out <- file.path(run, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cfg <- read.csv(file.path(root, "04_code/configs/GSE274832_SAMPLES.csv"),
                check.names = FALSE)
cfg$condition <- factor(cfg$phenotype, levels = c("NOR", "DOR"))
rownames(cfg) <- cfg$srr
files <- file.path(m02, "results/salmon", cfg$srr, "quant.sf")
names(files) <- cfg$srr
stopifnot(all(file.exists(files)), nrow(cfg) == 6L, all(table(cfg$condition) == 3L))

txmap <- read.delim(gzfile(file.path(root, "03_data/references/GENCODE_v50_GRCh38p14/gencode.v50.tx2gene.tsv.gz")),
                    check.names = FALSE)
txmap <- txmap[!duplicated(txmap$transcript_id_version), ]
tx2gene <- txmap[, c("transcript_id_version", "gene_id_version")]
gene_annot <- txmap[!duplicated(txmap$gene_id_version),
                    c("gene_id_version", "gene_id", "gene_name")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# Align the frozen sample table explicitly to the imported abundance columns,
# then create one complete DESeq2 object. Reduced sensitivity fits subset this
# object, keeping the abundance matrix and colData dimensions synchronized.
cfg <- cfg[match(colnames(txi$counts), cfg$srr), , drop = FALSE]
stopifnot(!anyNA(cfg$srr), identical(cfg$srr, colnames(txi$counts)))
rownames(cfg) <- cfg$srr
dds_all <- DESeqDataSetFromTximport(txi, colData = cfg, design = ~ condition)

fit_deseq <- function(keep_samples) {
  dds <- dds_all[, keep_samples]
  dds$condition <- droplevels(dds$condition)
  keep_gene <- rowSums(counts(dds) >= 10) >= min(3L, ncol(dds))
  dds <- dds[keep_gene, ]
  dds <- DESeq(dds, quiet = TRUE)
  rr <- results(dds, contrast = c("condition", "DOR", "NOR"), alpha = 0.05)
  list(dds = dds, result = rr)
}

primary <- fit_deseq(cfg$srr)
rr <- as.data.frame(primary$result)
rr$gene_id_version <- rownames(rr)
rr <- merge(gene_annot, rr, by = "gene_id_version", all.y = TRUE, sort = FALSE)
rr$direction <- ifelse(is.na(rr$log2FoldChange), "NA",
                       ifelse(rr$log2FoldChange > 0, "UP_IN_DOR",
                              ifelse(rr$log2FoldChange < 0, "DOWN_IN_DOR", "ZERO")))
rr$abs_log2FoldChange <- abs(rr$log2FoldChange)
rr <- rr[order(rr$pvalue, -rr$abs_log2FoldChange, na.last = TRUE), ]
write.csv(rr, file.path(out, "GSE274832_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv"), row.names = FALSE)

norm_counts <- counts(primary$dds, normalized = TRUE)
norm_df <- merge(gene_annot, data.frame(gene_id_version = rownames(norm_counts), norm_counts,
                                        check.names = FALSE), by = "gene_id_version", all.y = TRUE)
write.csv(norm_df, gzfile(file.path(out, "GSE274832_DESEQ2_NORMALIZED_COUNTS.csv.gz")), row.names = FALSE)

# Leave-one-sample-out sensitivity: no post-hoc exclusion; all six reduced fits are diagnostic.
loso_lfc <- list()
loso_summary <- list()
for (held in cfg$srr) {
  kept <- setdiff(cfg$srr, held)
  fit <- fit_deseq(kept)
  z <- as.data.frame(fit$result)
  z$gene_id_version <- rownames(z)
  col <- paste0("lfc_without_", held)
  loso_lfc[[col]] <- z[, c("gene_id_version", "log2FoldChange")]
  names(loso_lfc[[col]])[2] <- col
  common <- intersect(rr$gene_id_version[is.finite(rr$log2FoldChange)],
                      z$gene_id_version[is.finite(z$log2FoldChange)])
  full_v <- rr$log2FoldChange[match(common, rr$gene_id_version)]
  reduced_v <- z$log2FoldChange[match(common, z$gene_id_version)]
  loso_summary[[held]] <- data.frame(
    held_out_srr = held,
    held_out_phenotype = as.character(cfg[held, "phenotype"]),
    n_dor = sum(cfg[kept, "phenotype"] == "DOR"),
    n_nor = sum(cfg[kept, "phenotype"] == "NOR"),
    common_tested_genes = length(common),
    spearman_lfc_vs_full = cor(full_v, reduced_v, method = "spearman", use = "complete.obs"),
    same_direction_fraction = mean(sign(full_v) == sign(reduced_v), na.rm = TRUE)
  )
}
loso <- Reduce(function(x, y) merge(x, y, by = "gene_id_version", all = TRUE), loso_lfc)
full_lfc <- rr[, c("gene_id_version", "gene_id", "gene_name", "log2FoldChange", "pvalue", "padj")]
names(full_lfc)[4] <- "full_log2FoldChange"
loso <- merge(full_lfc, loso, by = "gene_id_version", all.x = TRUE)
lfc_cols <- grep("^lfc_without_", names(loso), value = TRUE)
loso$loso_median_log2FoldChange <- apply(loso[, lfc_cols], 1, median, na.rm = TRUE)
loso$loso_min_log2FoldChange <- apply(loso[, lfc_cols], 1, min, na.rm = TRUE)
loso$loso_max_log2FoldChange <- apply(loso[, lfc_cols], 1, max, na.rm = TRUE)
loso$loso_same_direction_fraction <- vapply(seq_len(nrow(loso)), function(i) {
  v <- as.numeric(loso[i, lfc_cols])
  v <- v[is.finite(v)]
  if (!length(v) || !is.finite(loso$full_log2FoldChange[i])) return(NA_real_)
  mean(sign(v) == sign(loso$full_log2FoldChange[i]))
}, numeric(1))
write.csv(loso, file.path(out, "GSE274832_LEAVE_ONE_SAMPLE_OUT_GENE_STABILITY.csv"), row.names = FALSE)
loso_qc <- do.call(rbind, loso_summary)
rownames(loso_qc) <- NULL
write.csv(loso_qc, file.path(out, "GSE274832_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv"), row.names = FALSE)

# Independent validation against the submitter-provided FPKM matrix.
official_path <- file.path(root, "03_data/intermediate/GSE274832_gene_fpkm.xls")
official <- read.delim(official_path, check.names = FALSE)
official_samples <- c("DOR1", "DOR2", "DOR3", "NOR1", "NOR2", "NOR3")
official_to_srr <- setNames(cfg$srr, official_samples)
tpm_gene <- read.delim(gzfile(file.path(m02, "results/expression_qc/GSE274832_GENE_TPM.tsv.gz")),
                       check.names = FALSE)
sample_validation <- lapply(official_samples, function(label) {
  srr <- official_to_srr[[label]]
  common <- merge(official[, c("gene_id", label)], tpm_gene[, c("gene_id", srr)], by = "gene_id")
  data.frame(official_label = label, srr = srr, common_genes = nrow(common),
             spearman_log_expression = cor(log2(common[[label]] + 1),
                                           log2(common[[srr]] + 1), method = "spearman"),
             pearson_log_expression = cor(log2(common[[label]] + 1),
                                          log2(common[[srr]] + 1), method = "pearson"))
})
sample_validation <- do.call(rbind, sample_validation)
write.csv(sample_validation, file.path(out, "GSE274832_OFFICIAL_FPKM_VALIDATION.csv"), row.names = FALSE)

official$official_log2FC <- log2((rowMeans(official[, c("DOR1", "DOR2", "DOR3")]) + 0.1) /
                                 (rowMeans(official[, c("NOR1", "NOR2", "NOR3")]) + 0.1))
effect_validation <- merge(rr[, c("gene_id", "gene_name", "log2FoldChange")],
                           official[, c("gene_id", "official_log2FC")], by = "gene_id")
effect_validation <- effect_validation[is.finite(effect_validation$log2FoldChange) &
                                       is.finite(effect_validation$official_log2FC), ]
write.csv(effect_validation, file.path(out, "GSE274832_OFFICIAL_VS_REPROCESSED_EFFECTS.csv"), row.names = FALSE)
effect_spearman <- cor(effect_validation$log2FoldChange, effect_validation$official_log2FC,
                       method = "spearman")
effect_direction <- mean(sign(effect_validation$log2FoldChange) == sign(effect_validation$official_log2FC))

# Figures.
plot_df <- rr[is.finite(rr$log2FoldChange) & is.finite(rr$pvalue), ]
plot_df$display <- "Other"
plot_df$display[!is.na(plot_df$padj) & plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) >= 1] <- "FDR<0.05 and |LFC|>=1"
p <- ggplot(plot_df, aes(log2FoldChange, -log10(pmax(pvalue, 1e-300)), color = display)) +
  geom_point(alpha = 0.55, size = 1) +
  scale_color_manual(values = c("Other" = "grey70", "FDR<0.05 and |LFC|>=1" = "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
  labs(title = "GSE274832: DOR minus NOR", x = "DESeq2 log2 fold change",
       y = "-log10 nominal p", color = NULL) + theme_bw(base_size = 12)
ggsave(file.path(out, "GSE274832_DESEQ2_VOLCANO.png"), p, width = 7.5, height = 6, dpi = 240)

p2 <- ggplot(loso_qc, aes(reorder(held_out_srr, spearman_lfc_vs_full), spearman_lfc_vs_full,
                           fill = held_out_phenotype)) +
  geom_col() + coord_cartesian(ylim = c(0, 1)) + coord_flip() +
  scale_fill_manual(values = c(DOR = "#D55E00", NOR = "#0072B2")) +
  labs(title = "Leave-one-sample-out effect stability", x = "Held-out sample",
       y = "Spearman correlation with full-model LFC", fill = NULL) + theme_bw(base_size = 12)
ggsave(file.path(out, "GSE274832_LOSO_GLOBAL_STABILITY.png"), p2, width = 7.5, height = 5, dpi = 240)

vst_mat <- assay(vst(primary$dds, blind = TRUE))
top_ids <- head(rr$gene_id_version[order(rr$pvalue, na.last = NA)], 30)
top_ids <- intersect(top_ids, rownames(vst_mat))
hm <- vst_mat[top_ids, , drop = FALSE]
rownames(hm) <- gene_annot$gene_name[match(rownames(hm), gene_annot$gene_id_version)]
ann <- data.frame(Phenotype = cfg[colnames(hm), "phenotype"], row.names = colnames(hm))
png(file.path(out, "GSE274832_TOP30_NOMINAL_GENES_HEATMAP.png"), width = 1800, height = 2200, res = 240)
pheatmap(hm, scale = "row", annotation_col = ann, show_colnames = TRUE,
         main = "Top 30 nominal genes (exploratory; row-scaled VST)")
dev.off()

n_fdr <- sum(!is.na(rr$padj) & rr$padj < 0.05)
n_fdr_lfc <- sum(!is.na(rr$padj) & rr$padj < 0.05 & abs(rr$log2FoldChange) >= 1)
n_nom_lfc <- sum(!is.na(rr$pvalue) & rr$pvalue < 0.05 & abs(rr$log2FoldChange) >= 1)
stable_fraction <- mean(loso$loso_same_direction_fraction >= 5 / 6, na.rm = TRUE)

summary <- c(
  "# GSE274832 M03 Within-Cohort Effect Summary", "",
  "## Model", "",
  "- Primary design: `~ condition`, contrast DOR minus NOR; 3 DOR vs 3 NOR.",
  "- Input: Salmon transcript estimates summarized with tximport `lengthScaledTPM`; inference by DESeq2.",
  sprintf("- Tested genes after count filter: %s.", format(nrow(rr), big.mark = ",")),
  "- All six M02-retained samples are used; no PCA-driven post-hoc deletion.", "",
  "## Results and stability", "",
  sprintf("- FDR < 0.05 genes: %d; FDR < 0.05 and |log2FC| >= 1: %d.", n_fdr, n_fdr_lfc),
  sprintf("- Nominal p < 0.05 and |log2FC| >= 1: %d (exploratory only).", n_nom_lfc),
  sprintf("- LOSO LFC-vs-full Spearman range: %.3f–%.3f.", min(loso_qc$spearman_lfc_vs_full), max(loso_qc$spearman_lfc_vs_full)),
  sprintf("- Fraction of tested genes retaining the full-model direction in at least 5/6 LOSO fits: %.3f.", stable_fraction),
  "", "## Independent submitter-table validation", "",
  sprintf("- Per-sample log-expression Spearman range (reprocessed TPM vs official FPKM): %.3f–%.3f.",
          min(sample_validation$spearman_log_expression), max(sample_validation$spearman_log_expression)),
  sprintf("- Gene-effect Spearman (DESeq2 LFC vs official-FPKM mean LFC): %.3f; direction agreement %.3f.",
          effect_spearman, effect_direction),
  "", "## Interpretation boundary", "",
  "- This n=3-per-group result is a cohort-specific effect estimate, not a standalone biomarker claim.",
  "- Any nominal gene list is exploratory. Final claims require direction/rank and pathway replication in independent cohorts.",
  "- Substantial within-DOR heterogeneity remains a prespecified limitation and is carried into cross-cohort synthesis."
)
writeLines(summary, file.path(out, "GSE274832_M03_EFFECT_SUMMARY.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(run, "logs/R_SESSION_INFO.txt"), useBytes = TRUE)
writeLines(paste("PASS_WITH_LIMITATION", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
           file.path(run, "GSE274832_M03_PASS_WITH_LIMITATION.txt"), useBytes = TRUE)
cat("PASS_WITH_LIMITATION: GSE274832 M03 primary and LOSO analyses completed\n")
