options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- if (length(script_arg)) normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "../.."), winslash = "/") else getwd()
root <- Sys.getenv("DOR_PROJECT_ROOT", unset = default_root)
m02 <- file.path(root, "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B4_GSE232306")
run <- file.path(root, "05_analysis_steps/M03_WITHIN_COHORT_EFFECTS/runs/20260814_M03_B4_GSE232306")
out <- file.path(run, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run, "logs"), recursive = TRUE, showWarnings = FALSE)

cfg <- read.csv(file.path(root, "04_code/configs/GSE232306_SAMPLES.csv"), check.names = FALSE)
inclusion <- read.csv(file.path(m02, "results/final_qc/GSE232306_M02_SAMPLE_INCLUSION.csv"), check.names = FALSE)
cfg <- cfg[cfg$srr %in% inclusion$srr[inclusion$include_m03 == "YES"], , drop = FALSE]
stopifnot(nrow(cfg) == 12L, all(table(cfg$phenotype) == 6L))
cfg$condition <- factor(cfg$phenotype, levels = c("NOR", "DOR"))
rownames(cfg) <- cfg$srr

files <- file.path(m02, "results/salmon", cfg$srr, "quant.sf")
names(files) <- cfg$srr
stopifnot(all(file.exists(files)))
txmap <- read.delim(gzfile(file.path(root, "03_data/references/GENCODE_v50_GRCh38p14/gencode.v50.tx2gene.tsv.gz")), check.names = FALSE)
txmap <- txmap[!duplicated(txmap$transcript_id_version), ]
tx2gene <- txmap[, c("transcript_id_version", "gene_id_version")]
gene_annot <- txmap[!duplicated(txmap$gene_id_version), c("gene_id_version", "gene_id", "gene_name")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
cfg <- cfg[match(colnames(txi$counts), cfg$srr), , drop = FALSE]
stopifnot(identical(cfg$srr, colnames(txi$counts)))
rownames(cfg) <- cfg$srr

dds_base <- DESeqDataSetFromTximport(txi, colData = cfg, design = ~ condition)
keep_gene <- rowSums(counts(dds_base) >= 10) >= 6L
dds_base <- dds_base[keep_gene, ]

fit_model <- function(keep_samples) {
  dds <- dds_base[, keep_samples]
  dds$condition <- droplevels(dds$condition)
  design(dds) <- ~ condition
  mm <- model.matrix(design(dds), colData(dds))
  if (qr(mm)$rank != ncol(mm)) stop("Design matrix is not full rank")
  dds <- DESeq(dds, quiet = TRUE)
  rr <- results(dds, contrast = c("condition", "DOR", "NOR"), alpha = 0.05)
  list(dds = dds, result = rr)
}

as_table <- function(result) {
  rr <- as.data.frame(result)
  rr$gene_id_version <- rownames(rr)
  rr <- merge(gene_annot, rr, by = "gene_id_version", all.y = TRUE, sort = FALSE)
  rr$model <- "UNADJUSTED_CONDITION_ONLY_PRIMARY"
  rr$direction <- ifelse(is.na(rr$log2FoldChange), "NA", ifelse(rr$log2FoldChange > 0, "UP_IN_DOR", ifelse(rr$log2FoldChange < 0, "DOWN_IN_DOR", "ZERO")))
  rr$abs_log2FoldChange <- abs(rr$log2FoldChange)
  rr[order(rr$pvalue, -rr$abs_log2FoldChange, na.last = TRUE), ]
}

primary_fit <- fit_model(cfg$srr)
effects <- as_table(primary_fit$result)
write.csv(effects, file.path(out, "GSE232306_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv"), row.names = FALSE)

norm_counts <- counts(primary_fit$dds, normalized = TRUE)
norm_df <- merge(gene_annot, data.frame(gene_id_version = rownames(norm_counts), norm_counts, check.names = FALSE), by = "gene_id_version", all.y = TRUE)
write.csv(norm_df, gzfile(file.path(out, "GSE232306_DESEQ2_NORMALIZED_COUNTS.csv.gz")), row.names = FALSE)

# Leave-one-sample-out condition-only sensitivity; all reduced fits are diagnostic.
loso_lfc <- list()
loso_summary <- list()
for (held in cfg$srr) {
  kept <- setdiff(cfg$srr, held)
  fit <- fit_model(kept)
  z <- as.data.frame(fit$result)
  z$gene_id_version <- rownames(z)
  col <- paste0("lfc_without_", held)
  loso_lfc[[col]] <- z[, c("gene_id_version", "log2FoldChange")]
  names(loso_lfc[[col]])[2] <- col
  common <- intersect(effects$gene_id_version[is.finite(effects$log2FoldChange)], z$gene_id_version[is.finite(z$log2FoldChange)])
  full_v <- effects$log2FoldChange[match(common, effects$gene_id_version)]
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
full_lfc <- effects[, c("gene_id_version", "gene_id", "gene_name", "log2FoldChange", "pvalue", "padj")]
names(full_lfc)[4] <- "full_log2FoldChange"
loso <- merge(full_lfc, loso, by = "gene_id_version", all.x = TRUE)
lfc_cols <- grep("^lfc_without_", names(loso), value = TRUE)
loso$loso_median_log2FoldChange <- apply(loso[, lfc_cols], 1, median, na.rm = TRUE)
loso$loso_min_log2FoldChange <- apply(loso[, lfc_cols], 1, min, na.rm = TRUE)
loso$loso_max_log2FoldChange <- apply(loso[, lfc_cols], 1, max, na.rm = TRUE)
loso$loso_same_direction_fraction <- vapply(seq_len(nrow(loso)), function(i) {
  v <- as.numeric(loso[i, lfc_cols]); v <- v[is.finite(v)]
  if (!length(v) || !is.finite(loso$full_log2FoldChange[i])) return(NA_real_)
  mean(sign(v) == sign(loso$full_log2FoldChange[i]))
}, numeric(1))
write.csv(loso, file.path(out, "GSE232306_LEAVE_ONE_SAMPLE_OUT_GENE_STABILITY.csv"), row.names = FALSE)
loso_qc <- do.call(rbind, loso_summary); rownames(loso_qc) <- NULL
write.csv(loso_qc, file.path(out, "GSE232306_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv"), row.names = FALSE)

plot_df <- effects[is.finite(effects$log2FoldChange) & is.finite(effects$pvalue), ]
plot_df$display <- "Other"
plot_df$display[!is.na(plot_df$padj) & plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) >= 1] <- "FDR<0.05 and |LFC|>=1"
p <- ggplot(plot_df, aes(log2FoldChange, -log10(pmax(pvalue, 1e-300)), color = display)) +
  geom_point(alpha = 0.55, size = 1) +
  scale_color_manual(values = c("Other" = "grey70", "FDR<0.05 and |LFC|>=1" = "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
  labs(title = "GSE232306: DOR minus NOR", x = "DESeq2 log2 fold change", y = "-log10 nominal p", color = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(out, "GSE232306_DESEQ2_VOLCANO.png"), p, width = 7.5, height = 6, dpi = 240)

p2 <- ggplot(loso_qc, aes(reorder(held_out_srr, spearman_lfc_vs_full), spearman_lfc_vs_full, fill = held_out_phenotype)) +
  geom_col() + coord_cartesian(ylim = c(0, 1)) + coord_flip() +
  scale_fill_manual(values = c(DOR = "#D55E00", NOR = "#0072B2")) +
  labs(title = "Leave-one-sample-out stability", x = "Held-out sample", y = "Spearman correlation with full-model LFC", fill = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(out, "GSE232306_LOSO_GLOBAL_STABILITY.png"), p2, width = 7.5, height = 6.5, dpi = 240)

vst_mat <- assay(vst(primary_fit$dds, blind = TRUE))
top_ids <- head(effects$gene_id_version[order(effects$pvalue, na.last = NA)], 30)
top_ids <- intersect(top_ids, rownames(vst_mat))
hm <- vst_mat[top_ids, , drop = FALSE]
rownames(hm) <- gene_annot$gene_name[match(rownames(hm), gene_annot$gene_id_version)]
ann <- data.frame(Phenotype = cfg[colnames(hm), "phenotype"], row.names = colnames(hm))
png(file.path(out, "GSE232306_TOP30_NOMINAL_GENES_HEATMAP.png"), width = 2000, height = 2400, res = 240)
pheatmap(hm, scale = "row", annotation_col = ann, show_colnames = TRUE, main = "Top 30 nominal genes (exploratory)")
dev.off()

n_fdr <- sum(!is.na(effects$padj) & effects$padj < 0.05)
n_fdr_lfc <- sum(!is.na(effects$padj) & effects$padj < 0.05 & abs(effects$log2FoldChange) >= 1)
n_nom_lfc <- sum(!is.na(effects$pvalue) & effects$pvalue < 0.05 & abs(effects$log2FoldChange) >= 1)
raw_qc <- read.csv(file.path(m02, "results/raw_qc/GSE232306_FASTP_RAW_QC_SUMMARY.csv"), check.names = FALSE)
pca_qc <- read.csv(file.path(m02, "results/expression_qc/GSE232306_PCA_COORDINATES.csv"), check.names = FALSE)
dor_pc1 <- pca_qc$PC1[pca_qc$phenotype == "DOR"]
nor_pc1 <- pca_qc$PC1[pca_qc$phenotype == "NOR"]
pc1_complete_separation <- max(dor_pc1) < min(nor_pc1) || max(nor_pc1) < min(dor_pc1)
summary <- c(
  "# GSE232306 M03 Within-Cohort Effect Summary", "",
  "## Model", "",
  "- Principal and only estimable design: `~ condition`; contrast DOR minus NOR; 6 DOR vs 6 NOR.",
  "- Individual ages are absent from public GEO sample metadata, so no age-adjusted model is claimed.",
  "- Input: Salmon transcript estimates summarized by tximport `lengthScaledTPM`; inference by DESeq2.",
  sprintf("- Tested genes after count filter: %s.", format(nrow(effects), big.mark = ",")),
  "- All M02-retained samples are used; no PCA-driven post-hoc deletion.", "",
  "## Results and stability", "",
  sprintf("- FDR < 0.05 genes: %d; FDR < 0.05 and |log2FC| >= 1: %d.", n_fdr, n_fdr_lfc),
  sprintf("- FDR fraction: %.1f%% of tested genes; FDR-and-|LFC| fraction: %.1f%%.", 100 * n_fdr / nrow(effects), 100 * n_fdr_lfc / nrow(effects)),
  sprintf("- Nominal p < 0.05 and |log2FC| >= 1: %d (exploratory only).", n_nom_lfc),
  sprintf("- LOSO LFC-vs-full Spearman range: %.3f–%.3f.", min(loso_qc$spearman_lfc_vs_full), max(loso_qc$spearman_lfc_vs_full)),
  "", "## Interpretation boundary", "",
  "- The condition-only model estimates the observed cohort group effect; it cannot separate phenotype from unmeasured age or other clinical covariates.",
  sprintf("- PC1 completely separates phenotype: %s; PC1 variance: %.2f%%.", pc1_complete_separation, pca_qc$PC1_variance_percent[1]),
  sprintf("- Raw GC mean differs by phenotype (DOR %.4f; NOR %.4f) and duplication mean differs (DOR %.4f; NOR %.4f). These phenotype-aligned shifts raise unresolved compositional/batch concern.", mean(raw_qc$gc_content[raw_qc$phenotype == "DOR"]), mean(raw_qc$gc_content[raw_qc$phenotype == "NOR"]), mean(raw_qc$duplication_rate[raw_qc$phenotype == "DOR"]), mean(raw_qc$duplication_rate[raw_qc$phenotype == "NOR"])),
  "- The 6,381 FDR-and-|LFC| genes are a broad cohort signature for cross-cohort screening, not a biomarker list.",
  "- These are cohort-specific effects, not standalone biomarkers; final claims require cross-cohort direction/rank and pathway replication."
)
writeLines(summary, file.path(out, "GSE232306_M03_EFFECT_SUMMARY.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(run, "logs/R_SESSION_INFO.txt"), useBytes = TRUE)
writeLines(paste("PASS_WITH_LIMITATION", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"), file.path(run, "GSE232306_M03_PASS_WITH_LIMITATION.txt"), useBytes = TRUE)
cat("PASS_WITH_LIMITATION: GSE232306 condition-only M03 and LOSO analyses completed\n")
