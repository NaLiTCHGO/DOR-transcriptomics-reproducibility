#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop(sprintf("Missing argument: %s", flag))
  args[[idx + 1L]]
}

project_root <- normalizePath(get_arg("--project-root"), winslash = "/", mustWork = TRUE)
run_dir <- normalizePath(get_arg("--run-dir"), winslash = "/", mustWork = FALSE)
results_dir <- file.path(run_dir, "results")
plots_dir <- file.path(results_dir, "plots")
logs_dir <- file.path(run_dir, "logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("GSE274832", "GSE193136", "GSE232306")
cohort_short <- c(GSE274832 = "274832", GSE193136 = "193136", GSE232306 = "232306")
meta_path <- file.path(
  project_root,
  "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_GENE_LEVEL_META_ANALYSIS.csv"
)
gmt_paths <- c(
  HALLMARK = file.path(project_root, "03_data/reference/gene_sets/MSigDB_v2026.1_Hs/h.all.v2026.1.Hs.symbols.gmt"),
  REACTOME = file.path(project_root, "03_data/reference/gene_sets/MSigDB_v2026.1_Hs/c2.cp.reactome.v2026.1.Hs.symbols.gmt")
)

for (path in c(meta_path, gmt_paths)) {
  if (!file.exists(path)) stop(sprintf("Required input missing: %s", path))
}

meta <- fread(meta_path)
required_meta <- c("gene_id", "gene_name", paste0("z_", cohorts), paste0("beta_", cohorts))
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0L) stop(sprintf("M04 meta table missing: %s", paste(missing_meta, collapse = ", ")))
if (anyDuplicated(meta$gene_id)) stop("M04 meta table has duplicate gene_id")

mapped <- meta[!is.na(gene_name) & !grepl("^ENSG[0-9]+$", gene_name)]
mapped_total <- nrow(mapped)
duplicate_symbol_rows <- sum(duplicated(mapped$gene_name))
duplicate_symbol_n <- uniqueN(mapped[duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE), gene_name])
setorder(mapped, gene_name, gene_id)
mapped <- mapped[!duplicated(gene_name)]
if (nrow(mapped) < 10000L) stop(sprintf("Mapped symbol universe below gate: %d", nrow(mapped)))
if (anyDuplicated(mapped$gene_name)) stop("Symbol collapse failed")

rank_table <- mapped[, c("gene_id", "gene_name", paste0("z_", cohorts)), with = FALSE]
setnames(rank_table, "gene_name", "gene_symbol")
fwrite(rank_table, file.path(results_dir, "M05_RANK_INPUT_SYMBOLS.csv"))

rank_audit <- rbindlist(lapply(cohorts, function(cohort) {
  z <- rank_table[[paste0("z_", cohort)]]
  data.table(
    cohort = cohort,
    m04_common_gene_rows = nrow(meta),
    rows_with_mapped_symbol_before_collapse = mapped_total,
    excluded_ensg_only_or_missing = nrow(meta) - mapped_total,
    duplicated_symbol_rows_removed = duplicate_symbol_rows,
    duplicated_symbol_names = duplicate_symbol_n,
    unique_rank_symbols = nrow(rank_table),
    finite_rank_values = sum(is.finite(z)),
    min_signed_wald_z = min(z, na.rm = TRUE),
    max_signed_wald_z = max(z, na.rm = TRUE),
    comparison = "DOR-minus-NOR",
    prefilter = "NONE_FULL_COMMON_MAPPED_UNIVERSE"
  )
}))
if (any(rank_audit$finite_rank_values != rank_audit$unique_rank_symbols)) stop("Non-finite rank metric")
fwrite(rank_audit, file.path(results_dir, "M05_RANK_INPUT_AUDIT.csv"))

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  parts <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) stop(sprintf("Malformed GMT: %s", path))
  names_vec <- vapply(parts, `[[`, character(1), 1L)
  if (anyDuplicated(names_vec)) stop(sprintf("Duplicate pathway names in GMT: %s", path))
  sets <- lapply(parts, function(x) unique(x[-c(1L, 2L)]))
  names(sets) <- names_vec
  sets
}

universe <- rank_table$gene_symbol
gene_sets <- list()
geneset_audit <- rbindlist(lapply(names(gmt_paths), function(collection) {
  source_sets <- read_gmt(gmt_paths[[collection]])
  intersected <- lapply(source_sets, intersect, y = universe)
  eligible <- intersected[lengths(intersected) >= 15L & lengths(intersected) <= 500L]
  gene_sets[[collection]] <<- eligible
  data.table(
    collection = collection,
    release = "MSigDB_v2026.1.Hs",
    local_path = substring(gmt_paths[[collection]], nchar(project_root) + 2L),
    source_sets = length(source_sets),
    eligible_sets_15_to_500 = length(eligible),
    sets_below_15 = sum(lengths(intersected) < 15L),
    sets_above_500 = sum(lengths(intersected) > 500L),
    mapped_universe_symbols = length(universe),
    collection_gene_coverage = length(intersect(unique(unlist(source_sets, use.names = FALSE)), universe)),
    min_eligible_size = min(lengths(eligible)),
    max_eligible_size = max(lengths(eligible))
  )
}))
if (geneset_audit[collection == "HALLMARK", eligible_sets_15_to_500] < 45L) stop("Hallmark eligible sets below gate")
if (geneset_audit[collection == "REACTOME", eligible_sets_15_to_500] < 900L) stop("Reactome eligible sets below gate")
fwrite(geneset_audit, file.path(results_dir, "M05_GENESET_AUDIT.csv"))

run_fgsea <- function(cohort, collection) {
  z_col <- paste0("z_", cohort)
  stats <- rank_table[[z_col]]
  names(stats) <- rank_table$gene_symbol
  lexical_rank <- rank(names(stats), ties.method = "first")
  jitter <- ((lexical_rank / (length(stats) + 1)) - 0.5) * 1e-12
  stats <- sort(stats + jitter, decreasing = TRUE)
  set.seed(20260814)
  ans <- as.data.table(fgseaMultilevel(
    pathways = gene_sets[[collection]],
    stats = stats,
    minSize = 15,
    maxSize = 500,
    eps = 1e-10,
    scoreType = "std",
    nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  if (nrow(ans) != length(gene_sets[[collection]])) {
    stop(sprintf("fgsea returned %d/%d pathways for %s %s", nrow(ans), length(gene_sets[[collection]]), cohort, collection))
  }
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  ans[, `:=`(
    cohort = cohort,
    collection = collection,
    rank_metric = "SIGNED_WALD_Z_FULL_COMMON_MAPPED_UNIVERSE",
    msigdb_release = "v2026.1.Hs"
  )]
  preferred_order <- c("collection", "cohort", "pathway", "pval", "padj", "ES", "NES", "log2err", "size", "leadingEdge", "rank_metric", "msigdb_release")
  setcolorder(ans, intersect(preferred_order, names(ans)))
  ans[]
}

fgsea_results <- rbindlist(lapply(names(gene_sets), function(collection) {
  rbindlist(lapply(cohorts, run_fgsea, collection = collection))
}))
fwrite(fgsea_results, file.path(results_dir, "M05_FGSEA_ALL_RESULTS.csv"))
fgsea_warning_audit <- fgsea_results[, .(
  pathway_n = .N,
  pvalue_at_eps_1e_10_n = sum(pval <= 1e-10, na.rm = TRUE),
  log2err_missing_n = sum(is.na(log2err)),
  minimum_reported_pvalue = min(pval, na.rm = TRUE)
), by = .(collection, cohort)]
fgsea_warning_audit[, interpretation := "Values at eps are upper-bounded/capped for ranking and FDR; exact smaller p-values are not claimed"]
fwrite(fgsea_warning_audit, file.path(results_dir, "M05_FGSEA_NUMERICAL_WARNING_AUDIT.csv"))

build_convergence <- function(collection_name) {
  sub <- fgsea_results[collection == collection_name]
  wide <- dcast(
    sub,
    pathway + size ~ cohort,
    value.var = c("NES", "pval", "padj")
  )
  nes_cols <- paste0("NES_", cohorts)
  p_cols <- paste0("pval_", cohorts)
  padj_cols <- paste0("padj_", cohorts)
  nes <- as.matrix(wide[, ..nes_cols])
  pmat <- as.matrix(wide[, ..p_cols])
  qmat <- as.matrix(wide[, ..padj_cols])
  positive_n <- rowSums(nes > 0)
  negative_n <- rowSums(nes < 0)
  same_direction <- positive_n == length(cohorts) | negative_n == length(cohorts)
  fdr_n <- rowSums(qmat < 0.05, na.rm = TRUE)
  median_abs <- apply(abs(nes), 1L, median)
  signed_z <- sign(nes) * qnorm(pmax(pmat / 2, 1e-300), lower.tail = FALSE)
  combined_z <- rowSums(signed_z) / sqrt(length(cohorts))
  combined_p <- 2 * pnorm(-abs(combined_z))
  wide[, `:=`(
    collection = collection_name,
    positive_nes_cohort_n = positive_n,
    negative_nes_cohort_n = negative_n,
    all_three_same_nes_direction = same_direction,
    consensus_direction = fifelse(positive_n > negative_n, "POSITIVE_IN_DOR", fifelse(negative_n > positive_n, "NEGATIVE_IN_DOR", "TIE")),
    fgsea_fdr_lt_0_05_cohort_n = fdr_n,
    median_nes = apply(nes, 1L, median),
    median_abs_nes = median_abs,
    mean_nes = rowMeans(nes),
    sd_nes = apply(nes, 1L, sd),
    range_nes = apply(nes, 1L, function(x) max(x) - min(x)),
    signed_stouffer_z = combined_z,
    signed_stouffer_pvalue = combined_p,
    signed_stouffer_padj = p.adjust(combined_p, method = "BH"),
    strict_pathway_convergence = same_direction & fdr_n >= 2L & median_abs >= 1.0,
    directional_supported_pathway = same_direction & fdr_n >= 1L
  )]
  setcolorder(wide, c("collection", "pathway", "size", nes_cols, p_cols, padj_cols, setdiff(names(wide), c("collection", "pathway", "size", nes_cols, p_cols, padj_cols))))
  wide[]
}

convergence <- rbindlist(lapply(names(gene_sets), build_convergence))
setorder(convergence, collection, -strict_pathway_convergence, -fgsea_fdr_lt_0_05_cohort_n, signed_stouffer_padj, -median_abs_nes)
fwrite(convergence, file.path(results_dir, "M05_PATHWAY_CONVERGENCE.csv"))
fwrite(convergence[strict_pathway_convergence == TRUE], file.path(results_dir, "M05_STRICT_CONVERGENT_PATHWAYS.csv"))

m04_pair_path <- file.path(project_root, "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_PAIRWISE_REPRODUCIBILITY.csv")
m04_pairs <- fread(m04_pair_path)
pair_rows <- list()
pair_index <- 1L
for (collection_name in names(gene_sets)) {
  sub <- fgsea_results[collection == collection_name]
  for (i in seq_len(length(cohorts) - 1L)) {
    for (j in (i + 1L):length(cohorts)) {
      a <- cohorts[[i]]
      b <- cohorts[[j]]
      xa <- sub[cohort == a, .(pathway, NES_a = NES, padj_a = padj)]
      xb <- sub[cohort == b, .(pathway, NES_b = NES, padj_b = padj)]
      pair <- merge(xa, xb, by = "pathway", all = FALSE)
      cand_a <- pair$padj_a < 0.05
      cand_b <- pair$padj_b < 0.05
      overlap <- cand_a & cand_b
      union <- cand_a | cand_b
      gene_row <- m04_pairs[cohort_a == a & cohort_b == b]
      pair_rows[[pair_index]] <- data.table(
        collection = collection_name,
        cohort_a = a,
        cohort_b = b,
        pathway_n = nrow(pair),
        pathway_nes_spearman_rho = cor(pair$NES_a, pair$NES_b, method = "spearman"),
        pathway_nes_direction_concordance = mean(sign(pair$NES_a) == sign(pair$NES_b)),
        fdr_pathway_a_n = sum(cand_a),
        fdr_pathway_b_n = sum(cand_b),
        fdr_overlap_n = sum(overlap),
        fdr_jaccard = ifelse(sum(union) > 0L, sum(overlap) / sum(union), NA_real_),
        fdr_overlap_direction_concordance = ifelse(sum(overlap) > 0L, mean(sign(pair$NES_a[overlap]) == sign(pair$NES_b[overlap])), NA_real_),
        m04_gene_log2fc_spearman_rho = gene_row$spearman_rho_log2fc,
        pathway_minus_gene_rho = cor(pair$NES_a, pair$NES_b, method = "spearman") - gene_row$spearman_rho_log2fc
      )
      pair_index <- pair_index + 1L
    }
  }
}
pairwise <- rbindlist(pair_rows)
fwrite(pairwise, file.path(results_dir, "M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv"))

split_edges <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1L]]
}
strict <- convergence[strict_pathway_convergence == TRUE]
leading_rows <- list()
if (nrow(strict) > 0L) {
  for (idx in seq_len(nrow(strict))) {
    collection_name <- strict$collection[[idx]]
    pathway_name <- strict$pathway[[idx]]
    edges <- lapply(cohorts, function(cohort_name) {
      split_edges(fgsea_results[collection == collection_name & pathway == pathway_name & cohort == cohort_name, leadingEdge][[1L]])
    })
    names(edges) <- cohorts
    union_edges <- Reduce(union, edges)
    intersection_edges <- Reduce(intersect, edges)
    jaccard <- function(x, y) {
      denom <- length(union(x, y))
      if (denom == 0L) NA_real_ else length(intersect(x, y)) / denom
    }
    leading_rows[[length(leading_rows) + 1L]] <- data.table(
      collection = collection_name,
      pathway = pathway_name,
      leading_edge_GSE274832_n = length(edges[["GSE274832"]]),
      leading_edge_GSE193136_n = length(edges[["GSE193136"]]),
      leading_edge_GSE232306_n = length(edges[["GSE232306"]]),
      leading_edge_union_n = length(union_edges),
      leading_edge_three_way_intersection_n = length(intersection_edges),
      leading_edge_three_way_intersection = paste(sort(intersection_edges), collapse = ";"),
      jaccard_GSE274832_GSE193136 = jaccard(edges[["GSE274832"]], edges[["GSE193136"]]),
      jaccard_GSE274832_GSE232306 = jaccard(edges[["GSE274832"]], edges[["GSE232306"]]),
      jaccard_GSE193136_GSE232306 = jaccard(edges[["GSE193136"]], edges[["GSE232306"]])
    )
  }
  leading_overlap <- rbindlist(leading_rows)
} else {
  leading_overlap <- data.table(
    collection = character(), pathway = character(),
    leading_edge_GSE274832_n = integer(), leading_edge_GSE193136_n = integer(), leading_edge_GSE232306_n = integer(),
    leading_edge_union_n = integer(), leading_edge_three_way_intersection_n = integer(), leading_edge_three_way_intersection = character(),
    jaccard_GSE274832_GSE193136 = numeric(), jaccard_GSE274832_GSE232306 = numeric(), jaccard_GSE193136_GSE232306 = numeric()
  )
}
fwrite(leading_overlap, file.path(results_dir, "M05_STRICT_PATHWAY_LEADING_EDGE_OVERLAP.csv"))

key_counts <- rbindlist(list(
  geneset_audit[, .(metric = paste0("eligible_sets_", collection), value = eligible_sets_15_to_500)],
  convergence[, .(
    metric = c(
      paste0("strict_convergent_", unique(collection)),
      paste0("all_three_same_direction_", unique(collection)),
      paste0("directional_supported_", unique(collection)),
      paste0("signed_stouffer_fdr_", unique(collection))
    ),
    value = c(
      sum(strict_pathway_convergence),
      sum(all_three_same_nes_direction),
      sum(directional_supported_pathway),
      sum(signed_stouffer_padj < 0.05)
    )
  ), by = collection][, .(metric, value)],
  fgsea_results[, .(metric = paste0("fgsea_fdr_", collection, "_", cohort), value = sum(padj < 0.05)), by = .(collection, cohort)][, .(metric, value)]
))
fwrite(key_counts, file.path(results_dir, "M05_KEY_COUNTS.csv"))

clean_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  gsub("_", " ", x)
}

plot_heatmap <- function(selected, filename, title, add_values = FALSE, plot_width = 8.2, wrap_width = NA_integer_) {
  nes_cols <- paste0("NES_", cohorts)
  long <- melt(selected, id.vars = c("collection", "pathway"), measure.vars = nes_cols, variable.name = "cohort", value.name = "NES")
  long[, cohort := sub("^NES_", "", cohort)]
  ordering <- selected[order(median_nes), pathway]
  long[, pathway_label := factor(clean_pathway(pathway), levels = clean_pathway(ordering))]
  p <- ggplot(long, aes(x = factor(cohort, levels = cohorts, labels = cohort_short), y = pathway_label, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2563EB", mid = "white", high = "#DC2626", midpoint = 0) +
    labs(title = title, x = "Cohort", y = NULL, fill = "NES") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 7.5))
  if (!is.na(wrap_width)) {
    p <- p + scale_y_discrete(labels = function(x) vapply(as.character(x), function(label) paste(strwrap(label, width = wrap_width), collapse = "\n"), FUN.VALUE = character(1)))
  }
  if (add_values) p <- p + geom_text(aes(label = sprintf("%.2f", NES)), size = 2.4)
  height <- max(6, 0.22 * nrow(selected) + 2.5)
  ggsave(file.path(plots_dir, filename), p, width = plot_width, height = height, dpi = 300, bg = "white")
}

hallmark <- convergence[collection == "HALLMARK"]
plot_heatmap(hallmark, "M05_HALLMARK_NES_HEATMAP.png", "Hallmark pathway NES across three DOR cohorts", FALSE)

reactome <- convergence[collection == "REACTOME"]
setorder(reactome, -strict_pathway_convergence, -fgsea_fdr_lt_0_05_cohort_n, signed_stouffer_padj, -median_abs_nes)
reactome_top <- head(reactome, 30L)
plot_heatmap(
  reactome_top,
  "M05_REACTOME_TOP30_NES_HEATMAP.png",
  sprintf("Reactome top 30; strict pathways prioritized (%d strict total)", sum(reactome$strict_pathway_convergence)),
  TRUE,
  plot_width = 15,
  wrap_width = 52L
)

repro_long <- melt(
  pairwise,
  id.vars = c("collection", "cohort_a", "cohort_b"),
  measure.vars = c("m04_gene_log2fc_spearman_rho", "pathway_nes_spearman_rho"),
  variable.name = "level", value.name = "rho"
)
repro_long[, pair := paste(cohort_short[cohort_a], cohort_short[cohort_b], sep = " vs ")]
repro_long[, level := factor(level, levels = c("m04_gene_log2fc_spearman_rho", "pathway_nes_spearman_rho"), labels = c("M04 gene log2FC", "M05 pathway NES"))]
p_repro <- ggplot(repro_long, aes(x = pair, y = rho, fill = level)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, color = "#64748B", linewidth = 0.5) +
  facet_wrap(~collection) +
  scale_fill_manual(values = c("#94A3B8", "#7C3AED")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(title = "Pathway-level versus gene-level cross-cohort rank reproducibility", x = NULL, y = "Spearman rho", fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plots_dir, "M05_PATHWAY_VS_GENE_REPRODUCIBILITY.png"), p_repro, width = 11, height = 5.2, dpi = 300, bg = "white")

scatter_rows <- list()
for (collection_name in names(gene_sets)) {
  sub <- fgsea_results[collection == collection_name]
  for (i in seq_len(length(cohorts) - 1L)) {
    for (j in (i + 1L):length(cohorts)) {
      a <- cohorts[[i]]; b <- cohorts[[j]]
      xa <- sub[cohort == a, .(pathway, NES_a = NES)]
      xb <- sub[cohort == b, .(pathway, NES_b = NES)]
      joined <- merge(xa, xb, by = "pathway")
      joined[, `:=`(collection = collection_name, pair = paste(cohort_short[a], cohort_short[b], sep = " vs "))]
      scatter_rows[[length(scatter_rows) + 1L]] <- joined
    }
  }
}
scatter <- rbindlist(scatter_rows)
p_scatter <- ggplot(scatter, aes(NES_a, NES_b)) +
  geom_point(size = 0.7, alpha = 0.35, color = "#334155") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8") +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "#CBD5E1", linewidth = 0.3) +
  facet_grid(collection ~ pair, scales = "free") +
  labs(title = "Pairwise pathway NES reproducibility", x = "Cohort A NES", y = "Cohort B NES") +
  theme_minimal(base_size = 10)
ggsave(file.path(plots_dir, "M05_PAIRWISE_NES_SCATTER.png"), p_scatter, width = 13, height = 7.2, dpi = 300, bg = "white")

count_plot <- convergence[, .(
  `All-three same direction` = sum(all_three_same_nes_direction),
  `Same direction + >=1 cohort FDR` = sum(directional_supported_pathway),
  `Strict convergence` = sum(strict_pathway_convergence)
), by = collection]
count_long <- melt(count_plot, id.vars = "collection", variable.name = "category", value.name = "pathway_n")
p_counts <- ggplot(count_long, aes(x = category, y = pathway_n, fill = collection)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = pathway_n), position = position_dodge(width = 0.9), vjust = -0.2, size = 3) +
  scale_fill_manual(values = c(HALLMARK = "#0EA5E9", REACTOME = "#F59E0B")) +
  labs(title = "Pathway convergence counts", x = NULL, y = "Pathways", fill = "Collection") +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plots_dir, "M05_CONVERGENCE_COUNTS.png"), p_counts, width = 10.5, height = 5.2, dpi = 300, bg = "white")

p_disp <- ggplot(convergence, aes(x = median_abs_nes, y = sd_nes, color = strict_pathway_convergence)) +
  geom_point(alpha = 0.55, size = 1.4) +
  facet_wrap(~collection, scales = "free") +
  scale_color_manual(values = c(`FALSE` = "#94A3B8", `TRUE` = "#DC2626"), labels = c("Other", "Strict convergence")) +
  labs(title = "Pathway NES magnitude and cross-cohort dispersion", x = "Median |NES|", y = "SD of NES across three cohorts", color = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(plots_dir, "M05_PATHWAY_NES_DISPERSION.png"), p_disp, width = 10.5, height = 5.2, dpi = 300, bg = "white")

pair_lines <- vapply(seq_len(nrow(pairwise)), function(i) {
  row <- pairwise[i]
  sprintf(
    "- %s, %s vs %s: pathway NES rho=%.4f, direction concordance=%.2f%%, FDR overlap=%d, pathway-minus-gene rho=%+.4f.",
    row$collection, row$cohort_a, row$cohort_b, row$pathway_nes_spearman_rho,
    100 * row$pathway_nes_direction_concordance, row$fdr_overlap_n, row$pathway_minus_gene_rho
  )
}, character(1))

collection_summary <- convergence[, .(
  pathway_n = .N,
  all_three_same_direction_n = sum(all_three_same_nes_direction),
  directional_supported_n = sum(directional_supported_pathway),
  strict_n = sum(strict_pathway_convergence),
  signed_stouffer_fdr_n = sum(signed_stouffer_padj < 0.05)
), by = collection]
collection_lines <- vapply(seq_len(nrow(collection_summary)), function(i) {
  row <- collection_summary[i]
  sprintf(
    "- %s: %d tested; %d all-three same direction; %d same direction with >=1 cohort FDR; %d strict; %d signed-Stouffer FDR<0.05.",
    row$collection, row$pathway_n, row$all_three_same_direction_n, row$directional_supported_n, row$strict_n, row$signed_stouffer_fdr_n
  )
}, character(1))

strict_display <- strict[order(collection, signed_stouffer_padj, -median_abs_nes)]
if (nrow(strict_display) > 0L) {
  strict_lines <- vapply(seq_len(min(nrow(strict_display), 30L)), function(i) {
    row <- strict_display[i]
    sprintf(
      "- %s / %s: NES %.3f, %.3f, %.3f; cohort FDR count=%d; signed-Stouffer FDR=%.3g.",
      row$collection, row$pathway, row$NES_GSE274832, row$NES_GSE193136, row$NES_GSE232306,
      row$fgsea_fdr_lt_0_05_cohort_n, row$signed_stouffer_padj
    )
  }, character(1))
} else {
  strict_lines <- "- None."
}

summary_text <- paste0(
  "# M05 Pathway and Module Convergence Summary\n\n",
  "`PRELIMINARY_MODULE_STATE: PASS_WITH_LIMITATION_PENDING_INDEPENDENT_OUTPUT_AUDIT`\n\n",
  "## Scope and frozen resources\n\n",
  "Three cohort-specific DOR-minus-NOR signed Wald-z rank vectors were built from the full M04 common mapped universe. No candidate, p-value, meta-FDR, FMNL1, or PGAP1 prefilter was used. After outcome-independent symbol collapse, each rank contains ", format(nrow(rank_table), big.mark = ","), " identical unique symbols.\n\n",
  "The frozen resources are Human MSigDB v2026.1.Hs Hallmark (primary compact collection) and Reactome (secondary granular collection). After intersecting the mapped universe and applying the prespecified 15-500 size range, ", geneset_audit[collection == "HALLMARK", eligible_sets_15_to_500], " Hallmark and ", geneset_audit[collection == "REACTOME", eligible_sets_15_to_500], " Reactome pathways were tested.\n\n",
  "fgsea reported ", sum(fgsea_warning_audit$pvalue_at_eps_1e_10_n), " pathway-cohort results at the prespecified eps=1e-10 floor, with the same number lacking finite log2err. These p-values are treated as capped/upper-bounded for ranking and FDR; no exact smaller p-value is claimed. See M05_FGSEA_NUMERICAL_WARNING_AUDIT.csv.\n\n",
  "## Pathway-level reproducibility\n\n",
  paste(pair_lines, collapse = "\n"), "\n\n",
  "Pathway-level rho must be compared directly with the M04 gene-level rho in the output table; an increase for one pair or collection does not by itself establish universal convergence.\n\n",
  "## Convergence counts\n\n",
  paste(collection_lines, collapse = "\n"), "\n\n",
  "Strict convergence requires identical NES direction in all three cohorts, fgsea FDR<0.05 in at least two cohorts, and median |NES|>=1.0. Signed Stouffer FDR is supportive prioritization only and is not a replacement for this rule.\n\n",
  "## Strict convergent pathways (up to 30 listed)\n\n",
  paste(strict_lines, collapse = "\n"), "\n\n",
  "## Interpretation boundary\n\n",
  "- Positive NES means enrichment toward the DOR-up end of the within-cohort rank; negative NES means enrichment toward the DOR-down end. It does not prove pathway activation or inhibition.\n",
  "- Hallmark is the primary compact result. Reactome is secondary and contains biologically overlapping terms, so raw pathway counts are not independent discoveries.\n",
  "- GSE232306 has phenotype-aligned GC, duplication, and PC1 signals with unavailable individual age and incomplete batch/clinical covariates. M06 must quantify whether any main pathway conclusion depends on this cohort.\n",
  "- M05 does not validate biomarkers, causal mechanisms, treatments, or clinical prediction.\n\n",
  "## Permitted claims\n\n",
  "- The analysis quantified pathway NES direction, ranking, FDR overlap, and cross-cohort dispersion using frozen full-rank inputs.\n",
  "- Pathways satisfying the exact strict rule may be called exploratory cross-cohort convergent pathways.\n",
  "- Gene-level and pathway-level reproducibility may be contrasted quantitatively.\n\n",
  "## Prohibited claims\n\n",
  "- Do not call positive/negative NES pathway activation/inhibition without orthogonal functional evidence.\n",
  "- Do not call a signed-Stouffer-only pathway replicated.\n",
  "- Do not omit collection redundancy, cohort limitations, or the required M06 leave-one-cohort-out step.\n"
)
writeLines(summary_text, file.path(results_dir, "M05_PATHWAY_CONVERGENCE_SUMMARY.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(logs_dir, "R_SESSION_INFO.txt"))
cat(sprintf("M05 fgsea completed: %d rank symbols; %d fgsea rows; %d strict pathways\n", nrow(rank_table), nrow(fgsea_results), nrow(strict)))
print(pairwise)
print(collection_summary)
