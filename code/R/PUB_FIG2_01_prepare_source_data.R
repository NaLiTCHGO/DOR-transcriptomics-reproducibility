#!/usr/bin/env Rscript

# Prepare direct plotting tables for manuscript Figure 2 from locked M02/M03
# results. This script performs display-only transformations and validation; it
# does not recompute PCA coordinates, DESeq2 models, P values, or FDR values.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing required argument: ", flag)
  args[[idx + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), winslash = "/", mustWork = TRUE)
run_dir <- normalizePath(arg_value("--run-dir"), winslash = "/", mustWork = TRUE)
source_dir <- file.path(run_dir, "source_data")
upstream_dir <- file.path(source_dir, "upstream_refs")
qc_dir <- file.path(run_dir, "qc")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upstream_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

cohort_spec <- tibble::tribble(
  ~cohort, ~pca_rel, ~de_rel, ~pca_rows, ~de_rows, ~n_dor, ~n_nor, ~model, ~expected_candidates,
  "GSE274832",
  "06_locked_results/modules/M02_COHORT_QC/v1_GSE274832/expression_qc/GSE274832_PCA_COORDINATES.csv",
  "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
  6L, 16012L, 3L, 3L, "condition-only", 30L,
  "GSE193136",
  "06_locked_results/modules/M02_COHORT_QC/v1_GSE193136/expression_qc/GSE193136_PCA_COORDINATES.csv",
  "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_DESEQ2_AGE_ADJUSTED_DOR_MINUS_NOR_ALL_GENES.csv",
  12L, 16656L, 6L, 6L, "age-adjusted primary", 22L,
  "GSE232306",
  "06_locked_results/modules/M02_COHORT_QC/v1_GSE232306/expression_qc/GSE232306_PCA_COORDINATES.csv",
  "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
  12L, 16064L, 6L, 6L, "condition-only; covariates unavailable", 6381L
)

required_pca <- c("srr", "phenotype", "PC1", "PC2", "PC1_variance_percent", "PC2_variance_percent")
required_de <- c("gene_id_version", "gene_id", "gene_name", "log2FoldChange", "pvalue", "padj")

pca_parts <- list()
de_parts <- list()
summary_parts <- list()
upstream_parts <- list()

for (i in seq_len(nrow(cohort_spec))) {
  spec <- cohort_spec[i, ]
  pca_path <- file.path(project_root, spec$pca_rel)
  de_path <- file.path(project_root, spec$de_rel)
  if (!file.exists(pca_path)) stop("Missing locked PCA table: ", pca_path)
  if (!file.exists(de_path)) stop("Missing locked DE table: ", de_path)

  pca <- read_csv(pca_path, na = c("", "NA"), show_col_types = FALSE)
  de <- read_csv(de_path, na = c("", "NA"), show_col_types = FALSE)
  missing_pca <- setdiff(required_pca, names(pca))
  missing_de <- setdiff(required_de, names(de))
  if (length(missing_pca)) stop(spec$cohort, " PCA missing columns: ", paste(missing_pca, collapse = ", "))
  if (length(missing_de)) stop(spec$cohort, " DE missing columns: ", paste(missing_de, collapse = ", "))
  if (nrow(pca) != spec$pca_rows) stop(spec$cohort, " unexpected PCA rows: ", nrow(pca))
  if (nrow(de) != spec$de_rows) stop(spec$cohort, " unexpected DE rows: ", nrow(de))
  if (anyDuplicated(pca$srr)) stop(spec$cohort, " PCA sample IDs are not unique")
  if (sum(pca$phenotype == "DOR") != spec$n_dor || sum(pca$phenotype == "NOR") != spec$n_nor) {
    stop(spec$cohort, " phenotype counts do not match the frozen design")
  }
  if (any(!is.finite(pca$PC1)) || any(!is.finite(pca$PC2))) stop(spec$cohort, " PCA contains non-finite coordinates")
  if (length(unique(pca$PC1_variance_percent)) != 1L || length(unique(pca$PC2_variance_percent)) != 1L) {
    stop(spec$cohort, " PCA variance percentages are not constant within cohort")
  }

  age <- if ("age_years" %in% names(pca)) as.numeric(pca$age_years) else rep(NA_real_, nrow(pca))
  pca_out <- tibble(
    cohort = spec$cohort,
    sample_id = pca$srr,
    gsm = if ("gsm" %in% names(pca)) pca$gsm else NA_character_,
    phenotype = factor(pca$phenotype, levels = c("NOR", "DOR")),
    age_years = age,
    PC1 = as.numeric(pca$PC1),
    PC2 = as.numeric(pca$PC2),
    PC1_variance_percent = as.numeric(pca$PC1_variance_percent),
    PC2_variance_percent = as.numeric(pca$PC2_variance_percent)
  )

  lfc <- as.numeric(de$log2FoldChange)
  pvalue <- as.numeric(de$pvalue)
  padj <- as.numeric(de$padj)
  eligible <- is.finite(lfc) & is.finite(pvalue) & pvalue >= 0 & pvalue <= 1
  if (!any(eligible)) stop(spec$cohort, " has no eligible volcano rows")
  positive_p <- pvalue[eligible & pvalue > 0]
  p_floor <- if (length(positive_p)) max(min(positive_p) / 10, 1e-300) else 1e-300
  neg_log10_p <- rep(NA_real_, length(pvalue))
  neg_log10_p[eligible] <- -log10(pmax(pvalue[eligible], p_floor))
  candidate <- eligible & !is.na(padj) & padj < 0.05 & abs(lfc) >= 1
  if (sum(candidate) != spec$expected_candidates) {
    stop(spec$cohort, " candidate count changed: expected ", spec$expected_candidates, ", observed ", sum(candidate))
  }

  y <- neg_log10_p[eligible]
  ax <- abs(lfc[eligible])
  y_q99 <- unname(quantile(y, 0.99, na.rm = TRUE, type = 7))
  y_max <- max(y, na.rm = TRUE)
  x_q99 <- unname(quantile(ax, 0.99, na.rm = TRUE, type = 7))
  x_max <- max(ax, na.rm = TRUE)
  extreme_y <- is.finite(y_q99) && y_max > max(60, y_q99 * 3)
  extreme_x <- is.finite(x_q99) && x_max > max(8, x_q99 * 3)
  suggested_y_cap <- if (extreme_y) max(12, y_q99 * 1.15) else NA_real_
  route <- if (extreme_x) {
    "EXTREME_X_FULL_RANGE_REVIEW"
  } else if (extreme_y) {
    "EXTREME_Y_CAPPED_NO_LABEL_REVIEW"
  } else {
    "STANDARD_NO_LABEL"
  }

  gene_label <- ifelse(is.na(de$gene_name) | de$gene_name == "", de$gene_id, de$gene_name)
  de_out <- tibble(
    cohort = spec$cohort,
    gene_id_version = de$gene_id_version,
    gene_id = de$gene_id,
    gene_name = gene_label,
    log2FoldChange = lfc,
    pvalue = pvalue,
    padj = padj,
    neg_log10_nominal_p = neg_log10_p,
    plot_eligible = eligible,
    candidate = candidate,
    display_status = factor(ifelse(candidate, "Meets FDR and effect threshold", "Other"),
                            levels = c("Other", "Meets FDR and effect threshold")),
    model = spec$model,
    display_route = route,
    suggested_y_cap = suggested_y_cap
  )

  pca_parts[[spec$cohort]] <- pca_out
  de_parts[[spec$cohort]] <- de_out
  summary_parts[[spec$cohort]] <- tibble(
    cohort = spec$cohort,
    model = spec$model,
    n_samples = nrow(pca_out),
    n_dor = sum(pca_out$phenotype == "DOR"),
    n_nor = sum(pca_out$phenotype == "NOR"),
    pc1_variance_percent = unique(pca_out$PC1_variance_percent),
    pc2_variance_percent = unique(pca_out$PC2_variance_percent),
    de_input_rows = nrow(de_out),
    de_plot_rows = sum(de_out$plot_eligible),
    de_excluded_rows = sum(!de_out$plot_eligible),
    candidate_count = sum(de_out$candidate),
    zero_nominal_p_rows = sum(de_out$pvalue == 0, na.rm = TRUE),
    nominal_p_display_floor = p_floor,
    max_neg_log10_nominal_p = y_max,
    max_abs_log2_fc = x_max,
    label_mode = "none",
    display_route = route,
    suggested_y_cap = suggested_y_cap
  )
  upstream_parts[[paste0(spec$cohort, "_PCA")]] <- tibble(
    panel_source = paste0(spec$cohort, "_PCA"), locked_file = pca_path,
    role = "Locked PCA coordinates; no PCA recomputation"
  )
  upstream_parts[[paste0(spec$cohort, "_DE")]] <- tibble(
    panel_source = paste0(spec$cohort, "_DE"), locked_file = de_path,
    role = "Locked DESeq2 all-gene table; display-only transformation"
  )
}

pca_all <- bind_rows(pca_parts)
de_all <- bind_rows(de_parts)
panel_summary <- bind_rows(summary_parts)
upstream <- bind_rows(upstream_parts)

if (nrow(pca_all) != 30L) stop("Combined PCA source must contain 30 rows")
if (nrow(de_all) != 48732L) stop("Combined DE source must contain 48,732 rows")

write_csv(pca_all, file.path(source_dir, "FIG2_PCA_DIRECT_SOURCE.csv"), na = "")
write_csv(de_all, file.path(source_dir, "FIG2_VOLCANO_DIRECT_SOURCE.csv"), na = "")
write_csv(panel_summary, file.path(source_dir, "FIG2_PANEL_SUMMARY.csv"), na = "")
write_csv(upstream, file.path(upstream_dir, "UPSTREAM_LOCKED_RESULT_REFERENCES.csv"), na = "")

report <- c(
  "# Figure 2 source-preparation validation",
  "",
  "- Status: `PASS`",
  "- Scientific recomputation: `NO`",
  "- PCA coordinates: read unchanged from locked M02 tables.",
  "- DE statistics: read unchanged from locked M03 all-gene tables.",
  "- Display-only transforms: eligibility flag, -log10 nominal P, candidate colour class, and envelope diagnostics.",
  sprintf("- Combined PCA rows: %s.", format(nrow(pca_all), big.mark = ",")),
  sprintf("- Combined DE rows: %s.", format(nrow(de_all), big.mark = ",")),
  sprintf("- Candidate counts: %s.", paste(paste0(panel_summary$cohort, "=", panel_summary$candidate_count), collapse = "; ")),
  "",
  "The direct plotting tables are authoritative for visual revisions of Figure 2; later layout changes must not rerun M02 or M03."
)
writeLines(report, file.path(qc_dir, "SOURCE_PREPARATION_VALIDATION.md"), useBytes = TRUE)
cat("PASS: Figure 2 direct plotting sources prepared from locked results\n")

