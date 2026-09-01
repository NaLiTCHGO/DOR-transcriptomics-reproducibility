options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(key) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) stop("Missing --", key, "=<path>")
  sub(paste0("^--", key, "="), "", hit[[1]])
}
project_root <- normalizePath(parse_arg("project-root"), winslash = "/", mustWork = TRUE)
run_root_raw <- parse_arg("run-root")
dir.create(run_root_raw, recursive = TRUE, showWarnings = FALSE)
run_root <- normalizePath(run_root_raw, winslash = "/", mustWork = TRUE)
source_root <- file.path(run_root, "source_data")
qc_root <- file.path(run_root, "qc")
dir.create(source_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

spec <- data.frame(
  figure = c(
    rep("FIG3_S4", 5), rep("FIG4_S5", 4), rep("FIG5_S6", 4),
    rep("S1", 4), rep("S2", 4), rep("S3", 4)
  ),
  source_rel = c(
    "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_GENE_LEVEL_META_ANALYSIS.csv",
    "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_PAIRWISE_REPRODUCIBILITY.csv",
    "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_CANDIDATE_OVERLAP_SUMMARY.csv",
    "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_LOCO_STABILITY.csv",
    "06_locked_results/modules/M04_REPRO_HETEROGENEITY/v1/results/M04_KEY_COUNTS.csv",
    "06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_PATHWAY_CONVERGENCE.csv",
    "06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv",
    "06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_STRICT_CONVERGENT_PATHWAYS.csv",
    "06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results/M05_FGSEA_ALL_RESULTS.csv",
    "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_COLLECTION_LOCO_SUMMARY.csv",
    "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv",
    "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_FROZEN_M05_STRICT_SURVIVAL.csv",
    "06_locked_results/modules/M06_LEAVE_ONE_COHORT_OUT/v1/results/M06_ROBUSTNESS_TIER_COUNTS.csv",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE274832/expression_qc/GSE274832_MARKER_EXPRESSION.csv",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE274832/expression_qc/GSE274832_SAMPLE_CORRELATION.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_TOP30_NOMINAL_GENES_HEATMAP.png",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE193136/expression_qc/GSE193136_MARKER_EXPRESSION.csv",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE193136/expression_qc/GSE193136_SAMPLE_CORRELATION.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_TOP30_AGE_ADJUSTED_NOMINAL_GENES_HEATMAP.png",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE232306/expression_qc/GSE232306_MARKER_EXPRESSION.csv",
    "06_locked_results/modules/M02_COHORT_QC/v1_GSE232306/expression_qc/GSE232306_SAMPLE_CORRELATION.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv",
    "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_TOP30_NOMINAL_GENES_HEATMAP.png"
  ),
  role = c(
    "gene-level effects and heterogeneity", "pairwise gene reproducibility", "candidate overlap",
    "gene-level LOCO stability", "M04 key counts", "pathway convergence", "pathway pairwise reproducibility",
    "strict pathways", "per-cohort fgsea NES", "LOCO collection summary", "universal LOCO core",
    "frozen strict-pathway survival", "robustness tiers",
    rep(c("marker-expression table", "sample-correlation matrix", "LOSO global stability", "locked Top-30 VST heatmap raster"), 3)
  ),
  stringsAsFactors = FALSE
)

spec$source_abs <- file.path(project_root, spec$source_rel)
if (any(!file.exists(spec$source_abs))) {
  stop("Missing locked sources: ", paste(spec$source_rel[!file.exists(spec$source_abs)], collapse = "; "))
}

spec$destination_rel <- file.path(spec$figure, basename(spec$source_abs))
spec$destination_abs <- file.path(source_root, spec$destination_rel)
for (i in seq_len(nrow(spec))) {
  dir.create(dirname(spec$destination_abs[[i]]), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(spec$source_abs[[i]], spec$destination_abs[[i]], overwrite = TRUE, copy.date = TRUE)) {
    stop("Copy failed: ", spec$source_abs[[i]])
  }
}
spec$bytes <- file.info(spec$destination_abs)$size
spec$row_count <- NA_integer_
csv_rows <- grepl("\\.csv$", spec$destination_abs, ignore.case = TRUE)
spec$row_count[csv_rows] <- vapply(spec$destination_abs[csv_rows], function(f) nrow(read_csv(f, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")), integer(1))

read_src <- function(group, name) read_csv(file.path(source_root, group, name), show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
meta <- read_src("FIG3_S4", "M04_GENE_LEVEL_META_ANALYSIS.csv")
pairwise <- read_src("FIG3_S4", "M04_PAIRWISE_REPRODUCIBILITY.csv")
conv <- read_src("FIG4_S5", "M05_PATHWAY_CONVERGENCE.csv")
strict <- read_src("FIG4_S5", "M05_STRICT_CONVERGENT_PATHWAYS.csv")
loco <- read_src("FIG5_S6", "M06_COLLECTION_LOCO_SUMMARY.csv")
core <- read_src("FIG5_S6", "M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv")
frozen <- read_src("FIG5_S6", "M06_FROZEN_M05_STRICT_SURVIVAL.csv")

checks <- data.frame(
  check = c(
    "M04 common genes", "M04 pairwise rows", "M05 pathway rows", "M05 strict Hallmark", "M05 strict Reactome",
    "M06 collection rotations", "M06 universal core", "M06 Hallmark core", "M06 Reactome core",
    "S1 marker/LOSO", "S2 marker/LOSO", "S3 marker/LOSO"
  ),
  observed = c(
    nrow(meta), nrow(pairwise), nrow(conv), sum(strict$collection == "HALLMARK"), sum(strict$collection == "REACTOME"),
    nrow(loco), sum(core$universal_loco_core), sum(core$universal_loco_core & core$collection == "HALLMARK"),
    sum(core$universal_loco_core & core$collection == "REACTOME"),
    paste(nrow(read_src("S1", "GSE274832_MARKER_EXPRESSION.csv")), nrow(read_src("S1", "GSE274832_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv")), sep = "/"),
    paste(nrow(read_src("S2", "GSE193136_MARKER_EXPRESSION.csv")), nrow(read_src("S2", "GSE193136_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv")), sep = "/"),
    paste(nrow(read_src("S3", "GSE232306_MARKER_EXPRESSION.csv")), nrow(read_src("S3", "GSE232306_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv")), sep = "/")
  ),
  expected = c("13993", "3", "1065", "8", "47", "6", "8", "1", "7", "126/6", "360/12", "360/12"),
  stringsAsFactors = FALSE
)
checks$status <- ifelse(as.character(checks$observed) == checks$expected, "PASS", "FAIL")
if (any(checks$status != "PASS")) stop("Source contract failure: ", paste(checks$check[checks$status != "PASS"], collapse = "; "))

write_csv(spec[, c("figure", "source_rel", "destination_rel", "role", "bytes", "row_count")], file.path(source_root, "SOURCE_FILE_REGISTER.csv"))
write_csv(checks, file.path(qc_root, "SOURCE_CONTRACT_CHECKS.csv"))
report <- c(
  "# Batch Source Preparation Validation",
  "",
  "- Figures: Figure 3-5 and Supplementary Figures S1-S6",
  "- Upstream authority: locked M02-M06 results only",
  "- Scientific/statistical recomputation: NO",
  "- Daily hash gate: omitted under the user-approved local lightweight policy",
  "- S1-S3 Top-30 heatmaps: retained as locked 240-dpi raster panels because the VST plotting matrices were not saved; no VST was recomputed",
  paste0("- Copied source files: ", nrow(spec)),
  paste0("- Contract checks: ", sum(checks$status == "PASS"), "/", nrow(checks), " PASS"),
  "- Final Gate: PASS"
)
writeLines(report, file.path(qc_root, "SOURCE_PREPARATION_VALIDATION.md"), useBytes = TRUE)
cat("PASS: batch source package prepared; files=", nrow(spec), "; checks=", nrow(checks), "\n", sep = "")
