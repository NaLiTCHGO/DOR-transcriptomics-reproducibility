#!/usr/bin/env Rscript

# Automated structural and scientific-contract QC for the Figure 2 v07 QC
# candidate. Human actual-render inspection remains a separate mandatory gate.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(magick)
})

arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing required argument: ", flag)
  args[[idx + 1L]]
}

run_dir <- normalizePath(arg_value("--run-dir"), winslash = "/", mustWork = TRUE)
source_dir <- file.path(run_dir, "source_data")
figure_dir <- file.path(run_dir, "figures")
qc_dir <- file.path(run_dir, "qc")
candidate <- file.path(figure_dir, "Figure_2_Cohort_QC_and_Effects_v07_QC.png")

checks <- list()
add_check <- function(check, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = check,
    result = if (isTRUE(pass)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

add_check("candidate_exists", file.exists(candidate), file.exists(candidate), TRUE)
if (!file.exists(candidate)) stop("Candidate PNG is missing")
info <- image_info(image_read(candidate))
add_check("pixel_width", info$width == 2008L, info$width, 2008L)
add_check("pixel_height", info$height == 2409L, info$height, 2409L)
add_check("colorspace", info$colorspace %in% c("sRGB", "RGB"), info$colorspace, "sRGB/RGB")
add_check("format", info$format == "PNG", info$format, "PNG")
add_check("nonempty_file", file.info(candidate)$size > 100000L, file.info(candidate)$size, ">100000 bytes")

pca <- read_csv(file.path(source_dir, "FIG2_PCA_DIRECT_SOURCE.csv"), show_col_types = FALSE)
volcano <- read_csv(file.path(source_dir, "FIG2_VOLCANO_DIRECT_SOURCE.csv"), show_col_types = FALSE)
summary <- read_csv(file.path(source_dir, "FIG2_PANEL_SUMMARY.csv"), show_col_types = FALSE)
warnings <- read_csv(file.path(qc_dir, "RUNTIME_WARNING_LOG.csv"), show_col_types = FALSE)

add_check("pca_rows", nrow(pca) == 30L, nrow(pca), 30L)
add_check("volcano_rows", nrow(volcano) == 48732L, nrow(volcano), 48732L)
add_check("cohort_count", n_distinct(pca$cohort) == 3L && n_distinct(volcano$cohort) == 3L,
          paste(n_distinct(pca$cohort), n_distinct(volcano$cohort), sep = "/"), "3/3")
expected_candidates <- c(GSE274832 = 30L, GSE193136 = 22L, GSE232306 = 6381L)
observed_candidates <- setNames(summary$candidate_count, summary$cohort)[names(expected_candidates)]
add_check("candidate_counts", identical(as.integer(observed_candidates), as.integer(expected_candidates)),
          paste(observed_candidates, collapse = "/"), paste(expected_candidates, collapse = "/"))
add_check("runtime_warnings", nrow(warnings) == 0L, nrow(warnings), 0L)
add_check("candidate_first_only_one_png", length(list.files(figure_dir, pattern = "\\.png$", full.names = TRUE)) == 1L,
          length(list.files(figure_dir, pattern = "\\.png$", full.names = TRUE)), 1L)
add_check("no_premature_final_formats",
          length(list.files(figure_dir, pattern = "\\.(tif|tiff|pdf|svg)$", ignore.case = TRUE)) == 0L,
          length(list.files(figure_dir, pattern = "\\.(tif|tiff|pdf|svg)$", ignore.case = TRUE)), 0L)

qc <- bind_rows(checks)
write_csv(qc, file.path(qc_dir, "AUTOMATED_QC_SUMMARY.csv"), na = "")
overall <- if (all(qc$result == "PASS")) "PASS" else "FAIL"
report <- c(
  "# Figure 2 v07 automated QC",
  "",
  paste0("- Automated gate: `", overall, "`"),
  "- Human actual-render gate: `PENDING`",
  "- Candidate status: `B_CANDIDATE_DRAFT`",
  "- Scientific/statistical values were not recomputed.",
  "",
  "| Check | Result | Observed | Expected |",
  "|---|---|---|---|",
  apply(qc, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")),
  "",
  "Automated PASS does not authorize submission formats or immutable release."
)
writeLines(report, file.path(qc_dir, "AUTOMATED_QC_REPORT.md"), useBytes = TRUE)
if (overall != "PASS") stop("Figure 2 automated QC failed")
cat("PASS: Figure 2 v07 automated QC; human visual QC remains pending\n")
