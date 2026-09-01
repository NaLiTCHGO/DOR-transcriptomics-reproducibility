options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(magick)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) stop("Missing --", key, "=<path>")
  sub(paste0("^--", key, "="), "", hit[[1]])
}
run_root <- normalizePath(get_arg("run-root"), winslash = "/", mustWork = TRUE)
figure_root <- file.path(run_root, "figures")
source_root <- file.path(run_root, "source_data")
qc_root <- file.path(run_root, "qc")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

expected <- data.frame(
  file = c(
    "main/Figure_3_Gene_Reproducibility_and_Heterogeneity_v09_QC.png",
    "main/Figure_4_Pathway_Convergence_v05_QC.png",
    "main/Figure_5_LOCO_Robustness_v08_QC.png",
    "supplement/Figure_S1_GSE274832_QC_and_Stability_v08_QC.png",
    "supplement/Figure_S2_GSE193136_QC_and_Stability_v08_QC.png",
    "supplement/Figure_S3_GSE232306_QC_and_Stability_v08_QC.png",
    "supplement/Figure_S4_Gene_Synthesis_Sensitivity_v05_QC.png",
    "supplement/Figure_S5_Pathway_Context_v07_QC.png",
    "supplement/Figure_S6_LOCO_Context_v04_QC.png"
  ),
  width = rep(2008L, 9),
  height = c(2409L, 2600L, 2200L, 2600L, 2600L, 2600L, 1150L, 2300L, 1150L),
  stringsAsFactors = FALSE
)
expected$path <- file.path(figure_root, expected$file)
if (any(!file.exists(expected$path))) stop("Missing candidates: ", paste(expected$file[!file.exists(expected$path)], collapse = "; "))

rows <- lapply(seq_len(nrow(expected)), function(i) {
  info <- image_info(image_read(expected$path[[i]]))
  data.frame(
    figure_file = expected$file[[i]],
    width_px = info$width[[1]], height_px = info$height[[1]],
    format = info$format[[1]], colorspace = info$colorspace[[1]],
    bytes = file.info(expected$path[[i]])$size,
    dimension_status = ifelse(info$width[[1]] == expected$width[[i]] && info$height[[1]] == expected$height[[i]], "PASS", "FAIL"),
    format_status = ifelse(toupper(info$format[[1]]) == "PNG", "PASS", "FAIL"),
    size_status = ifelse(file.info(expected$path[[i]])$size >= 50000, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
})
image_qc <- bind_rows(rows)

all_candidates <- list.files(figure_root, pattern = "_QC\\.png$", recursive = TRUE, full.names = TRUE)
forbidden <- list.files(run_root, pattern = "(600|\\.tiff?$|\\.pdf$)", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
warning_files <- file.path(qc_root, c(
  "RUNTIME_WARNING_LOG_VST_RECONSTRUCTION.csv",
  "RUNTIME_WARNING_LOG_MAIN.csv",
  "RUNTIME_WARNING_LOG_SUPP_COHORTS.csv",
  "RUNTIME_WARNING_LOG_SUPP_SYNTHESIS.csv"
))
if (any(!file.exists(warning_files))) stop("Missing runtime warning logs")
warning_counts <- vapply(warning_files, function(f) nrow(read_csv(f, show_col_types = FALSE)), integer(1))

meta <- read_csv(file.path(source_root, "FIG3_S4/M04_GENE_LEVEL_META_ANALYSIS.csv"), show_col_types = FALSE, progress = FALSE)
pairwise <- read_csv(file.path(source_root, "FIG3_S4/M04_PAIRWISE_REPRODUCIBILITY.csv"), show_col_types = FALSE)
conv <- read_csv(file.path(source_root, "FIG4_S5/M05_PATHWAY_CONVERGENCE.csv"), show_col_types = FALSE, progress = FALSE)
strict <- read_csv(file.path(source_root, "FIG4_S5/M05_STRICT_CONVERGENT_PATHWAYS.csv"), show_col_types = FALSE)
summary <- read_csv(file.path(source_root, "FIG5_S6/M06_COLLECTION_LOCO_SUMMARY.csv"), show_col_types = FALSE)
core <- read_csv(file.path(source_root, "FIG5_S6/M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv"), show_col_types = FALSE)
frozen <- read_csv(file.path(source_root, "FIG5_S6/M06_FROZEN_M05_STRICT_SURVIVAL.csv"), show_col_types = FALSE)

contracts <- data.frame(
  check = c(
    "Exactly nine batch candidates", "No premature final-release formats", "Zero runtime warnings",
    "VST reconstructions validated",
    "M04 shared genes", "M04 pairwise comparisons", "M04 high-I2 genes", "M04 strict exploratory genes",
    "M05 total pathways", "M05 strict Hallmark", "M05 strict Reactome",
    "M06 rotations", "M06 universal core", "M06 Hallmark core", "M06 Reactome core",
    "M06 frozen Hallmark rows", "All image dimensions", "All image formats", "All candidate files nontrivial"
  ),
  observed = c(
    length(all_candidates), length(forbidden), sum(warning_counts),
    nrow(read_csv(file.path(qc_root, "VST_HEATMAP_RECONSTRUCTION_VALIDATION.csv"), show_col_types = FALSE)),
    nrow(meta), nrow(pairwise), sum(meta$i2_percent >= 75), sum(meta$strict_exploratory_consensus),
    nrow(conv), sum(strict$collection == "HALLMARK"), sum(strict$collection == "REACTOME"),
    nrow(summary), sum(core$universal_loco_core), sum(core$universal_loco_core & core$collection == "HALLMARK"),
    sum(core$universal_loco_core & core$collection == "REACTOME"), sum(frozen$collection == "HALLMARK"),
    sum(image_qc$dimension_status == "PASS"), sum(image_qc$format_status == "PASS"), sum(image_qc$size_status == "PASS")
  ),
  expected = c(9, 0, 0, 3, 13993, 3, 5883, 2, 1065, 8, 47, 6, 8, 1, 7, 8, 9, 9, 9),
  stringsAsFactors = FALSE
)
contracts$status <- ifelse(as.numeric(contracts$observed) == as.numeric(contracts$expected), "PASS", "FAIL")

write_csv(image_qc, file.path(qc_root, "AUTOMATED_IMAGE_QC.csv"))
write_csv(contracts, file.path(qc_root, "AUTOMATED_QC_SUMMARY.csv"))
gate <- all(contracts$status == "PASS") &&
  all(image_qc$dimension_status == "PASS") && all(image_qc$format_status == "PASS") && all(image_qc$size_status == "PASS")
report <- c(
  "# Batch Automated QC Report",
  "",
  "- Scope: Figure 3-5 and Supplementary Figures S1-S6",
  paste0("- Expected candidates found: ", length(all_candidates), "/9"),
  paste0("- Contract checks: ", sum(contracts$status == "PASS"), "/", nrow(contracts), " PASS"),
  paste0("- Image structural checks: ", sum(image_qc$dimension_status == "PASS" & image_qc$format_status == "PASS" & image_qc$size_status == "PASS"), "/9 PASS"),
  paste0("- Runtime warnings: ", sum(warning_counts)),
  paste0("- Premature TIFF/PDF/600-ppi artifacts: ", length(forbidden)),
  "- Scientific claim/result recomputation: NO; S1-S3 panel-D display matrices were reconstructed and identity-checked against locked normalized counts.",
  "- Current candidate class: `B_CANDIDATE_DRAFT_UNIFIED_VISUAL_QC_PENDING`",
  paste0("- Final Gate: ", if (gate) "PASS" else "HOLD")
)
writeLines(report, file.path(qc_root, "AUTOMATED_QC_REPORT.md"), useBytes = TRUE)
if (!gate) stop("Batch automated QC failed; see ", file.path(qc_root, "AUTOMATED_QC_REPORT.md"))
cat("PASS: batch automated QC; 9/9 candidates and 19/19 contracts PASS\n")
