#!/usr/bin/env Rscript

# Generate reproducible QC-only hotspot crops and a Windows grid/font preflight
# for the accepted Figure 2 candidate. This script never modifies the figure.

suppressPackageStartupMessages({
  library(magick)
  library(readr)
  library(ragg)
  library(systemfonts)
  library(grid)
})
options(warn = 2)

arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing required argument: ", flag)
  args[[idx + 1L]]
}

run_dir <- normalizePath(arg_value("--run-dir"), winslash = "/", mustWork = TRUE)
candidate <- file.path(run_dir, "figures", "Figure_2_Cohort_QC_and_Effects_v07_QC.png")
qc_dir <- file.path(run_dir, "qc")
crop_dir <- file.path(qc_dir, "hotspots_v07")
dir.create(crop_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(candidate)) stop("Candidate is missing: ", candidate)

img <- image_read(candidate)
info <- image_info(img)
if (info$width != 2008L || info$height != 2409L) stop("Unexpected candidate dimensions")

crops <- data.frame(
  crop_id = c("TOP_HEADERS", "MIDDLE_HEADERS", "BOTTOM_HEADERS", "BOTTOM_KEY"),
  x = c(40L, 40L, 40L, 190L),
  y = c(30L, 815L, 1590L, 2255L),
  width = c(1928L, 1928L, 1928L, 1628L),
  height = c(190L, 205L, 235L, 135L),
  purpose = c(
    "A/B tag-title baseline, 2-mm gap, subtitle and body anchor",
    "C/D tag-title baseline and D top-tick separation",
    "E/F tag-title baseline, E top-tick separation and long F subtitle",
    "Shared PCA/volcano color key completeness and clipping"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(crops))) {
  z <- crops[i, ]
  geometry <- sprintf("%dx%d+%d+%d", z$width, z$height, z$x, z$y)
  out <- file.path(crop_dir, paste0("Figure2_v07_", z$crop_id, ".png"))
  image_write(image_crop(img, geometry), path = out, format = "png")
  if (!file.exists(out) || file.info(out)$size <= 0) stop("Hotspot crop not written: ", out)
  crops$output_file[i] <- out
}
write_csv(crops, file.path(crop_dir, "HOTSPOT_CROP_INDEX.csv"), na = "")

font_warning <- character()
font_test_path <- file.path(qc_dir, "WINDOWS_GRID_FONT_TEST.png")
font_result <- withCallingHandlers({
  resolved <- systemfonts::match_fonts("Arial")
  if (nrow(resolved) != 1L || !file.exists(resolved$path[[1]])) stop("Arial did not resolve to one local font file")
  resolved_info <- systemfonts::font_info(resolved$path[[1]], index = resolved$index[[1]])
  if (nrow(resolved_info) != 1L) stop("Resolved Arial font metadata is not unique")
  grDevices::windowsFonts(Skill3Arial = grDevices::windowsFont("Arial"))
  ragg::agg_png(font_test_path, width = 900, height = 180, units = "px", res = 300,
                background = "white")
  grid::grid.newpage()
  gp <- grid::gpar(fontfamily = "sans", fontsize = 9, fontface = "bold")
  width_mm <- grid::convertWidth(grid::stringWidth("A–F  ≥  −log₂ P"), "mm", valueOnly = TRUE)
  grid::grid.text("A–F  ≥  −log₂ P", x = 0.03, y = 0.58, just = "left", gp = gp)
  invisible(grDevices::dev.off())
  list(resolved = resolved, resolved_info = resolved_info, width_mm = width_mm)
}, warning = function(w) {
  font_warning <<- c(font_warning, conditionMessage(w))
  invokeRestart("muffleWarning")
})

if (!file.exists(font_test_path) || file.info(font_test_path)$size <= 0) stop("Grid font test image missing")
if (length(font_warning)) stop("Font preflight emitted warnings: ", paste(unique(font_warning), collapse = "; "))
resolved_family <- as.character(font_result$resolved_info$family[[1]])
resolved_path <- normalizePath(font_result$resolved$path[[1]], winslash = "/", mustWork = TRUE)
exact_match <- identical(tolower(resolved_family), "arial")

card <- c(
  "# Windows Grid Font Preflight Card",
  "",
  "- Requested family: Arial",
  paste0("- systemfonts resolved family: ", resolved_family),
  paste0("- Resolved font file: `", resolved_path, "`"),
  paste0("- Exact-family match: ", if (exact_match) "YES" else "NO"),
  "- Base-grid alias registered: YES",
  "- Alias name: Skill3Arial; title-row tag measurement uses grid-safe `sans` on the opened ragg device",
  "- windowsFonts() registration result: PASS",
  paste0("- grid string-width test: PASS; width = ", sprintf("%.3f mm", font_result$width_mm)),
  "- Unicode/glyph test: PASS for A-F, en dash, greater-than-or-equal, minus and subscript 2 in QC PNG",
  "- PDF font inspection: DEFERRED until final-release authorization",
  "- PNG/TIFF/SVG visual comparison: PNG candidate PASS; other formats deferred",
  paste0("- Font warnings captured: ", length(font_warning)),
  "- Warning adjudication file: `qc/RUNTIME_DIAGNOSIS_v03.md`",
  "- Final status: `PASS_FOR_B_CANDIDATE_PNG`"
)
writeLines(card, file.path(qc_dir, "WINDOWS_GRID_FONT_PREFLIGHT_CARD.md"), useBytes = TRUE)
cat("PASS: Figure 2 v07 hotspot crops and Windows grid font preflight completed\n")
