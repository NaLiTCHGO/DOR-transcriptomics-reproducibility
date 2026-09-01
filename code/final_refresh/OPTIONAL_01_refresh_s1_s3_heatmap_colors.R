options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(pheatmap)
  library(ragg)
  library(readr)
  library(grid)
})

args0 <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key) {
  hit <- grep(paste0("^--", key, "="), args0, value = TRUE)
  if (!length(hit)) stop("Missing --", key, "=<path>")
  sub(paste0("^--", key, "="), "", hit[[1]])
}

source_run <- normalizePath(get_arg("source-run"), winslash = "/", mustWork = TRUE)
run_root <- normalizePath(get_arg("run-root"), winslash = "/", mustWork = TRUE)
render_dpi <- 600L
output_root <- file.path(run_root, "reconstructed_heatmaps")
qc_root <- file.path(run_root, "qc")
manifest_root <- file.path(run_root, "manifests")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)

specs <- list(
  list(figure = "S1", cohort = "GSE274832",
       matrix = "GSE274832_TOP30_VST_MATRIX.csv.gz",
       annotation = "GSE274832_HEATMAP_ANNOTATION.csv",
       order = "GSE274832_LOCKED_COLUMN_ORDER.csv",
       output = "GSE274832_TOP30_VST_RECONSTRUCTED_ANGLE45.png",
       width = 3600L, height = 4400L),
  list(figure = "S2", cohort = "GSE193136",
       matrix = "GSE193136_TOP30_VST_MATRIX.csv.gz",
       annotation = "GSE193136_HEATMAP_ANNOTATION.csv",
       order = "GSE193136_LOCKED_COLUMN_ORDER.csv",
       output = "GSE193136_TOP30_AGE_ADJUSTED_VST_RECONSTRUCTED_ANGLE45.png",
       width = 4000L, height = 4800L),
  list(figure = "S3", cohort = "GSE232306",
       matrix = "GSE232306_TOP30_VST_MATRIX.csv.gz",
       annotation = "GSE232306_HEATMAP_ANNOTATION.csv",
       order = "GSE232306_LOCKED_COLUMN_ORDER.csv",
       output = "GSE232306_TOP30_VST_RECONSTRUCTED_ANGLE45.png",
       width = 4000L, height = 4800L)
)

unit_sum <- function(x, idx) {
  if (!length(idx)) return(grid::unit(0, "in"))
  sum(x[idx])
}

geometry_rows <- list()
validation_rows <- list()

render_one <- function(s) {
  src <- file.path(source_run, "reconstructed_heatmaps", s$figure)
  dst <- file.path(output_root, s$figure)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)

  m <- read.csv(gzfile(file.path(src, s$matrix)), check.names = FALSE)
  gene_name <- m$gene_name
  m$gene_name <- NULL
  hm <- as.matrix(m)
  storage.mode(hm) <- "double"
  rownames(hm) <- gene_name

  ann <- read.csv(file.path(src, s$annotation), check.names = FALSE)
  rownames(ann) <- ann$sample
  ann$sample <- NULL
  ann <- ann[colnames(hm), , drop = FALSE]
  if (anyNA(rownames(ann)) || !identical(rownames(ann), colnames(hm))) {
    stop("Annotation/sample order contract failed for ", s$cohort)
  }

  # Project-wide locked phenotype palette used by the main figures.
  annotation_colors <- list(Phenotype = c(DOR = "#D55E00", NOR = "#0072B2"))
  ph <- pheatmap(
    hm, scale = "row", annotation_col = ann,
    annotation_colors = annotation_colors,
    show_colnames = TRUE, show_rownames = TRUE, angle_col = "45",
    main = NA, border_color = "#A8A8A8", fontsize = 8.2,
    fontsize_col = 7.4, fontsize_row = 7.2,
    treeheight_row = 42, treeheight_col = 42, silent = TRUE
  )

  observed <- colnames(hm)[ph$tree_col$order]
  expected <- read.csv(file.path(src, s$order), check.names = FALSE)$sample
  if (!identical(observed, expected)) {
    stop("Clustered column order changed for ", s$cohort)
  }

  output_png <- file.path(dst, s$output)
  agg_png(output_png, width = s$width, height = s$height, units = "px",
          res = render_dpi, background = "white", scaling = 1)
  grid.newpage()
  grid.draw(ph$gtable)

  matrix_cell <- ph$gtable$layout[ph$gtable$layout$name == "matrix", , drop = FALSE]
  if (nrow(matrix_cell) != 1L) stop("Matrix geometry cell missing for ", s$cohort)
  total_width_in <- convertWidth(sum(ph$gtable$widths), "in", valueOnly = TRUE)
  left_in <- convertWidth(unit_sum(ph$gtable$widths, seq_len(matrix_cell$l - 1L)), "in", valueOnly = TRUE)
  right_in <- convertWidth(unit_sum(ph$gtable$widths, seq_len(matrix_cell$r)), "in", valueOnly = TRUE)
  device_width_in <- s$width / render_dpi
  centred_offset_in <- max(0, (device_width_in - total_width_in) / 2)
  matrix_left_px <- (centred_offset_in + left_in) * render_dpi
  matrix_right_px <- (centred_offset_in + right_in) * render_dpi
  invisible(dev.off())

  if (!file.exists(output_png) || file.info(output_png)$size < 50000) {
    stop("Output heatmap was not written for ", s$cohort)
  }

  # Preserve the exact archived data objects beside the refreshed raster.
  file.copy(file.path(src, s$matrix), file.path(dst, s$matrix), overwrite = TRUE)
  scaled_name <- sub("_VST_MATRIX", "_ROW_SCALED_VST_MATRIX", s$matrix, fixed = TRUE)
  file.copy(file.path(src, scaled_name), file.path(dst, scaled_name), overwrite = TRUE)
  file.copy(file.path(src, s$annotation), file.path(dst, s$annotation), overwrite = TRUE)
  file.copy(file.path(src, s$order), file.path(dst, s$order), overwrite = TRUE)

  geometry_rows[[length(geometry_rows) + 1L]] <<- data.frame(
    figure = s$figure, cohort = s$cohort,
    image_width_px = s$width, image_height_px = s$height,
    matrix_left_px = matrix_left_px, matrix_right_px = matrix_right_px,
    scientific_body_anchor_fraction = matrix_left_px / s$width
  )
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    figure = s$figure, cohort = s$cohort,
    genes = nrow(hm), samples = ncol(hm),
    matrix_values_reused = TRUE,
    clustered_column_order_identical = TRUE,
    DOR_color = annotation_colors$Phenotype[["DOR"]],
    NOR_color = annotation_colors$Phenotype[["NOR"]]
  )
}

for (s in specs) render_one(s)
write_csv(do.call(rbind, geometry_rows), file.path(output_root, "VST_HEATMAP_GEOMETRY.csv"))
write_csv(do.call(rbind, validation_rows), file.path(qc_root, "OPTIONAL_S1_S3_HEATMAP_COLOR_REFRESH_VALIDATION.csv"))
writeLines(capture.output(sessionInfo()), file.path(manifest_root, "R_SESSION_INFO_OPTIONAL_HEATMAP_REFRESH.txt"))
writeLines(c(
  "# Optional S1-S3 heatmap-colour refresh",
  "",
  paste0("- Source locked figure run: `", source_run, "`"),
  "- Statistical recomputation: `NO`",
  "- Matrix, annotation and locked clustered column order: reused unchanged",
  "- Display-only change: DOR #D55E00; NOR #0072B2",
  "- Output: independent 600-ppi heatmap rasters for the optional refreshed S1-S3 figures"
), file.path(manifest_root, "OPTIONAL_S1_S3_HEATMAP_COLOR_REFRESH_PROVENANCE.md"))
cat("PASS: optional S1-S3 heatmap colours refreshed without changing matrices or column order.\n")
