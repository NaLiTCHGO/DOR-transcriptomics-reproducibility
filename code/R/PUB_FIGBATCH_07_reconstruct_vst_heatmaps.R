options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
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
get_optional_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args0, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

project_root <- normalizePath(get_arg("project-root"), winslash = "/", mustWork = TRUE)
run_root <- normalizePath(get_arg("run-root"), winslash = "/", mustWork = TRUE)
release_mode <- identical(toupper(get_optional_arg("release-mode", "FALSE")), "TRUE")
render_dpi <- if (release_mode) 600 else 300
output_root <- file.path(run_root, "reconstructed_heatmaps")
qc_root <- file.path(run_root, "qc")
manifest_root <- file.path(run_root, "manifests")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)

warning_rows <- list()
record_warning <- function(w) {
  warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    class = paste(class(w), collapse = ";"),
    message = conditionMessage(w),
    stringsAsFactors = FALSE
  )
  invokeRestart("muffleWarning")
}

specs <- list(
  list(
    figure = "S1", cohort = "GSE274832", m02_run = "20260813_M02_B1_GSE274832",
    config = "GSE274832_SAMPLES.csv", inclusion = NA_character_, design = "~ condition",
    min_count_samples = 3L,
    locked_effect = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
    locked_normalized = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_DESEQ2_NORMALIZED_COUNTS.csv.gz",
    output_png = "GSE274832_TOP30_VST_RECONSTRUCTED_ANGLE45.png",
    width_px = 1800L, height_px = 2200L,
    expected_columns = c("SRR30244481", "SRR30244480", "SRR30244479", "SRR30244484", "SRR30244483", "SRR30244482")
  ),
  list(
    figure = "S2", cohort = "GSE193136", m02_run = "20260814_M02_B3_GSE193136",
    config = "GSE193136_SAMPLES.csv",
    inclusion = "results/final_qc/GSE193136_M02_SAMPLE_INCLUSION.csv",
    design = "~ age_centered + condition", min_count_samples = 6L,
    locked_effect = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_DESEQ2_AGE_ADJUSTED_DOR_MINUS_NOR_ALL_GENES.csv",
    locked_normalized = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_DESEQ2_NORMALIZED_COUNTS.csv.gz",
    output_png = "GSE193136_TOP30_AGE_ADJUSTED_VST_RECONSTRUCTED_ANGLE45.png",
    width_px = 2000L, height_px = 2400L,
    expected_columns = c("SRR17463945", "SRR17463943", "SRR17463944", "SRR17463948", "SRR17463947", "SRR17463946",
                         "SRR17463953", "SRR17463951", "SRR17463950", "SRR17463952", "SRR17463954", "SRR17463949")
  ),
  list(
    figure = "S3", cohort = "GSE232306", m02_run = "20260814_M02_B4_GSE232306",
    config = "GSE232306_SAMPLES.csv",
    inclusion = "results/final_qc/GSE232306_M02_SAMPLE_INCLUSION.csv",
    design = "~ condition", min_count_samples = 6L,
    locked_effect = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
    locked_normalized = "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_DESEQ2_NORMALIZED_COUNTS.csv.gz",
    output_png = "GSE232306_TOP30_VST_RECONSTRUCTED_ANGLE45.png",
    width_px = 2000L, height_px = 2400L,
    expected_columns = c("SRR24505880", "SRR24505874", "SRR24505882", "SRR24505876", "SRR24505884", "SRR24505878",
                         "SRR24505886", "SRR24505896", "SRR24505892", "SRR24505888", "SRR24505894", "SRR24505890")
  )
)

if (release_mode) {
  specs <- lapply(specs, function(s) {
    s$width_px <- as.integer(s$width_px * 2L)
    s$height_px <- as.integer(s$height_px * 2L)
    s
  })
}

txmap_path <- file.path(project_root, "03_data/references/GENCODE_v50_GRCh38p14/gencode.v50.tx2gene.tsv.gz")
txmap <- read.delim(gzfile(txmap_path), check.names = FALSE)
txmap <- txmap[!duplicated(txmap$transcript_id_version), ]
tx2gene <- txmap[, c("transcript_id_version", "gene_id_version")]
gene_annot <- txmap[!duplicated(txmap$gene_id_version), c("gene_id_version", "gene_id", "gene_name")]

unit_sum <- function(x, idx) {
  if (!length(idx)) return(grid::unit(0, "in"))
  sum(x[idx])
}

validation_rows <- list()
geometry_rows <- list()

reconstruct_one <- function(s) {
  m02_root <- file.path(project_root, "05_analysis_steps/M02_COHORT_QC/runs", s$m02_run)
  cfg <- read.csv(file.path(project_root, "04_code/configs", s$config), check.names = FALSE)
  if (!is.na(s$inclusion)) {
    inclusion <- read.csv(file.path(m02_root, s$inclusion), check.names = FALSE)
    cfg <- cfg[cfg$srr %in% inclusion$srr[inclusion$include_m03 == "YES"], , drop = FALSE]
  }
  cfg$condition <- factor(cfg$phenotype, levels = c("NOR", "DOR"))
  if (s$cohort == "GSE193136") {
    cfg$age_years <- as.numeric(cfg$age_years)
    cfg$age_centered <- cfg$age_years - mean(cfg$age_years)
  }
  rownames(cfg) <- cfg$srr

  quant_files <- file.path(m02_root, "results/salmon", cfg$srr, "quant.sf")
  names(quant_files) <- cfg$srr
  if (!all(file.exists(quant_files))) stop("Missing Salmon quant.sf for ", s$cohort)

  txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                  countsFromAbundance = "lengthScaledTPM")
  cfg <- cfg[match(colnames(txi$counts), cfg$srr), , drop = FALSE]
  if (anyNA(cfg$srr) || !identical(cfg$srr, colnames(txi$counts))) {
    stop("Sample-order contract failed for ", s$cohort)
  }
  rownames(cfg) <- cfg$srr

  dds <- DESeqDataSetFromTximport(txi, colData = cfg, design = as.formula(s$design))
  keep <- rowSums(counts(dds) >= 10) >= s$min_count_samples
  dds <- dds[keep, ]
  mm <- model.matrix(design(dds), colData(dds))
  if (qr(mm)$rank != ncol(mm)) stop("Design matrix is not full rank for ", s$cohort)
  dds <- DESeq(dds, quiet = TRUE)

  locked_effect <- read.csv(file.path(project_root, s$locked_effect), check.names = FALSE)
  # Match the original formal M03 plotting scripts exactly: ascending nominal
  # p value, excluding NA values. No secondary effect-size tie breaker is used.
  locked_effect <- locked_effect[order(locked_effect$pvalue, na.last = NA), ]
  top_ids <- head(locked_effect$gene_id_version[is.finite(locked_effect$pvalue)], 30L)

  vst_matrix <- assay(vst(dds, blind = TRUE))
  if (!all(top_ids %in% rownames(vst_matrix))) stop("Locked Top-30 IDs missing from VST matrix for ", s$cohort)
  hm <- vst_matrix[top_ids, , drop = FALSE]
  gene_labels <- gene_annot$gene_name[match(rownames(hm), gene_annot$gene_id_version)]
  if (anyNA(gene_labels)) stop("Top-30 gene label mapping failed for ", s$cohort)
  rownames(hm) <- gene_labels

  annotation <- data.frame(Phenotype = cfg[colnames(hm), "phenotype"], row.names = colnames(hm))
  if (s$cohort == "GSE193136") annotation$Age <- cfg[colnames(hm), "age_years"]
  annotation_colors <- list(Phenotype = c(DOR = "#00BFC4", NOR = "#F8766D"))

  ph <- pheatmap(
    hm, scale = "row", annotation_col = annotation,
    annotation_colors = annotation_colors,
    show_colnames = TRUE, show_rownames = TRUE, angle_col = "45",
    main = NA, border_color = "#A8A8A8", fontsize = 8.2,
    fontsize_col = 7.4, fontsize_row = 7.2,
    treeheight_row = 42, treeheight_col = 42, silent = TRUE
  )
  observed_columns <- colnames(hm)[ph$tree_col$order]
  if (!identical(observed_columns, s$expected_columns)) {
    stop("Reconstructed column order differs from locked raster for ", s$cohort,
         ": observed=", paste(observed_columns, collapse = ","))
  }

  cohort_root <- file.path(output_root, s$figure)
  dir.create(cohort_root, recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(cohort_root, s$output_png)
  agg_png(output_png, width = s$width_px, height = s$height_px, units = "px", res = render_dpi,
          background = "white", scaling = 1)
  device_open <- TRUE
  on.exit({
    if (isTRUE(device_open) && grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  grid.newpage()
  grid.draw(ph$gtable)

  matrix_cell <- ph$gtable$layout[ph$gtable$layout$name == "matrix", , drop = FALSE]
  if (nrow(matrix_cell) != 1L) stop("Could not locate pheatmap matrix cell for ", s$cohort)
  total_width_in <- convertWidth(sum(ph$gtable$widths), "in", valueOnly = TRUE)
  left_in <- convertWidth(unit_sum(ph$gtable$widths, seq_len(matrix_cell$l - 1L)), "in", valueOnly = TRUE)
  right_in <- convertWidth(unit_sum(ph$gtable$widths, seq_len(matrix_cell$r)), "in", valueOnly = TRUE)
  device_width_in <- s$width_px / render_dpi
  centred_offset_in <- max(0, (device_width_in - total_width_in) / 2)
  matrix_left_px <- (centred_offset_in + left_in) * render_dpi
  matrix_right_px <- (centred_offset_in + right_in) * render_dpi
  invisible(dev.off())
  device_open <- FALSE

  if (!file.exists(output_png) || file.info(output_png)$size < 50000) {
    stop("Reconstructed heatmap PNG was not written for ", s$cohort)
  }

  scaled_hm <- t(scale(t(hm)))
  vst_df <- data.frame(gene_name = rownames(hm), hm, check.names = FALSE)
  scaled_df <- data.frame(gene_name = rownames(scaled_hm), scaled_hm, check.names = FALSE)
  write.csv(vst_df, gzfile(file.path(cohort_root, paste0(s$cohort, "_TOP30_VST_MATRIX.csv.gz"))), row.names = FALSE)
  write.csv(scaled_df, gzfile(file.path(cohort_root, paste0(s$cohort, "_TOP30_ROW_SCALED_VST_MATRIX.csv.gz"))), row.names = FALSE)
  write.csv(data.frame(sample = rownames(annotation), annotation, check.names = FALSE),
            file.path(cohort_root, paste0(s$cohort, "_HEATMAP_ANNOTATION.csv")), row.names = FALSE)
  write.csv(data.frame(position = seq_along(observed_columns), sample = observed_columns),
            file.path(cohort_root, paste0(s$cohort, "_LOCKED_COLUMN_ORDER.csv")), row.names = FALSE)

  locked_norm <- read.csv(gzfile(file.path(project_root, s$locked_normalized)), check.names = FALSE)
  new_norm <- counts(dds, normalized = TRUE)
  common_ids <- intersect(rownames(new_norm), locked_norm$gene_id_version)
  locked_idx <- match(common_ids, locked_norm$gene_id_version)
  sample_cols <- colnames(new_norm)
  if (!all(sample_cols %in% names(locked_norm))) stop("Locked normalized-count sample columns missing for ", s$cohort)
  locked_values <- as.matrix(locked_norm[locked_idx, sample_cols, drop = FALSE])
  storage.mode(locked_values) <- "double"
  new_values <- new_norm[common_ids, sample_cols, drop = FALSE]
  max_abs_diff <- max(abs(new_values - locked_values), na.rm = TRUE)
  if (!is.finite(max_abs_diff) || max_abs_diff > 1e-6) {
    stop("Normalized-count reconstruction differs from locked result for ", s$cohort,
         "; max abs diff=", format(max_abs_diff, scientific = TRUE))
  }

  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    figure = s$figure, cohort = s$cohort, samples = ncol(hm), genes = nrow(hm),
    locked_top30_ids = all(top_ids == rownames(vst_matrix)[match(top_ids, rownames(vst_matrix))]),
    locked_column_order = identical(observed_columns, s$expected_columns),
    normalized_count_common_genes = length(common_ids),
    normalized_count_max_abs_diff = max_abs_diff,
    output_png = output_png, status = "PASS",
    stringsAsFactors = FALSE
  )
  geometry_rows[[length(geometry_rows) + 1L]] <<- data.frame(
    figure = s$figure, cohort = s$cohort, image_width_px = s$width_px,
    image_height_px = s$height_px, matrix_left_px = matrix_left_px,
    matrix_right_px = matrix_right_px,
    scientific_body_anchor_fraction = matrix_left_px / s$width_px,
    stringsAsFactors = FALSE
  )
}

withCallingHandlers({
  for (s in specs) reconstruct_one(s)
}, warning = record_warning)

warning_log <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
  timestamp = character(), class = character(), message = character(), stringsAsFactors = FALSE
)
write_csv(warning_log, file.path(qc_root, "RUNTIME_WARNING_LOG_VST_RECONSTRUCTION.csv"))
write_csv(do.call(rbind, validation_rows), file.path(qc_root, "VST_HEATMAP_RECONSTRUCTION_VALIDATION.csv"))
write_csv(do.call(rbind, geometry_rows), file.path(output_root, "VST_HEATMAP_GEOMETRY.csv"))

if (nrow(warning_log) > 0L) stop("VST reconstruction emitted warnings; see runtime warning log")

provenance <- c(
  "# VST heatmap display reconstruction provenance",
  "",
  "- Scope: Figure S1-S3 panel D only.",
  "- Reason: the original Top-30 VST plotting matrices were not serialized, although all Salmon quantifications, sample tables, reference map, M03 scripts, locked effect tables and locked normalized counts remain available.",
  "- Top-30 selection: read from the locked M03 effect tables; no newly calculated p value is used for selection.",
  "- Reconstruction: the original tximport lengthScaledTPM, DESeq2 design/filter and blind VST route is recreated in memory solely to rebuild the display matrix.",
  "- Input identity gate: reconstructed DESeq2 normalized counts must match the locked normalized-count tables with maximum absolute difference <= 1e-6.",
  "- Column-order gate: reconstructed pheatmap clustering must exactly match the column order visible in the locked original raster.",
  "- New archival outputs: unscaled Top-30 VST matrix, row-scaled plotting matrix, sample annotation, column order, geometry and a 45-degree native-label heatmap PNG for each cohort.",
  paste0("- Render mode: ", if (release_mode) "release 600 ppi" else "candidate 300 ppi", "."),
  "- Locked M03 results overwritten: NO.",
  "- Differential, LOSO, meta-analysis, pathway or LOCO results recomputed for claims: NO."
)
writeLines(provenance, file.path(manifest_root, "VST_HEATMAP_RECONSTRUCTION_PROVENANCE.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(manifest_root, "R_SESSION_INFO_VST_RECONSTRUCTION.txt"), useBytes = TRUE)
cat("PASS: S1-S3 Top-30 VST display matrices reconstructed at ", render_dpi,
    " ppi; locked normalized counts and column orders matched\n", sep = "")
