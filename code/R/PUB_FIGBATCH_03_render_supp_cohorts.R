options(stringsAsFactors = FALSE)

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
skill_root <- normalizePath(get_arg("skill-root"), winslash = "/", mustWork = TRUE)
release_mode <- identical(toupper(get_optional_arg("release-mode", "FALSE")), "TRUE")
source(file.path(project_root, "04_code/R/PUB_FIGBATCH_00_common.R"), local = FALSE)
theme_path <- activate_skill3_theme(skill_root)
suppressPackageStartupMessages(library(magick))

source_root <- file.path(run_root, "source_data")
reconstructed_root <- file.path(run_root, "reconstructed_heatmaps")
figure_root <- file.path(run_root, "figures", "supplement")
qc_root <- file.path(run_root, "qc")
manifest_root <- file.path(run_root, "manifests")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)

cohort_spec <- data.frame(
  figure = c("S1", "S2", "S3"),
  cohort = c("GSE274832", "GSE193136", "GSE232306"),
  heatmap = c(
    "GSE274832_TOP30_VST_RECONSTRUCTED_ANGLE45.png",
    "GSE193136_TOP30_AGE_ADJUSTED_VST_RECONSTRUCTED_ANGLE45.png",
    "GSE232306_TOP30_VST_RECONSTRUCTED_ANGLE45.png"
  ),
  title_suffix = c("unadjusted", "age-adjusted", "unadjusted"),
  stringsAsFactors = FALSE
)

geometry_path <- file.path(reconstructed_root, "VST_HEATMAP_GEOMETRY.csv")
if (!file.exists(geometry_path)) stop("Missing reconstructed VST geometry manifest: ", geometry_path)
heatmap_geometry <- read_csv(geometry_path, show_col_types = FALSE)

read_corr_long <- function(path, sample_levels) {
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "double"
  d <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(d) <- c("row_sample", "col_sample", "correlation")
  d$row_sample <- factor(d$row_sample, levels = rev(sample_levels), labels = rev(short_sample(sample_levels)))
  d$col_sample <- factor(d$col_sample, levels = sample_levels, labels = short_sample(sample_levels))
  d
}

reconstructed_raster_body <- function(path, geometry) {
  img <- image_read(path)
  info <- image_info(img)
  if (nrow(geometry) != 1L || info$width[[1]] != geometry$image_width_px[[1]] ||
      info$height[[1]] != geometry$image_height_px[[1]]) {
    stop("Reconstructed heatmap geometry contract failed: ", path)
  }
  rg <- grid::rasterGrob(as.raster(img), width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
                         interpolate = TRUE)
  p <- ggplot() +
    annotation_custom(rg, xmin = 0, xmax = info$width[[1]], ymin = 0, ymax = 1) +
    coord_cartesian(xlim = c(0, info$width[[1]]), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(0, 0, 0.8, 0, unit = "mm"))
  attr(p, "scientific_body_anchor_fraction") <- geometry$scientific_body_anchor_fraction[[1]]
  p
}

build_supp_cohort <- function(fig_id, cohort, heatmap_name, model_label) {
  group <- fig_id
  marker_path <- file.path(source_root, group, paste0(cohort, "_MARKER_EXPRESSION.csv"))
  corr_path <- file.path(source_root, group, paste0(cohort, "_SAMPLE_CORRELATION.csv"))
  loso_path <- file.path(source_root, group, paste0(cohort, "_LEAVE_ONE_SAMPLE_OUT_GLOBAL_QC.csv"))
  heatmap_path <- file.path(reconstructed_root, group, heatmap_name)
  if (!all(file.exists(c(marker_path, corr_path, loso_path, heatmap_path)))) stop("Incomplete source group: ", fig_id)
  markers <- read_csv(marker_path, show_col_types = FALSE, progress = FALSE)
  loso <- read_csv(loso_path, show_col_types = FALSE, progress = FALSE)
  sample_info <- markers %>% distinct(srr, phenotype) %>%
    mutate(phenotype = factor(phenotype, levels = c("DOR", "NOR"))) %>% arrange(phenotype, srr)
  sample_levels <- sample_info$srr
  gene_info <- markers %>% distinct(marker_group, gene_name) %>% arrange(marker_group, gene_name)
  gene_levels <- gene_info$gene_name
  markers <- markers %>% mutate(
    srr_short = factor(short_sample(srr), levels = short_sample(sample_levels)),
    gene_name = factor(gene_name, levels = rev(gene_levels)),
    marker_group = factor(marker_group, levels = unique(gene_info$marker_group))
  )

  p_marker <- ggplot(markers, aes(srr_short, gene_name, fill = log2_TPM_plus_1)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    facet_grid(marker_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = c("white", "#FEE08B", "#D55E00"), name = "log2(TPM+1)",
                         breaks = pretty(range(markers$log2_TPM_plus_1), n = 3)) +
    labs(x = "Sample", y = NULL) + theme_pub(7.0) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 6.2),
          axis.text.y = element_text(size = 6.1), strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 6.3), legend.position = "bottom",
          legend.box.spacing = grid::unit(0.75, "mm"), legend.margin = margin(0.5, 0, 0, 0, unit = "mm"),
          plot.margin = margin(2.5, 2.5, 1.5, 2.5, unit = "mm")) +
    guides(fill = guide_colourbar(direction = "horizontal", barwidth = grid::unit(24, "mm"),
                                  barheight = grid::unit(2.2, "mm"), title.position = "top"))

  corr_long <- read_corr_long(corr_path, sample_levels)
  corr_min <- min(corr_long$correlation, na.rm = TRUE)
  boundary <- sum(sample_info$phenotype == "DOR") + 0.5
  p_corr <- ggplot(corr_long, aes(col_sample, row_sample, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    geom_vline(xintercept = boundary, colour = "#222222", linewidth = 0.35) +
    geom_hline(yintercept = length(sample_levels) - boundary + 1, colour = "#222222", linewidth = 0.35) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = (corr_min + 1) / 2, limits = c(corr_min, 1), name = "Correlation",
                         breaks = round(seq(corr_min, 1, length.out = 3), 2)) +
    labs(x = "Sample", y = "Sample") + coord_equal() + theme_pub(7.0) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 6.0),
          axis.text.y = element_text(size = 6.0), legend.position = "bottom",
          legend.box.spacing = grid::unit(0.75, "mm"), legend.margin = margin(0.5, 0, 0, 0, unit = "mm"),
          plot.margin = margin(2.5, 2.5, 1.5, 2.5, unit = "mm")) +
    guides(fill = guide_colourbar(direction = "horizontal", barwidth = grid::unit(24, "mm"),
                                  barheight = grid::unit(2.2, "mm"), title.position = "top"))

  loso <- loso %>% mutate(
    held_out_phenotype = factor(held_out_phenotype, levels = c("DOR", "NOR")),
    sample_label = factor(short_sample(held_out_srr),
                          levels = rev(short_sample(held_out_srr[order(spearman_lfc_vs_full)])))
  )
  p_loso <- ggplot(loso, aes(sample_label, spearman_lfc_vs_full, fill = held_out_phenotype)) +
    geom_col(width = 0.72) +
    geom_point(aes(y = same_direction_fraction), shape = 21, fill = "white", colour = "#202020", size = 1.5, stroke = 0.35) +
    coord_flip(ylim = c(0, 1)) +
    scale_fill_manual(values = phenotype_palette, name = NULL) +
    labs(x = "Held-out sample", y = "Correlation / direction fraction") + theme_pub(7.0) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 6.2))

  geometry <- heatmap_geometry[heatmap_geometry$cohort == cohort, , drop = FALSE]
  p_heatmap <- reconstructed_raster_body(heatmap_path, geometry)
  a <- panel_block(p_marker, "A", "Marker-expression sanity check",
                   "Granulosa and off-target markers by sample")
  b <- panel_block(p_corr, "B", "Sample-correlation structure",
                   "DOR then NOR; dark rules mark the group boundary")
  c <- panel_block(p_loso, "C", "Leave-one-sample-out stability",
                   "Bars: LFC rho; open points: direction agreement")
  d <- panel_block(p_heatmap, "D", "Exploratory Top-30 VST heatmap",
                   paste0("Rebuilt ", model_label, " display; native 45-degree labels"))
  # coord_equal() preserves square correlation cells but gives the B block a
  # slightly shorter intrinsic height, so patchwork centres it about 2 mm
  # below A. Shift the complete B block upward inside the same allocation;
  # scientific geometry, labels and colour scale are unchanged.
  b <- plot_spacer() + inset_element(
    b, left = 0, bottom = 0.027, right = 1, top = 1.027,
    align_to = "full", clip = FALSE
  )
  # One shared 1/28-page spacer shifts C and D together. The C plot spine then
  # aligns with A, while the locked D matrix border (excluding its dendrogram)
  # aligns with the B plot spine. Scientific bodies keep their locked order.
  wrap_plots(list(a, b, c, d),
             design = "AAAAAAAAAAAAAABBBBBBBBBBBBBB\n#CCCCCCCCCCCCCCDDDDDDDDDDDDD",
             heights = c(1.0, 1.02)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

warning_csv <- file.path(qc_root, "RUNTIME_WARNING_LOG_SUPP_COHORTS.csv")
reconciliation <- file.path(qc_root, "FULL_RUNTIME_WARNING_RECONCILIATION_SUPP_COHORTS.md")
release_rows <- list()
capture_warnings({
  for (i in seq_len(nrow(cohort_spec))) {
    s <- cohort_spec[i, ]
    if (release_mode) {
      figure_id <- paste0("Figure_", s$figure)
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, figure_id), figure_id, 2008, 2600,
        function() build_supp_cohort(s$figure, s$cohort, s$heatmap, s$title_suffix)
      )
    } else {
      out <- file.path(figure_root, sprintf("Figure_%s_%s_QC_and_Stability_v08_QC.png", s$figure, s$cohort))
      save_candidate(out, 2008, 2600, function() build_supp_cohort(s$figure, s$cohort, s$heatmap, s$title_suffix))
    }
  }
}, warning_csv, reconciliation, "SUPPLEMENTARY_COHORT_FIGURES_S1_TO_S3")
if (release_mode) {
  write_csv(bind_rows(release_rows), file.path(manifest_root, "FINAL_EXPORT_MANIFEST_SUPP_COHORTS.csv"))
}

provenance <- c(
  "# Cohort-supplement renderer provenance",
  "",
  paste0("- Active Skill root: `", skill_root, "`"),
  paste0("- Theme module: `", theme_path, "`"),
  "- Active suite: `research_scientific_figure_suite v2.2-beta.9`",
  "- Panels A-C: redrawn from locked M02/M03 CSV tables",
  "- Panels A-B: sample labels use a 55-degree right-justified route; horizontal colour keys are pulled upward without changing scales.",
  "- Panel D: display reconstructed from preserved Salmon quantifications and the original tximport/DESeq2 blind-VST route; Top-30 IDs come from locked M03 effects; native pheatmap 45-degree column labels share the matrix layout and cannot drift from their columns.",
  "- Reconstruction identity gates: normalized counts match locked tables within 1e-6 and clustered sample order exactly matches the locked raster; matrices and orders are now archived under reconstructed_heatmaps/.",
  "- Row geometry: C/D shifted together by 1/28 page width so C aligns to A and the D matrix border (excluding dendrogram) aligns to B.",
  "- Scientific claim/result recomputation: `NO`; display VST matrices are identity-checked reconstructions.",
  paste0("- Output mode: `", if (release_mode) "RELEASE_FIGURE" else "CANDIDATE_300DPI", "`"),
  if (release_mode) {
    "- Outputs: independently rendered 300-dpi preview, true 600-ppi PNG, LZW TIFF and physical-size Cairo PDF for each figure."
  } else {
    "- Outputs: three 300-dpi RGB QC PNG candidates only"
  }
)
writeLines(provenance, file.path(manifest_root, "RENDERER_PROVENANCE_SUPP_COHORTS.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(manifest_root, "R_SESSION_INFO_SUPP_COHORTS.txt"), useBytes = TRUE)
cat("PASS: Figure S1-S3 figures rendered with zero captured warnings; mode=",
    if (release_mode) "RELEASE" else "CANDIDATE", "\n", sep = "")
