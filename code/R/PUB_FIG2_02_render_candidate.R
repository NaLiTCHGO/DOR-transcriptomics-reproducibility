#!/usr/bin/env Rscript

# Render the current Skill 3 manuscript Figure 2 from direct plotting tables.
# Default mode preserves the candidate-first 300-dpi QC route. Explicit
# --release-mode TRUE creates the independently rendered final format bundle.

arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing required argument: ", flag)
  args[[idx + 1L]]
}
optional_arg_value <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), winslash = "/", mustWork = TRUE)
run_dir <- normalizePath(arg_value("--run-dir"), winslash = "/", mustWork = TRUE)
skill_root <- normalizePath(arg_value("--skill-root"), winslash = "/", mustWork = TRUE)
release_mode <- identical(toupper(optional_arg_value("--release-mode", "FALSE")), "TRUE")
source_dir <- file.path(run_dir, "source_data")
figure_dir <- file.path(run_dir, "figures")
qc_dir <- file.path(run_dir, "qc")
manifest_dir <- file.path(run_dir, "manifests")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

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

write_warning_log <- function() {
  out <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
    timestamp = character(), class = character(), message = character(), stringsAsFactors = FALSE
  )
  readr::write_csv(out, file.path(qc_dir, "RUNTIME_WARNING_LOG.csv"), na = "")
  nrow(out)
}

failure_path <- file.path(qc_dir, "RENDER_FAILURE.md")

tryCatch(
  withCallingHandlers({
    suppressPackageStartupMessages({
      library(readr)
      library(dplyr)
      library(ggplot2)
      library(patchwork)
      library(gtable)
      library(ragg)
      library(scales)
    })

    theme_path <- file.path(skill_root, "R", "theme_research_publication.R")
    if (!file.exists(theme_path)) stop("Missing active Skill 3 theme: ", theme_path)
    source(theme_path, local = FALSE)
    source(file.path(project_root, "04_code/R/PUB_FIGBATCH_00_common.R"), local = FALSE)

    pca_path <- file.path(source_dir, "FIG2_PCA_DIRECT_SOURCE.csv")
    volcano_path <- file.path(source_dir, "FIG2_VOLCANO_DIRECT_SOURCE.csv")
    summary_path <- file.path(source_dir, "FIG2_PANEL_SUMMARY.csv")
    if (!all(file.exists(c(pca_path, volcano_path, summary_path)))) stop("Direct source package is incomplete")
    pca <- read_csv(pca_path, show_col_types = FALSE)
    volcano <- read_csv(volcano_path, show_col_types = FALSE)
    panel_summary <- read_csv(summary_path, show_col_types = FALSE)
    if (nrow(pca) != 30L || nrow(volcano) != 48732L || nrow(panel_summary) != 3L) {
      stop("Direct source row-count contract failed")
    }

    pca$phenotype <- factor(pca$phenotype, levels = c("NOR", "DOR"))
    volcano$display_status <- factor(
      volcano$display_status,
      levels = c("Other", "Meets FDR and effect threshold")
    )

    palette <- c(NOR = "#0072B2", DOR = "#D55E00")
    volcano_palette <- c("Other" = "#D9D9D9", "Meets FDR and effect threshold" = "#D55E00")

    panel_theme <- function() {
      theme_research_publication(base_size = 7.5, base_family = "Arial") +
        theme(
          plot.title = element_text(size = 8.5, face = "bold", hjust = 0, lineheight = 0.96,
                                    margin = margin(b = 0.7, unit = "mm")),
          plot.subtitle = element_text(size = 7.2, hjust = 0, lineheight = 1.02,
                                       margin = margin(b = 1.2, unit = "mm")),
          plot.title.position = "panel",
          axis.title = element_text(size = 7.5),
          axis.text = element_text(size = 7.0, colour = "black"),
          plot.margin = margin(4.8, 3.2, 2.2, 5.8, unit = "mm"),
          legend.position = "none"
        )
    }

    add_hanging_tag <- function(plot, tag, title) {
      gt <- ggplotGrob(plot)
      idx <- which(gt$layout$name == "title")
      if (length(idx) != 1L) stop("Expected exactly one title grob for panel ", tag)
      # Replace the independently measured ggplot title and overlaid panel tag
      # with one shared header grob. Identical font metrics and a common y/just
      # anchor make the tag and accession title truly collinear rather than
      # merely top-aligned.
      tag_grob <- grid::textGrob(
        tag,
        x = grid::unit(0, "npc") - grid::unit(2.0, "mm"),
        y = grid::unit(0.50, "npc"),
        just = c("right", "centre"),
        gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 8.5,
                        col = "#202020")
      )
      title_grob <- grid::textGrob(
        title,
        x = grid::unit(0, "npc"),
        y = grid::unit(0.50, "npc"),
        just = c("left", "centre"),
        gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 8.5,
                        col = "#202020")
      )
      gt$grobs[[idx]] <- grid::grobTree(tag_grob, title_grob)
      gt$layout$clip[idx] <- "off"
      patchwork::wrap_ggplot_grob(gt)
    }

    make_pca <- function(cohort) {
      d <- filter(pca, .data$cohort == .env$cohort)
      s <- filter(panel_summary, .data$cohort == .env$cohort)
      if (nrow(d) == 0L || nrow(s) != 1L) stop("PCA source mismatch for ", cohort)
      xpad <- max(diff(range(d$PC1)) * 0.08, 0.5)
      ypad <- max(diff(range(d$PC2)) * 0.10, 0.5)
      ggplot(d, aes(PC1, PC2)) +
        geom_hline(yintercept = 0, linewidth = 0.25, colour = "#E6E6E6") +
        geom_vline(xintercept = 0, linewidth = 0.25, colour = "#E6E6E6") +
        geom_point(aes(fill = phenotype), shape = 21, size = 2.35, stroke = 0.4,
                   colour = "white", alpha = 0.94) +
        scale_fill_manual(values = palette, drop = FALSE) +
        coord_cartesian(
          xlim = range(d$PC1) + c(-xpad, xpad),
          ylim = range(d$PC2) + c(-ypad, ypad),
          clip = "off"
        ) +
        labs(
          title = cohort,
          subtitle = sprintf("Expression PCA; %d DOR + %d NOR", s$n_dor, s$n_nor),
          x = sprintf("PC1 (%.1f%%)", s$pc1_variance_percent),
          y = sprintf("PC2 (%.1f%%)", s$pc2_variance_percent)
        ) +
        panel_theme()
    }

    make_volcano <- function(cohort) {
      d <- filter(volcano, .data$cohort == .env$cohort, .data$plot_eligible)
      s <- filter(panel_summary, .data$cohort == .env$cohort)
      if (nrow(d) == 0L || nrow(s) != 1L) stop("Volcano source mismatch for ", cohort)
      d$y_plot <- d$neg_log10_nominal_p
      d$capped <- FALSE
      cap_note <- NULL
      y_upper <- max(d$y_plot, na.rm = TRUE) * 1.06
      if (s$display_route == "EXTREME_Y_CAPPED_NO_LABEL_REVIEW") {
        y_cap <- s$suggested_y_cap
        if (!is.finite(y_cap) || y_cap <= 0) stop("Invalid y cap for ", cohort)
        d$capped <- d$neg_log10_nominal_p > y_cap
        d$y_plot <- pmin(d$neg_log10_nominal_p, y_cap)
        y_upper <- y_cap * 1.14
        cap_note <- sprintf("%d points above y=%.1f", sum(d$capped), y_cap)
      }
      regular <- filter(d, !.data$capped)
      capped <- filter(d, .data$capped)
      subtitle <- sprintf("%s; %s candidates", s$model, comma(s$candidate_count, accuracy = 1))
      p <- ggplot() +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3,
                   colour = "#7A7A7A") +
        geom_point(
          data = filter(regular, .data$display_status == "Other"),
          aes(log2FoldChange, y_plot, colour = display_status),
          size = 0.42, alpha = 0.42, stroke = 0
        ) +
        geom_point(
          data = filter(regular, .data$display_status != "Other"),
          aes(log2FoldChange, y_plot, colour = display_status),
          size = 0.58, alpha = 0.72, stroke = 0
        ) +
        scale_colour_manual(values = volcano_palette, drop = FALSE) +
        coord_cartesian(ylim = c(0, y_upper), clip = "off") +
        labs(
          title = cohort,
          subtitle = subtitle,
          x = expression(log[2]~fold~change),
          y = expression(-log[10]~nominal~italic(P))
        ) +
        panel_theme()
      if (nrow(capped)) {
        p <- p +
          geom_point(
            data = capped,
            aes(log2FoldChange, y_plot, fill = display_status),
            shape = 24, size = 1.45, stroke = 0.25, colour = "black"
          ) +
          scale_fill_manual(values = volcano_palette, guide = "none") +
          annotate("text", x = -Inf, y = y_upper, label = cap_note, hjust = -0.02, vjust = 1.1,
                   family = "Arial", size = 2.15, colour = "#333333")
      }
      p
    }

    legend_strip <- ggplot() +
      annotate("text", x = 0.02, y = 0.68, label = "PCA:", hjust = 0,
               family = "Arial", fontface = "bold", size = 2.55) +
      annotate("point", x = 0.11, y = 0.68, shape = 21, size = 2.8, fill = palette[["DOR"]],
               colour = "white", stroke = 0.35) +
      annotate("text", x = 0.135, y = 0.68, label = "DOR", hjust = 0,
               family = "Arial", size = 2.45) +
      annotate("point", x = 0.23, y = 0.68, shape = 21, size = 2.8, fill = palette[["NOR"]],
               colour = "white", stroke = 0.35) +
      annotate("text", x = 0.255, y = 0.68, label = "NOR", hjust = 0,
               family = "Arial", size = 2.45) +
      # The volcano key belongs to the B/D/F (right-column) panels. Keep the
      # whole key visually centred beneath that column instead of the page.
      annotate("text", x = 0.616, y = 0.68, label = "Volcano:", hjust = 0,
               family = "Arial", fontface = "bold", size = 2.55) +
      annotate("point", x = 0.740, y = 0.68, shape = 16, size = 2.0,
               colour = volcano_palette[["Meets FDR and effect threshold"]]) +
      annotate("text", x = 0.765, y = 0.68, label = "FDR < 0.05 and |log2FC| >= 1", hjust = 0,
               family = "Arial", size = 2.45) +
      # Keep the grey point and its label outside the threshold-text footprint.
      annotate("point", x = 1.075, y = 0.68, shape = 16, size = 2.0,
               colour = volcano_palette[["Other"]]) +
      annotate("text", x = 1.10, y = 0.68, label = "Other", hjust = 0,
               family = "Arial", size = 2.45) +
      coord_cartesian(xlim = c(0, 1.13), ylim = c(0, 1), clip = "off") +
      theme_void(base_family = "Arial") +
      theme(plot.margin = margin(0, 4, 0, 4, unit = "mm"))

    build_composite <- function() {
      panels <- list(
        add_hanging_tag(make_pca("GSE274832"), "A", "GSE274832"),
        add_hanging_tag(make_volcano("GSE274832"), "B", "GSE274832"),
        add_hanging_tag(make_pca("GSE193136"), "C", "GSE193136"),
        add_hanging_tag(make_volcano("GSE193136"), "D", "GSE193136"),
        add_hanging_tag(make_pca("GSE232306"), "E", "GSE232306"),
        add_hanging_tag(make_volcano("GSE232306"), "F", "GSE232306")
      )
      (panels[[1]] + panels[[2]]) /
        (panels[[3]] + panels[[4]]) /
        (panels[[5]] + panels[[6]]) /
        legend_strip +
        plot_layout(heights = c(1, 1, 1, 0.12), widths = c(1, 1))
    }

    if (release_mode) {
      release_manifest <- save_release_bundle(
        file.path(figure_dir, "main", "Figure_2"), "Figure_2",
        2008L, 2409L, build_composite
      )
      readr::write_csv(release_manifest, file.path(manifest_dir, "FINAL_EXPORT_MANIFEST_FIG2.csv"), na = "")
    } else {
      output_path <- file.path(figure_dir, "Figure_2_Cohort_QC_and_Effects_v07_QC.png")
      # Open the verified ragg device before gtable construction. This ensures
      # all Arial string metrics are resolved by the target renderer rather than
      # by R's fallback PostScript font database.
      agg_png(output_path, width = 2008, height = 2409, units = "px", res = 300,
              background = "white", scaling = 1)
      print(build_composite())
      invisible(dev.off())
      if (!file.exists(output_path) || file.info(output_path)$size <= 0) stop("Candidate PNG was not written")
    }

    provenance <- c(
      "# Figure 2 renderer provenance",
      "",
      paste0("- Active Skill root: `", skill_root, "`"),
      paste0("- Theme module: `", theme_path, "`"),
      "- Active suite version: `research_scientific_figure_suite v2.2-beta.9`",
      paste0("- R: `", R.version.string, "`"),
      paste0("- ggplot2: `", as.character(packageVersion("ggplot2")), "`"),
      paste0("- patchwork: `", as.character(packageVersion("patchwork")), "`"),
      paste0("- ragg: `", as.character(packageVersion("ragg")), "`"),
      "- Statistical recomputation: `NO`",
      paste0("- Export mode: ", if (release_mode) "RELEASE_FIGURE; independently rendered 300-ppi preview, true 600-ppi PNG, true 600-ppi LZW TIFF and exact-size PDF." else "CANDIDATE; v07 2008 x 2409 px, 300-dpi RGB QC PNG only."),
      "- User-review repair: the grey Other point/text pair is separated from the terminal threshold digit while the key remains centred beneath F.",
      "- Header route: one shared title-row grob; tag right edge fixed 2.0 mm left of title, with identical font metrics and a common vertical anchor."
    )
    writeLines(provenance, file.path(manifest_dir, "RENDERER_PROVENANCE.md"), useBytes = TRUE)
    writeLines(capture.output(sessionInfo()), file.path(manifest_dir, "R_SESSION_INFO.txt"), useBytes = TRUE)
  }, warning = record_warning),
  error = function(e) {
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
    warning_count <- write_warning_log()
    writeLines(c(
      "# Figure 2 render failure",
      "",
      paste0("- Error: `", conditionMessage(e), "`"),
      paste0("- Captured warnings: ", warning_count),
      "- Candidate status: `FAIL_NOT_AUTHORITATIVE`"
    ), failure_path, useBytes = TRUE)
    stop(e)
  }
)

warning_count <- write_warning_log()
reconciliation <- c(
  "# Full Runtime Warning Reconciliation Card",
  "",
  "- Figure: Figure 2",
  paste0("- Artifact: ", if (release_mode) "v07 final release bundle" else "v07 QC PNG"),
  "- Capture started before data-to-plot preparation: YES",
  "- Plot construction covered: YES",
  "- Grob/gtable construction covered: YES",
  "- Font-metric calculation covered: YES",
  "- Composite assembly covered: YES",
  "- All export devices covered: YES",
  paste0("- Captured warning count: ", warning_count),
  paste0("- RUNTIME_WARNING_LOG.csv rows: ", warning_count),
  paste0("- Claim allowed: ", if (warning_count == 0L) "ZERO_RUNTIME_WARNINGS" else "WARNINGS_REQUIRE_ADJUDICATION"),
  paste0("- Final Gate: ", if (warning_count == 0L) "PASS" else "HOLD")
)
writeLines(reconciliation, file.path(qc_dir, "FULL_RUNTIME_WARNING_RECONCILIATION_CARD.md"), useBytes = TRUE)
cat("PASS: Figure 2 v07 ", if (release_mode) "release bundle" else "QC candidate",
    " rendered; captured warnings = ", warning_count, "\n", sep = "")
