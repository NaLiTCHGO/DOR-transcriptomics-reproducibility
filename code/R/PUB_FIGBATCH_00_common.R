options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(ragg)
  library(gtable)
})

parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    out[[gsub("-", "_", key)]] <- value
  }
  out
}

require_cli_path <- function(x, name) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument --", gsub("_", "-", name), "=<path>")
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

cohort_order <- c("GSE274832", "GSE193136", "GSE232306")
cohort_short <- c(GSE274832 = "274832", GSE193136 = "193136", GSE232306 = "232306")
phenotype_palette <- c(DOR = "#D55E00", NOR = "#0072B2")
collection_palette <- c(HALLMARK = "#0072B2", REACTOME = "#E69F00")

activate_skill3_theme <- function(skill_root) {
  theme_path <- file.path(skill_root, "R", "theme_research_publication.R")
  if (!file.exists(theme_path)) stop("Missing active Skill 3 theme: ", theme_path)
  source(theme_path, local = FALSE)
  invisible(theme_path)
}

clean_pathway <- function(x) {
  x <- sub("^(HALLMARK|REACTOME)_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

short_sample <- function(x) sub("^SRR", "", x)

theme_pub <- function(base_size = 7.5) {
  base_theme <- if (exists("theme_research_publication", mode = "function")) {
    theme_research_publication(base_size = base_size, base_family = "Arial")
  } else {
    theme_minimal(base_size = base_size, base_family = "Arial")
  }
  base_theme +
    theme(
      text = element_text(colour = "#202020"),
      axis.title = element_text(size = 7.5),
      axis.text = element_text(size = 7.0, colour = "#202020"),
      strip.text = element_text(size = 7.2, face = "bold", colour = "#202020"),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7.0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.22, colour = "#E9E9E9"),
      plot.margin = margin(2.5, 2.5, 2.5, 2.5, unit = "mm")
    )
}

as_body_gtable <- function(body) {
  if (inherits(body, "patchwork")) return(patchwork::patchworkGrob(body))
  if (inherits(body, "ggplot")) return(ggplot2::ggplotGrob(body))
  if (inherits(body, "gtable")) return(body)
  stop("Unsupported panel body class: ", paste(class(body), collapse = "/"))
}

scientific_body_columns <- function(gt) {
  is_panel <- grepl("panel", gt$layout$name, ignore.case = TRUE) &
    !grepl("panel-area|panel-guide|panel-strip", gt$layout$name, ignore.case = TRUE)
  cells <- gt$layout[is_panel, , drop = FALSE]
  if (!nrow(cells)) {
    stop(
      "Could not locate a rendered scientific-body panel cell; layout names: ",
      paste(unique(gt$layout$name), collapse = "; ")
    )
  }
  c(left = min(cells$l), right = max(cells$r))
}

panel_block <- function(body, tag, title, subtitle = NULL, header_height = NULL,
                        body_anchor_fraction = NULL) {
  # CONTROLLED_LAYOUT_REPAIR: construct each header from the rendered panel
  # columns. Title/subtitle x=0 is therefore the scientific-body left edge,
  # independent of y-axis text, heatmap row labels or outer plot margins.
  gt <- as_body_gtable(body)
  body_cols <- scientific_body_columns(gt)
  panel_left <- as.integer(body_cols[["left"]])
  panel_right <- as.integer(body_cols[["right"]])
  if (panel_left <= 1L) {
    gt <- gtable::gtable_add_cols(gt, grid::unit(6.0, "mm"), pos = 0)
    panel_left <- panel_left + 1L
  }

  has_subtitle <- !is.null(subtitle) && nzchar(subtitle)
  header_rows_mm <- if (has_subtitle) c(4.6, 3.5) else 4.8
  gt <- gtable::gtable_add_rows(gt, grid::unit(header_rows_mm, "mm"), pos = 0)

  if (is.null(body_anchor_fraction)) {
    body_anchor_fraction <- attr(body, "scientific_body_anchor_fraction", exact = TRUE)
  }
  if (is.null(body_anchor_fraction)) body_anchor_fraction <- 0
  if (!is.numeric(body_anchor_fraction) || length(body_anchor_fraction) != 1L ||
      !is.finite(body_anchor_fraction) || body_anchor_fraction < 0 || body_anchor_fraction >= 1) {
    stop("Invalid scientific-body anchor fraction for panel ", tag)
  }

  if (body_anchor_fraction > 0) {
    # Raster-backed heatmaps may contain a dendrogram/outer image region inside
    # the ggplot panel. Anchor the header to the locked matrix border instead.
    header_l <- panel_left
    header_r <- panel_right
    title_x <- grid::unit(body_anchor_fraction, "npc")
    tag_x <- title_x - grid::unit(2.0, "mm")
  } else {
    header_l <- panel_left
    header_r <- ncol(gt)
    title_x <- grid::unit(0, "npc")
    tag_x <- grid::unit(1, "npc") - grid::unit(2.0, "mm")
  }

  grid_family <- "sans"
  tag_grob <- grid::textGrob(
    tag,
    x = tag_x,
    y = grid::unit(0.50, "npc"),
    just = c("right", "centre"),
    gp = grid::gpar(fontfamily = grid_family, fontface = "bold", fontsize = 9.0, col = "#202020")
  )
  title_grob <- grid::textGrob(
    title,
    x = title_x,
    y = grid::unit(0.50, "npc"),
    just = c("left", "centre"),
    gp = grid::gpar(fontfamily = grid_family, fontface = "bold", fontsize = 8.5, col = "#202020")
  )
  if (body_anchor_fraction > 0) {
    gt <- gtable::gtable_add_grob(
      gt, tag_grob, t = 1, l = header_l, b = 1, r = header_r,
      z = Inf, clip = "off", name = paste0("panel_tag_", tag)
    )
  } else {
    gt <- gtable::gtable_add_grob(
      gt, tag_grob, t = 1, l = 1, b = 1, r = panel_left - 1L,
      z = Inf, clip = "off", name = paste0("panel_tag_", tag)
    )
  }
  gt <- gtable::gtable_add_grob(
    gt, title_grob, t = 1, l = header_l, b = 1, r = header_r,
    z = Inf, clip = "off", name = paste0("panel_title_", tag)
  )
  if (has_subtitle) {
    subtitle_grob <- grid::textGrob(
      subtitle,
      x = title_x,
      y = grid::unit(0.52, "npc"),
      just = c("left", "centre"),
      gp = grid::gpar(fontfamily = grid_family, fontsize = 7.2, col = "#404040")
    )
    gt <- gtable::gtable_add_grob(
      gt, subtitle_grob, t = 2, l = header_l, b = 2, r = header_r,
      z = Inf, clip = "off", name = paste0("panel_subtitle_", tag)
    )
  }
  patchwork::wrap_elements(full = gt)
}

save_candidate <- function(path, width_px, height_px, build_fn) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  agg_png(path, width = width_px, height = height_px, units = "px", res = 300,
          background = "white", scaling = 1)
  ok <- FALSE
  on.exit({
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  candidate <- build_fn()
  print(candidate)
  invisible(dev.off())
  ok <- file.exists(path) && file.info(path)$size > 0
  if (!ok) stop("Candidate was not written: ", path)
  invisible(path)
}

save_release_bundle <- function(folder, figure_id, width_px_300, height_px_300, build_fn) {
  # RELEASE_FIGURE: render every final format independently from the frozen
  # R object/code path. No raster upscaling or DPI-metadata-only conversion is
  # permitted. Physical dimensions remain identical to the approved 300-dpi
  # candidate: width_px_300/300 by height_px_300/300 inches.
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    preview = file.path(folder, paste0(figure_id, "_preview_300dpi.png")),
    png600 = file.path(folder, paste0(figure_id, "_submission_600dpi.png")),
    tiff600 = file.path(folder, paste0(figure_id, "_submission_600dpi.tiff")),
    pdf = file.path(folder, paste0(figure_id, "_submission.pdf"))
  )
  width_in <- width_px_300 / 300
  height_in <- height_px_300 / 300
  # Cairo on Windows serializes the MediaBox in whole PostScript points. Use
  # the nearest whole-point dimensions so the final page is the closest exact
  # representation of the approved physical size (170.04 mm rather than the
  # truncated 169.69 mm for a nominal 170-mm page).
  pdf_width_pt <- round(width_in * 72)
  pdf_height_pt <- round(height_in * 72)

  render_once <- function(path, open_device) {
    open_device()
    device_open <- TRUE
    on.exit({
      if (isTRUE(device_open) && grDevices::dev.cur() > 1L) {
        try(grDevices::dev.off(), silent = TRUE)
      }
    }, add = TRUE)
    figure_object <- build_fn()
    print(figure_object)
    invisible(grDevices::dev.off())
    device_open <- FALSE
    if (!file.exists(path) || file.info(path)$size <= 0) {
      stop("Release artifact was not written: ", path)
    }
    invisible(path)
  }

  render_once(paths[["preview"]], function() {
    ragg::agg_png(paths[["preview"]], width = width_px_300, height = height_px_300,
                  units = "px", res = 300, background = "white", scaling = 1)
  })
  render_once(paths[["png600"]], function() {
    ragg::agg_png(paths[["png600"]], width = width_px_300 * 2L, height = height_px_300 * 2L,
                  units = "px", res = 600, background = "white", scaling = 1)
  })
  render_once(paths[["tiff600"]], function() {
    ragg::agg_tiff(paths[["tiff600"]], width = width_px_300 * 2L, height = height_px_300 * 2L,
                   units = "px", res = 600, background = "white", scaling = 1,
                   compression = "lzw", bitsize = 8)
  })
  render_once(paths[["pdf"]], function() {
    grDevices::cairo_pdf(paths[["pdf"]], width = pdf_width_pt / 72, height = pdf_height_pt / 72,
                         family = "Arial", onefile = FALSE, bg = "white")
  })

  data.frame(
    figure_id = figure_id,
    physical_width_mm = width_in * 25.4,
    physical_height_mm = height_in * 25.4,
    pdf_width_pt = pdf_width_pt,
    pdf_height_pt = pdf_height_pt,
    preview_width_px = width_px_300,
    preview_height_px = height_px_300,
    submission_width_px = width_px_300 * 2L,
    submission_height_px = height_px_300 * 2L,
    preview_png = paths[["preview"]],
    submission_png = paths[["png600"]],
    submission_tiff = paths[["tiff600"]],
    submission_pdf = paths[["pdf"]],
    stringsAsFactors = FALSE
  )
}

capture_warnings <- function(expr, warning_csv, reconciliation_md, label) {
  warnings <- character()
  error_message <- NULL
  error_calls <- character()
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      error_calls <<- vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), character(1))
    }
  )
  dir.create(dirname(warning_csv), recursive = TRUE, showWarnings = FALSE)
  warning_df <- data.frame(
    warning_id = seq_along(warnings),
    message = warnings,
    stringsAsFactors = FALSE
  )
  write_csv(warning_df, warning_csv)
  card <- c(
    "# Full Runtime Warning Reconciliation Card",
    "",
    paste0("- Render group: ", label),
    "- Capture covers data-to-plot preparation, grob construction, font metrics, composite assembly and export: YES",
    paste0("- Captured warning count: ", length(warnings)),
    paste0("- Error: ", if (is.null(error_message)) "NONE" else paste0("`", error_message, "`")),
    paste0("- Final Gate: ", if (length(warnings) == 0L && is.null(error_message)) "PASS" else "HOLD"),
    if (length(error_calls)) c("", "## Error call stack", "", paste0("- `", error_calls, "`")) else NULL
  )
  writeLines(card, reconciliation_md, useBytes = TRUE)
  if (!is.null(error_message)) stop(error_message)
  if (length(warnings) > 0L) stop(label, " produced ", length(warnings), " runtime warnings; see ", warning_csv)
  invisible(TRUE)
}
