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
source(file.path(project_root, "04_code/R/PUB_FIGBATCH_00_common.R"), local = FALSE)
theme_path <- activate_skill3_theme(skill_root)
figure_set <- unique(strsplit(toupper(get_optional_arg("figure-set", "ALL")), ",", fixed = TRUE)[[1]])
figure_set <- trimws(figure_set)
release_mode <- identical(toupper(get_optional_arg("release-mode", "FALSE")), "TRUE")
allowed_figure_sets <- c("ALL", "FIG3", "FIG4", "FIG5")
if (!length(figure_set) || any(!figure_set %in% allowed_figure_sets)) {
  stop("Unsupported --figure-set: ", paste(figure_set, collapse = ","))
}
render_requested <- function(id) "ALL" %in% figure_set || id %in% figure_set

source_root <- file.path(run_root, "source_data")
figure_root <- file.path(run_root, "figures", "main")
qc_root <- file.path(run_root, "qc")
manifest_root <- file.path(run_root, "manifests")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)

read_src <- function(group, name) read_csv(file.path(source_root, group, name), show_col_types = FALSE, progress = FALSE)

build_figure3 <- function() {
  meta <- read_src("FIG3_S4", "M04_GENE_LEVEL_META_ANALYSIS.csv")
  pairwise <- read_src("FIG3_S4", "M04_PAIRWISE_REPRODUCIBILITY.csv")

  matrix_rows <- list()
  k <- 0L
  for (metric in c("spearman_rho_log2fc", "direction_concordance_all_nonzero")) {
    for (a in cohort_order) for (b in cohort_order) {
      value <- if (a == b) {
        1
      } else {
        hit <- pairwise[(pairwise$cohort_a == a & pairwise$cohort_b == b) |
                          (pairwise$cohort_a == b & pairwise$cohort_b == a), , drop = FALSE]
        if (nrow(hit) != 1L) stop("Pairwise matrix lookup failed: ", a, " / ", b)
        hit[[metric]][[1]]
      }
      k <- k + 1L
      matrix_rows[[k]] <- data.frame(metric = metric, row = a, col = b, value = value)
    }
  }
  matrix_long <- bind_rows(matrix_rows) %>%
    mutate(
      row = factor(row, levels = rev(cohort_order), labels = rev(unname(cohort_short[cohort_order]))),
      col = factor(col, levels = cohort_order, labels = unname(cohort_short[cohort_order]))
    )

  p_rho <- ggplot(filter(matrix_long, metric == "spearman_rho_log2fc"), aes(col, row, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.3f", value)), family = "Arial", size = 6.8 / ggplot2::.pt) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), name = expression(rho)) +
    labs(x = NULL, y = NULL, title = "Spearman rank correlation") +
    coord_equal() + theme_pub() +
    guides(fill = guide_colourbar(barwidth = grid::unit(2.3, "mm"), barheight = grid::unit(19, "mm"),
                                  title.position = "top")) +
    theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel", panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "right",
          legend.box.spacing = grid::unit(0.4, "mm"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(2.5, 0.7, 2.5, 2.5, unit = "mm"))
  p_dir <- ggplot(filter(matrix_long, metric == "direction_concordance_all_nonzero"), aes(col, row, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = percent(value, accuracy = 0.1)), family = "Arial", size = 6.8 / ggplot2::.pt) +
    scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1), labels = percent,
                        name = "Direction") +
    labs(x = NULL, y = NULL, title = "Effect-direction concordance") +
    coord_equal() + theme_pub() +
    guides(fill = guide_colourbar(barwidth = grid::unit(2.3, "mm"), barheight = grid::unit(19, "mm"),
                                  title.position = "top")) +
    theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel", panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "right",
          legend.box.spacing = grid::unit(0.4, "mm"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(2.5, 0.7, 2.5, 2.5, unit = "mm"))
  # Fixed-aspect matrices otherwise inherit more horizontal room than they can
  # use, leaving an artificial blank band before each right-side colour key.
  # A larger trailing spacer narrows both complete matrix+key children without
  # changing any font, scale, cell or legend dimension.
  panel_a_body <- p_rho + p_dir + plot_spacer() + plot_layout(widths = c(1, 1, 0.40))

  pair_rows <- list()
  for (i in seq_len(nrow(pairwise))) {
    a <- pairwise$cohort_a[[i]]
    b <- pairwise$cohort_b[[i]]
    pair_rows[[i]] <- data.frame(
      x = meta[[paste0("beta_", a)]], y = meta[[paste0("beta_", b)]],
      pair = paste(cohort_short[[a]], "vs", cohort_short[[b]]), stringsAsFactors = FALSE
    )
  }
  scatter <- bind_rows(pair_rows)
  scatter$pair <- factor(scatter$pair, levels = vapply(seq_len(nrow(pairwise)), function(i) {
    paste(cohort_short[[pairwise$cohort_a[[i]]]], "vs", cohort_short[[pairwise$cohort_b[[i]]]])
  }, character(1)))
  stats <- pairwise %>%
    mutate(pair = factor(paste(cohort_short[cohort_a], "vs", cohort_short[cohort_b]), levels = levels(scatter$pair)),
           label = sprintf("rho = %.3f; direction = %.1f%%", spearman_rho_log2fc, 100 * direction_concordance_all_nonzero),
           x = -Inf, y = Inf)
  lim <- as.numeric(quantile(abs(c(scatter$x, scatter$y)), 0.99, na.rm = TRUE))
  p_scatter <- ggplot(scatter, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "#D0D0D0", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "#D0D0D0", linewidth = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#808080", linewidth = 0.3) +
    geom_point(size = 0.32, alpha = 0.16, colour = "#334155", stroke = 0) +
    geom_text(data = stats, aes(x, y, label = label), inherit.aes = FALSE,
              hjust = -0.02, vjust = 1.15, family = "Arial", size = 6.8 / ggplot2::.pt) +
    facet_wrap(~pair, nrow = 1) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    labs(x = "Cohort A DOR-minus-NOR log2FC", y = "Cohort B DOR-minus-NOR log2FC") +
    theme_pub() + theme(panel.grid.minor = element_blank())

  order_i2 <- c("LOW_LT25", "MODERATE_25_TO_LT50", "SUBSTANTIAL_50_TO_LT75", "HIGH_GE75")
  count_i2 <- meta %>% count(heterogeneity_category) %>%
    mutate(heterogeneity_category = factor(heterogeneity_category, levels = order_i2,
                                           labels = c("<25", "25-<50", "50-<75", ">=75")))
  p_i2_hist <- ggplot(meta, aes(i2_percent)) +
    geom_histogram(breaks = seq(0, 100, by = 5), fill = "#56B4E9", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 75, linetype = "dashed", colour = "#D55E00", linewidth = 0.4) +
    labs(x = expression(I^2~"(%)"), y = "Genes", title = "Distribution") + theme_pub() +
    theme(plot.title = element_text(size = 7.5, face = "bold"))
  p_i2_count <- ggplot(count_i2, aes(heterogeneity_category, n, fill = heterogeneity_category)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = comma(n)), vjust = -0.25, family = "Arial", size = 6.7 / ggplot2::.pt) +
    scale_fill_manual(values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = expression(I^2~"category (%)"), y = "Genes", title = "Categories") + theme_pub() +
    theme(plot.title = element_text(size = 7.5, face = "bold"),
          axis.text.x = element_text(size = 6.2, angle = 18, hjust = 1))
  panel_c_body <- p_i2_hist + p_i2_count

  strict <- meta %>% filter(strict_exploratory_consensus) %>% arrange(random_padj, i2_percent)
  if (nrow(strict) != 2L) stop("Expected two strict exploratory consensus genes; found ", nrow(strict))
  effect_long <- bind_rows(lapply(cohort_order, function(c) {
    data.frame(gene_name = strict$gene_name, cohort = c, log2FC = strict[[paste0("beta_", c)]],
               i2_percent = strict$i2_percent)
  })) %>%
    mutate(
      cohort = factor(cohort, levels = cohort_order, labels = unname(cohort_short[cohort_order])),
      gene_label = gene_name
    )
  gene_levels <- unique(effect_long$gene_label)
  effect_long$gene_label <- factor(effect_long$gene_label, levels = rev(gene_levels))
  fill_lim <- max(abs(effect_long$log2FC))
  p_effect <- ggplot(effect_long, aes(cohort, gene_label, fill = log2FC)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", log2FC)), family = "Arial", size = 7 / ggplot2::.pt) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-fill_lim, fill_lim), name = "log2FC") +
    labs(x = "Cohort", y = NULL) + theme_pub() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(size = 6.6, angle = 45, hjust = 1, vjust = 1))

  a <- panel_block(panel_a_body, "A", "Three-cohort gene-effect reproducibility matrices",
                   "Diagonal values are self-comparisons; off-diagonals use 13,993 shared genes",
                   body_anchor_fraction = 0.076)
  b <- panel_block(p_scatter, "B", "Pairwise gene-effect correspondence",
                   "Display clipped at the pooled 99th percentile of absolute log2FC")
  c <- panel_block(panel_c_body, "C", "Gene-level heterogeneity",
                   "I2 is descriptive with only three cohorts; 5,883 genes are at or above 75%")
  d <- panel_block(p_effect, "D", "Exploratory consensus effects",
                   "Exploratory only; not biomarkers")
  # Centre the three-facet B module at 9/11 of the page width. This removes its
  # previous visual dominance while retaining all point and text sizes. The
  # 6:5 bottom split preserves the established C:D width balance.
  wrap_plots(list(a, b, c, d), design = "AAAAAAAAAAA\n#BBBBBBBBB#\nCCCCCCDDDDD",
             heights = c(0.78, 1.0, 0.92)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

build_figure4 <- function() {
  convergence <- read_src("FIG4_S5", "M05_PATHWAY_CONVERGENCE.csv")
  pairwise <- read_src("FIG4_S5", "M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv")
  strict <- read_src("FIG4_S5", "M05_STRICT_CONVERGENT_PATHWAYS.csv")
  fgsea <- read_src("FIG4_S5", "M05_FGSEA_ALL_RESULTS.csv")

  repro <- bind_rows(
    pairwise %>% transmute(collection, cohort_a, cohort_b, level = "Gene log2FC", rho = m04_gene_log2fc_spearman_rho),
    pairwise %>% transmute(collection, cohort_a, cohort_b, level = "Pathway NES", rho = pathway_nes_spearman_rho)
  ) %>% mutate(pair = paste(cohort_short[cohort_a], "vs", cohort_short[cohort_b]),
               level = factor(level, levels = c("Gene log2FC", "Pathway NES")))
  p_repro <- ggplot(repro, aes(pair, rho, fill = level)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#777777") +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    facet_wrap(~collection, nrow = 1) +
    scale_fill_manual(values = c("Gene log2FC" = "#BDBDBD", "Pathway NES" = "#0072B2"), name = NULL) +
    coord_cartesian(ylim = c(-1, 1)) + labs(x = NULL, y = "Spearman rho") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 22, hjust = 1))

  strict_h <- strict %>% filter(collection == "HALLMARK") %>% arrange(median_nes)
  if (nrow(strict_h) != 8L) stop("Expected eight strict Hallmark pathways")
  h_long <- bind_rows(lapply(cohort_order, function(c) {
    data.frame(pathway = strict_h$pathway, cohort = c, NES = strict_h[[paste0("NES_", c)]])
  })) %>%
    mutate(pathway_label = factor(clean_pathway(pathway), levels = clean_pathway(strict_h$pathway)),
           cohort = factor(cohort, levels = cohort_order, labels = unname(cohort_short[cohort_order])))
  lim_h <- max(abs(h_long$NES))
  p_h <- ggplot(h_long, aes(cohort, pathway_label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", NES)), family = "Arial", size = 6.5 / ggplot2::.pt) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim_h, lim_h), name = "NES") +
    labs(x = "Cohort", y = NULL) + theme_pub() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  counts <- convergence %>% group_by(collection) %>% summarise(
    `All-three same direction` = sum(all_three_same_nes_direction),
    `Direction + >=1 cohort FDR` = sum(directional_supported_pathway),
    `Strict convergence` = sum(strict_pathway_convergence), .groups = "drop"
  )
  counts_long <- bind_rows(lapply(names(counts)[-1], function(nm) {
    data.frame(collection = counts$collection, category = nm, pathway_n = counts[[nm]])
  }))
  counts_long$category <- factor(counts_long$category,
                                 levels = c("All-three same direction", "Direction + >=1 cohort FDR", "Strict convergence"))
  p_counts <- ggplot(counts_long, aes(category, pathway_n, fill = collection)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_text(aes(label = pathway_n), position = position_dodge(width = 0.75), vjust = -0.25,
              family = "Arial", size = 6.5 / ggplot2::.pt) +
    scale_fill_manual(values = collection_palette, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = NULL, y = "Pathways") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, hjust = 1))

  scatter_parts <- list(); idx <- 0L
  for (collection in c("HALLMARK", "REACTOME")) {
    sub <- fgsea %>% filter(.data$collection == collection)
    for (i in seq_len(length(cohort_order) - 1L)) for (j in (i + 1L):length(cohort_order)) {
      a <- cohort_order[[i]]; b <- cohort_order[[j]]
      xa <- sub %>% filter(cohort == a) %>% select(pathway, NES_a = NES)
      xb <- sub %>% filter(cohort == b) %>% select(pathway, NES_b = NES)
      idx <- idx + 1L
      scatter_parts[[idx]] <- inner_join(xa, xb, by = "pathway") %>%
        mutate(collection = collection, pair = paste(cohort_short[[a]], "vs", cohort_short[[b]]))
    }
  }
  nes_scatter <- bind_rows(scatter_parts)
  p_scatter <- ggplot(nes_scatter, aes(NES_a, NES_b)) +
    geom_hline(yintercept = 0, colour = "#D4D4D4", linewidth = 0.22) +
    geom_vline(xintercept = 0, colour = "#D4D4D4", linewidth = 0.22) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#888888", linewidth = 0.28) +
    geom_point(size = 0.34, alpha = 0.28, colour = "#334155", stroke = 0) +
    facet_grid(collection ~ pair, scales = "free") +
    labs(x = "Cohort A NES", y = "Cohort B NES") + theme_pub(7.2)

  a <- panel_block(p_repro, "A", "Pathway versus gene rank reproducibility",
                   "Improvement is comparison-specific")
  b <- panel_block(p_h, "B", "Strict Hallmark convergence",
                   "Negative NES = DOR-down enrichment")
  c <- panel_block(p_counts, "C", "Pathway convergence counts",
                   "Strict rule: same NES direction, FDR < 0.05 in at least two cohorts and median |NES| at least 1")
  d <- panel_block(p_scatter, "D", "Pairwise pathway NES correspondence",
                   "All 49 Hallmark and 1,016 Reactome pathways are shown")
  wrap_plots(list(a, c, b, d), design = "AC\nBB\nDD",
             heights = c(0.78, 1.0, 1.08), widths = c(1, 1)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

build_figure5 <- function() {
  summary <- read_src("FIG5_S6", "M06_COLLECTION_LOCO_SUMMARY.csv") %>%
    mutate(held_out_cohort = factor(held_out_cohort, levels = cohort_order,
                                    labels = paste("Hold out", unname(cohort_short[cohort_order]))))
  core <- read_src("FIG5_S6", "M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv") %>% filter(universal_loco_core)
  frozen <- read_src("FIG5_S6", "M06_FROZEN_M05_STRICT_SURVIVAL.csv") %>% filter(collection == "HALLMARK")
  if (nrow(core) != 8L || nrow(frozen) != 8L) stop("M06 core/frozen source contract failed")

  funnel <- bind_rows(
    summary %>% transmute(collection, held_out_cohort, metric = "Pair-selected", n = retained_pair_selected_n),
    summary %>% transmute(collection, held_out_cohort, metric = "Held-out direction", n = held_out_direction_replication_n),
    summary %>% transmute(collection, held_out_cohort, metric = "Held-out strict", n = held_out_strict_replication_n)
  )
  funnel$metric <- factor(funnel$metric, levels = c("Pair-selected", "Held-out direction", "Held-out strict"))
  p_funnel <- ggplot(funnel, aes(held_out_cohort, n, fill = metric)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7) +
    geom_text(aes(label = n), position = position_dodge(width = 0.78), vjust = -0.2,
              family = "Arial", size = 5.5 / ggplot2::.pt) +
    facet_wrap(~collection, scales = "free_y") +
    scale_fill_manual(values = c("#BDBDBD", "#56B4E9", "#D55E00"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "Pathways") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, hjust = 1))

  rates <- bind_rows(
    summary %>% transmute(collection, held_out_cohort, metric = "Direction", rate = direction_replication_rate_among_selected),
    summary %>% transmute(collection, held_out_cohort, metric = "Strict", rate = strict_replication_rate_among_selected)
  )
  rates$metric <- factor(rates$metric, levels = c("Direction", "Strict"))
  p_rates <- ggplot(rates, aes(held_out_cohort, rate, fill = metric)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    geom_text(aes(y = rate + ifelse(metric == "Strict", 0.075, 0.025),
                  label = percent(rate, accuracy = 1)),
              position = position_dodge(width = 0.72), vjust = 0.5, colour = "#202020",
              family = "Arial", size = 5.4 / ggplot2::.pt) +
    facet_wrap(~collection) +
    scale_fill_manual(values = c(Direction = "#777777", Strict = "#0072B2"), name = NULL) +
    scale_y_continuous(limits = c(0, 1.12), labels = percent) +
    labs(x = NULL, y = "Replication among pair-selected") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, hjust = 1))

  core <- core %>% arrange(collection, pathway)
  display_core_pathway <- function(x) {
    cleaned <- clean_pathway(x)
    replacements <- c(
      "Auf1 Hnrnp D0 Binds and Destabilizes Mrna" = "AUF1/HNRNPD destabilizes mRNA",
      "Eukaryotic Translation Initiation" = "Translation initiation",
      "Nervous System Development" = "Nervous-system development",
      "Neutrophil Degranulation" = "Neutrophil degranulation",
      "Regulation of Ras by Gaps" = "RAS regulation by GAPs",
      "Rna Processing" = "RNA processing",
      "Transcriptional Regulation by Runx2" = "RUNX2 transcriptional regulation"
    )
    hit <- match(cleaned, names(replacements))
    cleaned[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
    cleaned
  }
  core_long <- bind_rows(lapply(cohort_order, function(c) {
    data.frame(collection = core$collection, pathway = core$pathway, cohort = c, NES = core[[paste0("NES_", c)]])
  })) %>% mutate(
    label = paste0(ifelse(collection == "HALLMARK", "H: ", "R: "), display_core_pathway(pathway)),
    cohort = factor(cohort, levels = cohort_order, labels = unname(cohort_short[cohort_order]))
  )
  label_order <- unique(core_long$label)
  core_long$label <- factor(core_long$label, levels = rev(label_order))
  lim_core <- max(abs(core_long$NES))
  p_core <- ggplot(core_long, aes(cohort, label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", NES)), family = "Arial", size = 6.2 / ggplot2::.pt) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim_core, lim_core), name = "NES") +
    labs(x = "Cohort", y = NULL) + theme_pub() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.2),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  survival <- bind_rows(lapply(cohort_order, function(c) {
    data.frame(pathway = frozen$pathway, held_out = c,
               survives = as.integer(frozen[[paste0("survives_without_", c)]]))
  })) %>% mutate(
    label = clean_pathway(pathway),
    held_out = factor(held_out, levels = cohort_order, labels = paste("Without", unname(cohort_short[cohort_order])))
  )
  frozen_order <- frozen %>% arrange(desc(surviving_rotation_n), pathway) %>% pull(pathway)
  survival$label <- factor(survival$label, levels = rev(clean_pathway(frozen_order)))
  p_survival <- ggplot(survival, aes(held_out, label, fill = factor(survives))) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = ifelse(survives == 1, "PASS", "-")), family = "Arial", size = 6.3 / ggplot2::.pt) +
    scale_fill_manual(values = c(`0` = "#F0F0F0", `1` = "#009E73"), guide = "none") +
    labs(x = NULL, y = NULL) + theme_pub() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 6.2))

  a <- panel_block(p_funnel, "A", "LOCO selection and replication",
                   "The held-out cohort is excluded from selection")
  b <- panel_block(p_rates, "B", "Held-out replication rates",
                   "Conditional on retained-pair selection")
  c <- panel_block(p_core, "C", "Universal internal LOCO core",
                   "1 Hallmark + 7 Reactome pass all rotations")
  d <- panel_block(p_survival, "D", "Frozen Hallmark survival",
                   "Only P53 survives all rotations")
  # Allocate one extra twentieth of the lower row to D. This moves D left and
  # enlarges its scientific body, correcting the right-heavy lower-row centre
  # without changing labels, values, order or typography.
  wrap_plots(list(a, b, c, d),
             design = "AAAAAAAAAABBBBBBBBBB\nCCCCCCCCCDDDDDDDDDDD",
             heights = c(1, 1.05)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

warning_csv <- file.path(qc_root, "RUNTIME_WARNING_LOG_MAIN.csv")
reconciliation <- file.path(qc_root, "FULL_RUNTIME_WARNING_RECONCILIATION_MAIN.md")
release_rows <- list()
capture_warnings({
  if (render_requested("FIG3")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_3"), "Figure_3", 2008, 2409, build_figure3)
    } else {
      save_candidate(file.path(figure_root, "Figure_3_Gene_Reproducibility_and_Heterogeneity_v09_QC.png"),
                     2008, 2409, build_figure3)
    }
  }
  if (render_requested("FIG4")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_4"), "Figure_4", 2008, 2600, build_figure4)
    } else {
      save_candidate(file.path(figure_root, "Figure_4_Pathway_Convergence_v05_QC.png"),
                     2008, 2600, build_figure4)
    }
  }
  if (render_requested("FIG5")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_5"), "Figure_5", 2008, 2200, build_figure5)
    } else {
      save_candidate(file.path(figure_root, "Figure_5_LOCO_Robustness_v08_QC.png"),
                     2008, 2200, build_figure5)
    }
  }
}, warning_csv, reconciliation, paste0("MAIN_FIGURES_", paste(figure_set, collapse = "_")))
if (release_mode) {
  write_csv(bind_rows(release_rows), file.path(manifest_root, "FINAL_EXPORT_MANIFEST_MAIN.csv"))
}

provenance <- c(
  "# Main-result renderer provenance",
  "",
  paste0("- Active Skill root: `", skill_root, "`"),
  paste0("- Theme module: `", theme_path, "`"),
  "- Active suite: `research_scientific_figure_suite v2.2-beta.9`",
  paste0("- R: `", R.version.string, "`"),
  paste0("- ggplot2: `", as.character(packageVersion("ggplot2")), "`"),
  paste0("- patchwork: `", as.character(packageVersion("patchwork")), "`"),
  paste0("- ragg: `", as.character(packageVersion("ragg")), "`"),
  "- Scientific/statistical recomputation: `NO`",
  "- Optional Figure 4C wording refresh: the support text now states the exact FDR < 0.05 criterion used in the locked strict-convergence definition.",
  paste0("- Figure-set route: `", paste(figure_set, collapse = ","), "`"),
  "- Figure 3 v09: A child widths tightened to bring fixed-aspect matrix colour keys closer; B centred at 9/11 page width with typography unchanged.",
  "- Figure 5 v08: lower row reallocated 9:11 so D moves left and widens; display labels, values and typography unchanged.",
  paste0("- Output mode: `", if (release_mode) "RELEASE_FIGURE" else "CANDIDATE_300DPI", "`"),
  if (release_mode) {
    "- Outputs: independently rendered 300-dpi preview, true 600-ppi PNG, LZW TIFF and physical-size Cairo PDF."
  } else {
    "- Outputs: selected 300-dpi RGB QC PNG candidates only; no final-release formats"
  }
)
writeLines(provenance, file.path(manifest_root, "RENDERER_PROVENANCE_MAIN.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(manifest_root, "R_SESSION_INFO_MAIN.txt"), useBytes = TRUE)
cat("PASS: selected main figures rendered with zero captured warnings; route=",
    paste(figure_set, collapse = ","), "; mode=", if (release_mode) "RELEASE" else "CANDIDATE", "\n", sep = "")
