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
if (!length(figure_set) || any(!figure_set %in% c("ALL", "S4", "S5", "S6"))) {
  stop("Unsupported --figure-set: ", paste(figure_set, collapse = ","))
}
render_requested <- function(id) "ALL" %in% figure_set || id %in% figure_set

source_root <- file.path(run_root, "source_data")
figure_root <- file.path(run_root, "figures", "supplement")
qc_root <- file.path(run_root, "qc")
manifest_root <- file.path(run_root, "manifests")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)
read_src <- function(group, name) read_csv(file.path(source_root, group, name), show_col_types = FALSE, progress = FALSE)

build_s4 <- function() {
  overlap <- read_src("FIG3_S4", "M04_CANDIDATE_OVERLAP_SUMMARY.csv") %>%
    filter(candidate_membership != "NONE") %>% arrange(gene_n) %>%
    mutate(label = gsub("GSE", "", candidate_membership), label = gsub("\\+", " + ", label))
  loco <- read_src("FIG3_S4", "M04_LOCO_STABILITY.csv") %>%
    mutate(held_out = factor(paste("Without", cohort_short[held_out_cohort]),
                             levels = paste("Without", unname(cohort_short[cohort_order]))))

  p_overlap <- ggplot(overlap, aes(gene_n, reorder(label, gene_n))) +
    geom_segment(aes(x = 0.8, xend = gene_n, yend = reorder(label, gene_n)), colour = "#9ECAE1", linewidth = 1.1) +
    geom_point(size = 2.2, colour = "#0072B2") +
    geom_text(data = overlap %>% filter(gene_n != 9), aes(label = comma(gene_n)),
              hjust = -0.25, family = "Arial", size = 6.8 / ggplot2::.pt) +
    # Keep the single-digit count on the point's horizontal centreline and
    # explicitly nudge its x position rightward so its glyph cannot cover the
    # marker. A multiplicative nudge is required because the axis is log10.
    geom_text(data = overlap %>% filter(gene_n == 9),
              aes(x = gene_n * 1.55, label = comma(gene_n)),
              hjust = 0, vjust = 0.5, family = "Arial", size = 6.8 / ggplot2::.pt) +
    scale_x_log10(expand = expansion(mult = c(0.02, 0.22))) +
    labs(x = "Within-cohort thresholded genes (log scale)", y = "Candidate membership") + theme_pub() +
    theme(panel.grid.major.y = element_blank())

  loco_long <- bind_rows(
    loco %>% transmute(held_out, metric = "Rank correlation", value = spearman_vs_full_random_beta),
    loco %>% transmute(held_out, metric = "Direction agreement", value = direction_agreement_vs_full_random)
  )
  p_loco <- ggplot(loco_long, aes(held_out, value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(width = 0.72), vjust = -0.25,
              family = "Arial", size = 6.7 / ggplot2::.pt) +
    scale_fill_manual(values = c("Rank correlation" = "#CC79A7", "Direction agreement" = "#009E73"), name = NULL) +
    scale_y_continuous(limits = c(0, 1.08)) + labs(x = NULL, y = "Agreement with full model") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, hjust = 1))

  a <- panel_block(p_overlap, "A", "Thresholded candidate overlap",
                   "13,993 common genes; no three-cohort overlap")
  b <- panel_block(p_loco, "B", "Gene-level leave-one-cohort-out stability",
                   "Retained-pair effects versus the full model")
  a + plot_spacer() + b + plot_layout(widths = c(0.86, 0.06, 1.08)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

build_s5 <- function() {
  convergence <- read_src("FIG4_S5", "M05_PATHWAY_CONVERGENCE.csv")
  p_disp <- ggplot(convergence, aes(median_abs_nes, sd_nes, colour = strict_pathway_convergence)) +
    geom_point(alpha = 0.52, size = 0.85, stroke = 0) +
    facet_wrap(~collection, scales = "free") +
    scale_colour_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = "#D55E00"),
                        labels = c("Other", "Strict convergence"), name = NULL) +
    labs(x = "Median |NES|", y = "SD of NES across cohorts") + theme_pub() +
    theme(legend.position = "bottom", legend.box.spacing = grid::unit(0.75, "mm"),
          legend.margin = margin(0.5, 0, 0, 0, unit = "mm"),
          plot.margin = margin(2.5, 2.5, 1.5, 2.5, unit = "mm"))

  reactome_top <- convergence %>% filter(collection == "REACTOME") %>%
    arrange(desc(strict_pathway_convergence), desc(fgsea_fdr_lt_0_05_cohort_n),
            signed_stouffer_padj, desc(median_abs_nes)) %>% slice_head(n = 30)
  display_reactome_pathway <- function(x) {
    cleaned <- clean_pathway(x)
    replacements <- c(
      "Rrna Processing" = "rRNA processing",
      "Auf1 Hnrnp D0 Binds and Destabilizes Mrna" = "AUF1/HNRNPD destabilizes mRNA",
      "Eukaryotic Translation Initiation" = "Translation initiation",
      "Eukaryotic Translation Elongation" = "Translation elongation",
      "Transcriptional Regulation by Runx2" = "RUNX2 transcriptional regulation",
      "Regulation of Ras by Gaps" = "RAS regulation by GAPs",
      "Ribosome Associated Quality Control" = "Ribosome-associated quality control",
      "Nonsense Mediated Decay Nmd" = "Nonsense-mediated decay (NMD)",
      "Response of Eif2ak4 Gcn2 to Amino Acid Deficiency" = "GCN2 response to amino-acid deficiency",
      "Activation of the Mrna Upon Binding of the Cap Binding Complex and Eifs and Subsequent Binding to 43s" = "mRNA cap binding and 43S recruitment",
      "Regulation of Expression of Slits and Robos" = "SLIT/ROBO expression regulation",
      "Ribosome Quality Control Rqc Complex Extracts and Degrades Nascent Peptide" = "RQC degradation of nascent peptides",
      "Srp Dependent Cotranslational Protein Targeting to Membrane" = "SRP-dependent cotranslational targeting",
      "Sars Cov 1 Modulates Host Translation Machinery" = "SARS-CoV-1 modulation of host translation",
      "The Role of Gtse1 in G2 m Progression after G2 Checkpoint" = "GTSE1 in G2/M progression",
      "Influenza Viral Rna Transcription and Replication" = "Influenza RNA transcription/replication",
      "Regulation of Cholesterol Biosynthesis by Srebp Srebf" = "SREBP regulation of cholesterol synthesis",
      "Nuclear Events Mediated by Nfe2l2" = "NFE2L2-mediated nuclear events",
      "Ubiquitin Dependent Degradation of Cyclin d" = "Ubiquitin-dependent cyclin D degradation",
      "Activation of Gene Expression by Srebf Srebp" = "SREBP-dependent gene activation"
    )
    hit <- match(cleaned, names(replacements))
    cleaned[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
    cleaned
  }
  reactome_display_levels <- display_reactome_pathway(reactome_top$pathway)
  if (anyDuplicated(reactome_display_levels)) stop("Display-only Reactome labels are not unique")
  reactome_long <- bind_rows(lapply(cohort_order, function(c) {
    data.frame(pathway = reactome_top$pathway, cohort = c, NES = reactome_top[[paste0("NES_", c)]])
  })) %>% mutate(
    label = factor(display_reactome_pathway(pathway), levels = rev(reactome_display_levels)),
    cohort = factor(cohort, levels = cohort_order, labels = unname(cohort_short[cohort_order]))
  )
  lim <- max(abs(reactome_long$NES))
  p_heat <- ggplot(reactome_long, aes(cohort, label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = sprintf("%.2f", NES)), family = "Arial", size = 5.8 / ggplot2::.pt) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), name = "NES") +
    labs(x = "Cohort", y = NULL) + theme_pub(7.0) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.3),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  a <- panel_block(p_disp, "A", "Pathway NES magnitude and dispersion",
                   "Strict pathways are highlighted; dispersion remains substantial for many terms")
  b <- panel_block(p_heat, "B", "Reactome top-30 pathway context",
                   "Frozen ordering; Reactome is secondary context")
  wrap_plots(list(a, b), design = "A\nB", heights = c(0.58, 1.42)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

build_s6 <- function() {
  summary <- read_src("FIG5_S6", "M06_COLLECTION_LOCO_SUMMARY.csv") %>%
    mutate(held_out = factor(paste("Hold out", cohort_short[held_out_cohort]),
                             levels = paste("Hold out", unname(cohort_short[cohort_order]))))
  tiers <- read_src("FIG5_S6", "M06_ROBUSTNESS_TIER_COUNTS.csv")
  context <- bind_rows(
    summary %>% transmute(collection, held_out, metric = "Retained-pair NES rho", value = retained_pair_global_nes_rho),
    summary %>% transmute(collection, held_out, metric = "Direction concordance", value = retained_pair_global_direction_concordance)
  )
  p_context <- ggplot(context, aes(held_out, value, fill = collection)) +
    geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.25) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = collection_palette, name = NULL) +
    labs(x = NULL, y = "Global retained-pair context") + theme_pub() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, hjust = 1))

  tier_order <- c("ZERO_ROTATIONS", "ONE_ROTATION_ONLY", "TWO_ROTATIONS", "UNIVERSAL_ALL_THREE_ROTATIONS")
  grid <- expand.grid(collection = c("HALLMARK", "REACTOME"), robustness_tier = tier_order, stringsAsFactors = FALSE)
  tier_complete <- left_join(grid, tiers[, c("collection", "robustness_tier", "pathway_n")],
                             by = c("collection", "robustness_tier")) %>%
    mutate(pathway_n = ifelse(is.na(pathway_n), 0, pathway_n),
           robustness_tier = factor(robustness_tier, levels = tier_order,
                                    labels = c("0", "1", "2", "3 (universal)")))
  p_tier <- ggplot(tier_complete, aes(robustness_tier, pathway_n, fill = collection)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    geom_text(aes(label = pathway_n), position = position_dodge(width = 0.72), vjust = -0.25,
              family = "Arial", size = 6.7 / ggplot2::.pt) +
    scale_fill_manual(values = collection_palette, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = "Successful retained-pair rotations", y = "Frozen strict pathways") + theme_pub() +
    theme(legend.position = "bottom")

  a <- panel_block(p_context, "A", "Global retained-pair pathway context",
                   "Global ranks and directions are not uniformly reproducible across rotations")
  b <- panel_block(p_tier, "B", "Frozen strict-pathway stability tiers",
                   "Universal = survival in all three rotations")
  wrap_plots(list(a, b), design = "AB", widths = c(1.12, 0.88)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
}

warning_csv <- file.path(qc_root, "RUNTIME_WARNING_LOG_SUPP_SYNTHESIS.csv")
reconciliation <- file.path(qc_root, "FULL_RUNTIME_WARNING_RECONCILIATION_SUPP_SYNTHESIS.md")
release_rows <- list()
capture_warnings({
  if (render_requested("S4")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_S4"), "Figure_S4", 2008, 1150, build_s4)
    } else {
      save_candidate(file.path(figure_root, "Figure_S4_Gene_Synthesis_Sensitivity_v05_QC.png"), 2008, 1150, build_s4)
    }
  }
  if (render_requested("S5")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_S5"), "Figure_S5", 2008, 2300, build_s5)
    } else {
      save_candidate(file.path(figure_root, "Figure_S5_Pathway_Context_v07_QC.png"), 2008, 2300, build_s5)
    }
  }
  if (render_requested("S6")) {
    if (release_mode) {
      release_rows[[length(release_rows) + 1L]] <- save_release_bundle(
        file.path(figure_root, "Figure_S6"), "Figure_S6", 2008, 1150, build_s6)
    } else {
      save_candidate(file.path(figure_root, "Figure_S6_LOCO_Context_v04_QC.png"), 2008, 1150, build_s6)
    }
  }
}, warning_csv, reconciliation, "SUPPLEMENTARY_SYNTHESIS_FIGURES_S4_TO_S6")
if (release_mode) {
  write_csv(bind_rows(release_rows), file.path(manifest_root, "FINAL_EXPORT_MANIFEST_SUPP_SYNTHESIS.csv"))
}

provenance <- c(
  "# Synthesis-supplement renderer provenance",
  "",
  paste0("- Active Skill root: `", skill_root, "`"),
  paste0("- Theme module: `", theme_path, "`"),
  "- Active suite: `research_scientific_figure_suite v2.2-beta.9`",
  "- Scientific/statistical recomputation: `NO`",
  "- S4 sources: locked M04 tables; S5 sources: locked M05 tables; S6 sources: locked M06 tables",
  paste0("- Figure-set route: `", paste(figure_set, collapse = ","), "`"),
  "- S4 v05: support text shortened and a controlled inter-panel gutter added to prevent header collision.",
  "- S5 v07: legend remains close but clear of the x-axis title; display-only pathway abbreviations reduce left-label width while locked pathway IDs/order/NES remain unchanged.",
  paste0("- Output mode: `", if (release_mode) "RELEASE_FIGURE" else "CANDIDATE_300DPI", "`"),
  if (release_mode) {
    "- Outputs: independently rendered 300-dpi preview, true 600-ppi PNG, LZW TIFF and physical-size Cairo PDF."
  } else {
    "- Outputs: selected 300-dpi RGB QC PNG candidate(s) only"
  }
)
writeLines(provenance, file.path(manifest_root, "RENDERER_PROVENANCE_SUPP_SYNTHESIS.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(manifest_root, "R_SESSION_INFO_SUPP_SYNTHESIS.txt"), useBytes = TRUE)
cat("PASS: selected synthesis figure(s) rendered with zero captured warnings; route=",
    paste(figure_set, collapse = ","), "; mode=", if (release_mode) "RELEASE" else "CANDIDATE", "\n", sep = "")
