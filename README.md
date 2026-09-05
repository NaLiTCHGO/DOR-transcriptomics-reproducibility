# DOR transcriptomics minimal reproducible code

Release `v1.0.1` contains the canonical analysis and accepted-figure code required to reproduce the reported DOR cross-cohort workflow. It intentionally excludes raw sequencing data, generated results, manuscript files, author contact information, workstation metadata, historical snapshots, failed experiments, and local tool caches.

## Scope

- cohort-specific FASTQ/QC/Salmon/DESeq2 workflows for GSE274832, GSE193136 and GSE232306;
- cross-cohort gene-effect synthesis and heterogeneity;
- Hallmark/Reactome convergence and no-leakage pathway LOCO;
- final reporting and accepted Figure 1-Figure 5 / Figure S1-Figure S6 generation and QC code;
- frozen sample/configuration, software-version, parameter and random-seed manifests.

## Use

Set the project root explicitly when running PowerShell entry points, or set `DOR_PROJECT_ROOT` for R scripts. External executables and the two figure skills are not bundled; provide their paths in the relevant command-line parameters. Public raw data are identified in `PUBLIC_ACCESSIONS.csv`.

The package starts from public accessions and does not redistribute FASTQ, reference genomes, Salmon indices, large matrices, or journal manuscript content.

## Accepted figure entry points

- Figure 1: `code/python/FIGURE1_04_generate_final_submission.py` (compatibility entry point) or `figure1_reproducibility/code/generate_Figure1.py --package-root figure1_reproducibility`.
- Figure 3-Figure 5: `code/R/PUB_FIGBATCH_02_render_main.R`.
- Figure S4-Figure S6: `code/R/PUB_FIGBATCH_04_render_supp_synthesis.R`.

The files `code/final_refresh/Figure_1_and_5_refresh.py` and `code/final_refresh/OPTIONAL_02_render_figure4_wording.R` are retained only as explicit retirement guards. They are not accepted rendering entry points.

## Version note

`v1.0.1` is a reproducibility-consistency patch. It changes figure-generation entry points, labels, documentation and citation metadata; it does not change the scientific analyses or reported numerical results. See `CHANGELOG.md`.
