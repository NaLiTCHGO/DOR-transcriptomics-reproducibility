# DOR transcriptomics minimal reproducible code

This repository-ready package contains only the canonical analysis and accepted-figure code required to reproduce the reported DOR cross-cohort workflow. It intentionally excludes raw sequencing data, generated results, manuscript files, author contact information, workstation metadata, historical snapshots, failed experiments, and local tool caches.

## Scope

- cohort-specific FASTQ/QC/Salmon/DESeq2 workflows for GSE274832, GSE193136 and GSE232306;
- cross-cohort gene-effect synthesis and heterogeneity;
- Hallmark/Reactome convergence and no-leakage pathway LOCO;
- final reporting and accepted Figure 1-Figure 5 / Figure S1-Figure S6 generation and QC code, including the final language-polished refresh scripts;
- frozen sample/configuration, software-version, parameter and random-seed manifests.

## Use

Set the project root explicitly when running PowerShell entry points, or set `DOR_PROJECT_ROOT` for R scripts. External executables and the two figure skills are not bundled; provide their paths in the relevant command-line parameters. Public raw data are identified in `PUBLIC_ACCESSIONS.csv`.

The package starts from public accessions and does not redistribute FASTQ, reference genomes, Salmon indices, large matrices, or journal manuscript content.
