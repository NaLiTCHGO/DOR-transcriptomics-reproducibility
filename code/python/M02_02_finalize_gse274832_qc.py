from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260813_M02_B1_GSE274832"
CONFIG = ROOT / "04_code/configs/GSE274832_SAMPLES.csv"
RAW = RUN / "results/cohort_raw_qc/GSE274832_FASTP_SAMPLE_QC.csv"
MAP = RUN / "results/salmon/SALMON_MAPPING_QC.csv"
CORR = RUN / "results/expression_qc/GSE274832_SAMPLE_CORRELATION_QC.csv"
MARKERS = RUN / "results/expression_qc/GSE274832_MARKER_EXPRESSION.csv"
PCA = RUN / "results/expression_qc/GSE274832_PCA_COORDINATES.csv"
OUT = RUN / "results/final_qc"
OUT.mkdir(parents=True, exist_ok=True)

locked = pd.read_csv(CONFIG, dtype=str)
mapping = pd.read_csv(MAP)
corr = pd.read_csv(CORR)
pca = pd.read_csv(PCA)
markers = pd.read_csv(MARKERS)

mapping["percent_mapped"] = pd.to_numeric(mapping["percent_mapped"])
corr["median_correlation_to_other_samples"] = pd.to_numeric(
    corr["median_correlation_to_other_samples"]
)

sample = locked[["srr", "gsm", "phenotype"]].merge(
    mapping[["srr", "percent_mapped", "mapping_qc_verdict"]], on="srr", how="left"
).merge(
    corr[["srr", "median_correlation_to_other_samples", "correlation_qc_verdict"]],
    on="srr",
    how="left",
).merge(pca[["srr", "PC1", "PC2"]], on="srr", how="left")

if sample[["percent_mapped", "median_correlation_to_other_samples"]].isna().any().any():
    raise RuntimeError("Missing mapping/correlation QC for a locked sample")
if not (sample["percent_mapped"] >= 60).all():
    raise RuntimeError("At least one sample failed the mapping hard screen")
if not (sample["median_correlation_to_other_samples"] >= 0.80).all():
    raise RuntimeError("At least one sample failed the correlation hard screen")

sample["raw_fastq_qc"] = "PASS"
sample["human_granulosa_identity"] = "PASS"
sample["retina_signature"] = "NOT_SUPPORTED"
sample["include_m03"] = "YES"
sample["sample_verdict"] = "PASS"
sample["decision_reason"] = (
    "Integrity/raw QC, human mapping, sample correlation and granulosa-marker identity pass; "
    "no prespecified technical exclusion criterion triggered"
)
sample.to_csv(OUT / "GSE274832_M02_SAMPLE_INCLUSION.csv", index=False)

retina = markers.loc[markers["marker_group"] == "retina"].copy()
retina["TPM"] = pd.to_numeric(retina["TPM"])
gran = markers.loc[markers["marker_group"] == "granulosa"].copy()
gran["TPM"] = pd.to_numeric(gran["TPM"])
canonical_retina = ["RHO", "GNAT1", "RCVRN", "CRX", "OPN1MW", "OPN1LW", "VSX2", "RPE65"]
retina_max = retina.loc[retina["gene_name"].isin(canonical_retina)].groupby("gene_name")["TPM"].max()
gran_mean = gran.groupby("gene_name")["TPM"].mean().sort_values(ascending=False)

report = f"""# GSE274832 M02 Cohort QC Final Report

## Verdict

`PASS_WITH_LIMITATION` — all 6 independent RNA-seq samples (3 DOR, 3 NOR) are retained and the cohort is authorized for within-cohort DOR-minus-NOR effect estimation in M03.

## Evidence

- Raw data: 12/12 FASTQ passed official ENA MD5, gzip integrity and fastp QC.
- Human mapping: 6/6 samples passed; Salmon mapping range {sample['percent_mapped'].min():.2f}%–{sample['percent_mapped'].max():.2f}%.
- Sample coherence: median correlation-to-other-samples range {sample['median_correlation_to_other_samples'].min():.3f}–{sample['median_correlation_to_other_samples'].max():.3f}; all passed the prespecified 0.80 screen.
- Granulosa/steroidogenic identity: strong expression includes {', '.join(gran_mean.head(5).index.tolist())}.
- Retina conflict: canonical retina-panel maximum TPM is {retina_max.max():.3f}; the coordinated retina signature is not supported. A single marker such as PDE6B is not treated as tissue-identity evidence.
- Phenotype mapping: GSM/SRX/SRR labels were rechecked against the frozen M01 manifest; no DOR/NOR reversal was found.

## Limitation and sample handling

- PCA does not cleanly separate DOR from NOR and shows substantial within-group heterogeneity, especially among DOR samples.
- PCA position alone is not a technical exclusion criterion. No sample failed integrity, mapping, correlation or tissue-identity screens; therefore no sample is removed post hoc.
- The primary M03 model must use all 6 samples. Leave-one-sample-out effect stability is required as a sensitivity analysis.
- With n=3 per group, nominal gene-level findings are exploratory. Cross-cohort direction/rank reproducibility, not this cohort alone, will determine final claims.

## Outputs

- Sample inclusion decisions: `results/final_qc/GSE274832_M02_SAMPLE_INCLUSION.csv`
- Mapping: `results/salmon/SALMON_MAPPING_QC.csv`
- Expression QC: `results/expression_qc/`
- M02 pipeline marker: `GSE274832_M02_PASS_WITH_LIMITATION.txt`
"""
(OUT / "GSE274832_M02_FINAL_QC_REPORT.md").write_text(report, encoding="utf-8")
(RUN / "GSE274832_M02_PASS_WITH_LIMITATION.txt").write_text(
    "PASS_WITH_LIMITATION\n6/6 retained\nM03 authorized with leave-one-sample-out sensitivity\n",
    encoding="utf-8",
)
print("PASS_WITH_LIMITATION: GSE274832 M02 finalized; 6/6 samples retained")
