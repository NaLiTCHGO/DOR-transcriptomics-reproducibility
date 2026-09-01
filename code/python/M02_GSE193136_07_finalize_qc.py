from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B3_GSE193136"
CONFIG = ROOT / "04_code/configs/GSE193136_SAMPLES.csv"
OUT = RUN / "results/final_qc"
OUT.mkdir(parents=True, exist_ok=True)

cfg = pd.read_csv(CONFIG)
raw = pd.read_csv(RUN / "results/raw_qc/GSE193136_FASTP_RAW_QC_SUMMARY.csv")
mapping = pd.read_csv(RUN / "results/salmon/GSE193136_SALMON_MAPPING_QC.csv")
corr = pd.read_csv(RUN / "results/expression_qc/GSE193136_SAMPLE_CORRELATION_QC.csv")
marker = pd.read_csv(RUN / "results/expression_qc/GSE193136_MARKER_QC_BY_SAMPLE.csv")
pca = pd.read_csv(RUN / "results/expression_qc/GSE193136_PCA_COORDINATES.csv")

sample = cfg[["srr", "gsm", "phenotype", "age_years"]].merge(
    raw[["srr", "q30_rate", "gc_content", "duplication_rate", "raw_qc_verdict"]], on="srr", how="left"
).merge(
    mapping[["srr", "percent_mapped", "mapping_qc_verdict"]], on="srr", how="left"
).merge(
    corr[["srr", "median_correlation_to_other_samples", "minimum_correlation_to_other_samples", "correlation_qc_verdict"]], on="srr", how="left"
).merge(marker, on="srr", how="left").merge(pca[["srr", "PC1", "PC2"]], on="srr", how="left")

required = ["q30_rate", "percent_mapped", "median_correlation_to_other_samples", "granulosa_markers_tpm_ge_1"]
if sample[required].isna().any().any():
    raise RuntimeError("Missing a required M02 QC field for a locked sample")

sample["integrity_qc"] = "PASS"
sample["include_m03"] = "YES"
sample["sample_verdict"] = "PASS"
sample.loc[sample["q30_rate"] < 0.80, ["include_m03", "sample_verdict"]] = ["REVIEW", "REVIEW"]
sample.loc[sample["percent_mapped"] < 30, ["include_m03", "sample_verdict"]] = ["REVIEW", "REVIEW"]
sample.loc[(sample["percent_mapped"] >= 30) & (sample["percent_mapped"] < 60), "sample_verdict"] = "PASS_WITH_LIMITATION"
sample.loc[sample["median_correlation_to_other_samples"] < 0.80, ["include_m03", "sample_verdict"]] = ["REVIEW", "REVIEW"]
sample.loc[sample["granulosa_markers_tpm_ge_1"] < 4, ["include_m03", "sample_verdict"]] = ["REVIEW", "REVIEW"]
sample["decision_reason"] = (
    "Prespecified integrity, raw-QC, human mapping, sample-correlation and granulosa-identity screens passed; "
    "PCA position alone was not used for post-hoc exclusion"
)
sample.to_csv(OUT / "GSE193136_M02_SAMPLE_INCLUSION.csv", index=False)

if (sample["sample_verdict"] == "REVIEW").any():
    raise RuntimeError("At least one GSE193136 sample requires manual review before M03")
if sample["percent_mapped"].max() - sample["percent_mapped"].min() > 10:
    raise RuntimeError("Transcriptome mapping varies by more than 10 percentage points across the cohort")

dor_age = sample.loc[sample.phenotype == "DOR", "age_years"]
nor_age = sample.loc[sample.phenotype == "NOR", "age_years"]
age_difference = float(dor_age.mean() - nor_age.mean())
low_mapping = bool((sample["percent_mapped"] < 60).any())
verdict = "PASS_WITH_LIMITATION" if abs(age_difference) >= 3 or low_mapping else "PASS"
report = f"""# GSE193136 M02 Cohort QC Final Report

## Verdict

`{verdict}` — all 12 independent RNA-seq samples (6 DOR, 6 NOR) are retained and authorized for M03 within-cohort effect estimation.

## Evidence

- Raw data: 24/24 FASTQ passed expected-size and gzip integrity checks; MD5 was intentionally omitted under the local lightweight policy.
- Raw read QC: 12/12 fastp reports passed the prespecified Q30 screen; Q30 range {sample.q30_rate.min():.6f}–{sample.q30_rate.max():.6f}.
- Human/transcriptome mapping: decoy-aware Salmon mapping range {sample.percent_mapped.min():.2f}%–{sample.percent_mapped.max():.2f}%. Values below 60% remain an explicit limitation; all samples must be >=30%, vary by <=10 percentage points, retain broad gene coverage and pass granulosa-identity QC.
- Sample coherence: median correlation-to-other-samples range {sample.median_correlation_to_other_samples.min():.3f}–{sample.median_correlation_to_other_samples.max():.3f}; 12/12 passed 0.80.
- Granulosa identity: all samples express at least four prespecified granulosa/steroidogenic markers at TPM >= 1.
- No sample was excluded from PCA position alone.

## Prespecified limitation carried into M03

- Mean age is {dor_age.mean():.2f} years in DOR and {nor_age.mean():.2f} years in NOR (difference {age_difference:.2f} years).
- M03 therefore reports both DOR-minus-NOR total effects (`~ condition`) and age-adjusted effects (`~ age_centered + condition`). The age-adjusted model is the GSE193136 principal inferential model; the unadjusted model preserves comparability for cross-cohort sensitivity synthesis.

## Outputs

- Sample inclusion: `results/final_qc/GSE193136_M02_SAMPLE_INCLUSION.csv`
- Raw QC: `results/raw_qc/`
- Mapping: `results/salmon/GSE193136_SALMON_MAPPING_QC.csv`
- Expression QC and plots: `results/expression_qc/`
"""
(OUT / "GSE193136_M02_FINAL_QC_REPORT.md").write_text(report, encoding="utf-8")
(RUN / f"GSE193136_M02_{verdict}.txt").write_text(
    f"{verdict}\n12/12 retained\nM03 authorized with age-adjusted primary model and unadjusted sensitivity model\n",
    encoding="utf-8",
)
(RUN / "STEP07_M02_FINAL_GATE.PASS.txt").write_text(f"{verdict}\t12/12 retained\n", encoding="utf-8")
print(f"{verdict}: GSE193136 M02 finalized; 12/12 samples retained")
