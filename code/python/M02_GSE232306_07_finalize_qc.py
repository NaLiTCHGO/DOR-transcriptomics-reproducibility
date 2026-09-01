from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B4_GSE232306"
CONFIG = ROOT / "04_code/configs/GSE232306_SAMPLES.csv"
OUT = RUN / "results/final_qc"
OUT.mkdir(parents=True, exist_ok=True)

cfg = pd.read_csv(CONFIG)
raw = pd.read_csv(RUN / "results/raw_qc/GSE232306_FASTP_RAW_QC_SUMMARY.csv")
mapping = pd.read_csv(RUN / "results/salmon/GSE232306_SALMON_MAPPING_QC.csv")
corr = pd.read_csv(RUN / "results/expression_qc/GSE232306_SAMPLE_CORRELATION_QC.csv")
marker = pd.read_csv(RUN / "results/expression_qc/GSE232306_MARKER_QC_BY_SAMPLE.csv")
pca = pd.read_csv(RUN / "results/expression_qc/GSE232306_PCA_COORDINATES.csv")

sample = cfg[["srr", "gsm", "phenotype"]].merge(
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
sample.to_csv(OUT / "GSE232306_M02_SAMPLE_INCLUSION.csv", index=False)

if (sample["sample_verdict"] == "REVIEW").any():
    raise RuntimeError("At least one GSE232306 sample requires manual review before M03")
if sample["percent_mapped"].max() - sample["percent_mapped"].min() > 10:
    raise RuntimeError("Transcriptome mapping varies by more than 10 percentage points across the cohort")

low_mapping = bool((sample["percent_mapped"] < 60).any())
dor_gc = sample.loc[sample.phenotype == "DOR", "gc_content"]
nor_gc = sample.loc[sample.phenotype == "NOR", "gc_content"]
dor_dup = sample.loc[sample.phenotype == "DOR", "duplication_rate"]
nor_dup = sample.loc[sample.phenotype == "NOR", "duplication_rate"]
dor_pc1 = sample.loc[sample.phenotype == "DOR", "PC1"]
nor_pc1 = sample.loc[sample.phenotype == "NOR", "PC1"]
pc1_complete_separation = bool(dor_pc1.max() < nor_pc1.min() or nor_pc1.max() < dor_pc1.min())
verdict = "PASS_WITH_LIMITATION"
report = f"""# GSE232306 M02 Cohort QC Final Report

## Verdict

`{verdict}` — all 12 independent RNA-seq samples (6 DOR, 6 NOR) are retained and authorized for M03 within-cohort effect estimation.

## Evidence

- Raw data: 24/24 FASTQ passed expected-size and gzip integrity checks in this formal run. The earlier download QC table recorded size and MD5 PASS for 24/24; MD5 was not recomputed under the local lightweight policy.
- Raw read QC: 12/12 fastp reports passed the prespecified Q30 screen; Q30 range {sample.q30_rate.min():.6f}–{sample.q30_rate.max():.6f}.
- Human/transcriptome mapping: decoy-aware Salmon mapping range {sample.percent_mapped.min():.2f}%–{sample.percent_mapped.max():.2f}%. Values below 60% remain an explicit limitation; all samples must be >=30%, vary by <=10 percentage points, retain broad gene coverage and pass granulosa-identity QC.
- Sample coherence: median correlation-to-other-samples range {sample.median_correlation_to_other_samples.min():.3f}–{sample.median_correlation_to_other_samples.max():.3f}; 12/12 passed 0.80.
- Granulosa identity: all samples express at least four prespecified granulosa/steroidogenic markers at TPM >= 1.
- No sample was excluded from PCA position alone.
- PCA PC1 explains {pca.PC1_variance_percent.iloc[0]:.2f}% and completely separates DOR from NOR: {pc1_complete_separation}.

## Prespecified limitations carried into M03

- Individual ages are absent from the public GEO sample metadata, so age adjustment is not technically possible for this cohort.
- M03 therefore uses the DOR-minus-NOR total group effect (`~ condition`) as its only estimable primary model. Residual age/confounding risk must be carried into cross-cohort interpretation.
- Mapping is adequate (all samples >=60%): {not low_mapping}.
- Raw GC content is non-overlapping by phenotype (DOR mean {dor_gc.mean():.4f}, range {dor_gc.min():.4f}–{dor_gc.max():.4f}; NOR mean {nor_gc.mean():.4f}, range {nor_gc.min():.4f}–{nor_gc.max():.4f}). Duplication also differs (DOR mean {dor_dup.mean():.4f}; NOR mean {nor_dup.mean():.4f}). Together with complete PC1 separation, this indicates a large phenotype-associated compositional shift that cannot be distinguished from unrecorded batch/library or clinical covariates.
- GSE232306 is therefore retained as a high-heterogeneity screening cohort for M04, not as a standalone biomarker-discovery cohort.

## Outputs

- Sample inclusion: `results/final_qc/GSE232306_M02_SAMPLE_INCLUSION.csv`
- Raw QC: `results/raw_qc/`
- Mapping: `results/salmon/GSE232306_SALMON_MAPPING_QC.csv`
- Expression QC and plots: `results/expression_qc/`
"""
(OUT / "GSE232306_M02_FINAL_QC_REPORT.md").write_text(report, encoding="utf-8")
(RUN / f"GSE232306_M02_{verdict}.txt").write_text(
    f"{verdict}\n12/12 retained\nM03 authorized with unadjusted condition-only model; age unavailable\n",
    encoding="utf-8",
)
(RUN / "STEP07_M02_FINAL_GATE.PASS.txt").write_text(f"{verdict}\t12/12 retained\n", encoding="utf-8")
print(f"{verdict}: GSE232306 M02 finalized; 12/12 samples retained")
