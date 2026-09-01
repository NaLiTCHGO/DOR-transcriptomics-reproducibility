from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B4_GSE232306"
CONFIG = ROOT / "04_code/configs/GSE232306_SAMPLES.csv"
FASTP = RUN / "results/fastp"
OUT = RUN / "results/raw_qc"
OUT.mkdir(parents=True, exist_ok=True)

cfg = pd.read_csv(CONFIG, dtype=str)
rows = []
for row in cfg.itertuples(index=False):
    path = FASTP / f"{row.srr}.fastp.json"
    if not path.exists():
        raise FileNotFoundError(path)
    obj = json.loads(path.read_text(encoding="utf-8"))
    before = obj["summary"]["before_filtering"]
    q30 = float(before["q30_rate"])
    n_rate = float(before.get("n_content", 0.0)) / max(float(before["total_bases"]), 1.0)
    verdict = "PASS" if q30 >= 0.80 and n_rate <= 0.05 else "REVIEW"
    rows.append(
        {
            "srr": row.srr,
            "gsm": row.gsm,
            "phenotype": row.phenotype,
            "total_reads": int(before["total_reads"]),
            "total_bases": int(before["total_bases"]),
            "q20_rate": float(before["q20_rate"]),
            "q30_rate": q30,
            "gc_content": float(before["gc_content"]),
            "read1_mean_length": float(before["read1_mean_length"]),
            "read2_mean_length": float(before["read2_mean_length"]),
            "duplication_rate": float(obj["duplication"]["rate"]),
            "n_base_fraction": n_rate,
            "raw_qc_verdict": verdict,
        }
    )

summary = pd.DataFrame(rows).sort_values("srr")
if len(summary) != 12:
    raise RuntimeError(f"Expected 12 samples; found {len(summary)}")
summary.to_csv(OUT / "GSE232306_FASTP_RAW_QC_SUMMARY.csv", index=False)

group = summary.groupby("phenotype").agg(
    n=("srr", "size"),
    reads_median=("total_reads", "median"),
    q30_min=("q30_rate", "min"),
    q30_max=("q30_rate", "max"),
    gc_min=("gc_content", "min"),
    gc_max=("gc_content", "max"),
).reset_index()
group.to_csv(OUT / "GSE232306_FASTP_GROUP_SUMMARY.csv", index=False)

status = "PASS" if (summary["raw_qc_verdict"] == "PASS").all() else "PASS_WITH_REVIEW"
(RUN / "STEP03_RAW_QC.PASS.txt").write_text(
    f"{status}\t12/12 raw QC summarized\n", encoding="utf-8"
)
print(f"{status}: GSE232306 raw QC summarized for 12/12 samples; no per-sample age values are public")
print(summary[["srr", "phenotype", "total_reads", "q30_rate", "gc_content", "raw_qc_verdict"]].to_string(index=False))
