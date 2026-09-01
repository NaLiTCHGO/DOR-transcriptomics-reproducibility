import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
OUT = ROOT / "05_analysis_steps/M03_WITHIN_COHORT_EFFECTS/runs/20260814_M03_B1_GSE274832/results"

effects = pd.read_csv(OUT / "GSE274832_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv")
loso = pd.read_csv(OUT / "GSE274832_LEAVE_ONE_SAMPLE_OUT_GENE_STABILITY.csv")
official = pd.read_csv(OUT / "GSE274832_OFFICIAL_VS_REPROCESSED_EFFECTS.csv")

audited = effects.merge(
    loso[
        [
            "gene_id_version",
            "loso_same_direction_fraction",
            "loso_min_log2FoldChange",
            "loso_max_log2FoldChange",
        ]
    ],
    on="gene_id_version",
).merge(official[["gene_id", "official_log2FC"]], on="gene_id", how="left")

candidates = audited.loc[
    (audited["padj"] < 0.05) & (audited["log2FoldChange"].abs() >= 1)
].copy()
candidates["official_same_direction"] = (
    candidates["log2FoldChange"] * candidates["official_log2FC"] > 0
)
candidates["loso_direction_class"] = pd.cut(
    candidates["loso_same_direction_fraction"],
    [-0.01, 0.499, 0.832, 1.001],
    labels=["UNSTABLE_<3/6", "INTERMEDIATE_3-4/6", "STABLE_5-6/6"],
)
candidates = candidates.sort_values(
    ["loso_same_direction_fraction", "padj"], ascending=[False, True]
)
candidates.to_csv(OUT / "GSE274832_FDR_CANDIDATE_STABILITY.csv", index=False)

counts = candidates["loso_direction_class"].value_counts(dropna=False)
stable = int((candidates["loso_same_direction_fraction"] >= 5 / 6).sum())
all_six = int((candidates["loso_same_direction_fraction"] == 1).sum())
official_same = int(candidates["official_same_direction"].sum())

summary = f"""# GSE274832 FDR Candidate Stability Audit

- FDR < 0.05 and |log2FC| >= 1 candidates: {len(candidates)}
- Retain direction in at least 5/6 leave-one-sample-out fits: {stable}/{len(candidates)}
- Retain direction in all 6 leave-one-sample-out fits: {all_six}/{len(candidates)}
- Same direction in the submitter-provided FPKM table: {official_same}/{len(candidates)}
- Stability classes: {counts.to_dict()}

These are cohort-level candidates, not final biomarkers. Cross-cohort replication remains required.
"""
(OUT / "GSE274832_FDR_CANDIDATE_STABILITY_SUMMARY.md").write_text(
    summary, encoding="utf-8"
)

print(summary)
