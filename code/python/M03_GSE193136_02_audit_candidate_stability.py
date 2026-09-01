from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M03_WITHIN_COHORT_EFFECTS/runs/20260814_M03_B3_GSE193136"
OUT = RUN / "results"

effects = pd.read_csv(OUT / "GSE193136_DESEQ2_AGE_ADJUSTED_DOR_MINUS_NOR_ALL_GENES.csv")
unadjusted = pd.read_csv(OUT / "GSE193136_DESEQ2_UNADJUSTED_DOR_MINUS_NOR_ALL_GENES.csv")
loso = pd.read_csv(OUT / "GSE193136_LEAVE_ONE_SAMPLE_OUT_GENE_STABILITY.csv")

audited = effects.merge(
    loso[["gene_id_version", "loso_same_direction_fraction", "loso_min_log2FoldChange", "loso_max_log2FoldChange"]],
    on="gene_id_version",
).merge(
    unadjusted[["gene_id_version", "log2FoldChange"]].rename(columns={"log2FoldChange": "unadjusted_log2FoldChange"}),
    on="gene_id_version", how="left",
)
audited["unadjusted_same_direction"] = audited["log2FoldChange"] * audited["unadjusted_log2FoldChange"] > 0
candidates = audited.loc[(audited["padj"] < 0.05) & (audited["log2FoldChange"].abs() >= 1)].copy()
candidates["loso_direction_class"] = pd.cut(
    candidates["loso_same_direction_fraction"], [-0.01, 0.749, 0.916, 1.001],
    labels=["UNSTABLE_<9/12", "INTERMEDIATE_9-11/12", "STABLE_12/12"],
)
candidates = candidates.sort_values(["loso_same_direction_fraction", "padj"], ascending=[False, True])
candidates.to_csv(OUT / "GSE193136_FDR_CANDIDATE_STABILITY.csv", index=False)

stable_11 = int((candidates["loso_same_direction_fraction"] >= 11 / 12).sum())
all_12 = int((candidates["loso_same_direction_fraction"] == 1).sum())
same_unadjusted = int(candidates["unadjusted_same_direction"].sum())
summary = f"""# GSE193136 FDR Candidate Stability Audit

- Age-adjusted FDR < 0.05 and |log2FC| >= 1 candidates: {len(candidates)}
- Retain direction in at least 11/12 leave-one-sample-out fits: {stable_11}/{len(candidates)}
- Retain direction in all 12 leave-one-sample-out fits: {all_12}/{len(candidates)}
- Same direction in the unadjusted model: {same_unadjusted}/{len(candidates)}

These are GSE193136 cohort-level candidates. Cross-cohort replication remains required.
"""
(OUT / "GSE193136_FDR_CANDIDATE_STABILITY_SUMMARY.md").write_text(summary, encoding="utf-8")
(RUN / "STEP09_M03_CANDIDATE_AUDIT.PASS.txt").write_text("PASS\tFDR candidate stability audited\n", encoding="utf-8")
print(summary)

