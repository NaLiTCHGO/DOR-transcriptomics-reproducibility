#!/usr/bin/env python3
"""Independent structural and scientific-boundary audit for M04 outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    results = run_dir / "results"
    plots = results / "plots"
    checks: list[dict] = []

    required = [
        "M04_INPUT_AUDIT.csv",
        "M04_PAIRWISE_REPRODUCIBILITY.csv",
        "M04_TOPK_OVERLAP.csv",
        "M04_GENE_LEVEL_META_ANALYSIS.csv",
        "M04_CANDIDATE_OVERLAP_SUMMARY.csv",
        "M04_LOCO_STABILITY.csv",
        "M04_GENE_LEVEL_LOCO_EFFECTS.csv",
        "M04_KEY_COUNTS.csv",
        "M04_REPRO_HETEROGENEITY_SUMMARY.md",
    ]
    for name in required:
        path = results / name
        checks.append(
            {
                "check": f"required_output::{name}",
                "pass": path.is_file() and path.stat().st_size > 0,
                "detail": str(path),
            }
        )

    required_plots = [
        "M04_PAIRWISE_EFFECT_SCATTER.png",
        "M04_REPRODUCIBILITY_MATRICES.png",
        "M04_CANDIDATE_OVERLAP.png",
        "M04_HETEROGENEITY_DISTRIBUTION.png",
        "M04_TOP_CROSS_COHORT_EFFECT_HEATMAP.png",
        "M04_LOCO_STABILITY.png",
    ]
    for name in required_plots:
        path = plots / name
        checks.append(
            {
                "check": f"required_plot::{name}",
                "pass": path.is_file() and path.stat().st_size >= 10000,
                "detail": f"{path}; bytes={path.stat().st_size if path.is_file() else 0}",
            }
        )

    try:
        audit = pd.read_csv(results / "M04_INPUT_AUDIT.csv")
        pairwise = pd.read_csv(results / "M04_PAIRWISE_REPRODUCIBILITY.csv")
        meta = pd.read_csv(results / "M04_GENE_LEVEL_META_ANALYSIS.csv")
        loco = pd.read_csv(results / "M04_LOCO_STABILITY.csv")
        gene_loco = pd.read_csv(results / "M04_GENE_LEVEL_LOCO_EFFECTS.csv")
        summary_text = (results / "M04_REPRO_HETEROGENEITY_SUMMARY.md").read_text(encoding="utf-8")

        checks.extend(
            [
                {"check": "exactly_three_input_cohorts", "pass": len(audit) == 3, "detail": f"rows={len(audit)}"},
                {
                    "check": "common_universe_at_least_10000",
                    "pass": len(meta) >= 10000 and audit["common_valid_genes"].nunique() == 1,
                    "detail": f"meta_rows={len(meta)}; audit_common={audit['common_valid_genes'].tolist()}",
                },
                {
                    "check": "unique_unversioned_ensembl_gene_ids",
                    "pass": (not meta["gene_id"].duplicated().any())
                    and meta["gene_id"].astype(str).str.match(r"^ENSG[0-9]+$").all(),
                    "detail": f"duplicates={int(meta['gene_id'].duplicated().sum())}",
                },
                {"check": "three_pairwise_results", "pass": len(pairwise) == 3, "detail": f"rows={len(pairwise)}"},
                {
                    "check": "finite_pairwise_metrics",
                    "pass": np.isfinite(pairwise[["spearman_rho_log2fc", "direction_concordance_all_nonzero"]].to_numpy()).all(),
                    "detail": "Spearman and direction metrics",
                },
                {
                    "check": "finite_meta_effects_and_heterogeneity",
                    "pass": np.isfinite(meta[["fixed_beta", "random_beta", "i2_percent", "tau2_dl"]].to_numpy()).all(),
                    "detail": "fixed/random beta, I2, tau2",
                },
                {"check": "three_loco_results", "pass": len(loco) == 3, "detail": f"rows={len(loco)}"},
                {
                    "check": "gene_level_loco_matches_meta",
                    "pass": len(gene_loco) == len(meta) and not gene_loco["gene_id"].duplicated().any(),
                    "detail": f"loco_rows={len(gene_loco)}; meta_rows={len(meta)}",
                },
                {
                    "check": "comparison_direction_frozen",
                    "pass": set(audit["comparison"]) == {"DOR-minus-NOR"},
                    "detail": ";".join(sorted(set(audit["comparison"]))),
                },
                {
                    "check": "no_megamatrix_or_biomarker_claim",
                    "pass": "No sample-level expression matrix" in summary_text
                    and "not clinical biomarkers" in summary_text
                    and "Prohibited claims" in summary_text,
                    "detail": "summary contains required boundaries",
                },
            ]
        )
    except Exception as exc:
        checks.append({"check": "table_level_audit", "pass": False, "detail": repr(exc)})

    table = pd.DataFrame(checks)
    table.to_csv(results / "M04_AUDIT_CHECKS.csv", index=False)
    failed = table.loc[~table["pass"].astype(bool)]
    status = "PASS_WITH_LIMITATION" if failed.empty else "FAIL"
    report_lines = [
        "# M04 Independent Output Audit",
        "",
        f"`STATUS: {status}`",
        "",
        f"- Checks passed: {int(table['pass'].sum())}/{len(table)}.",
        f"- Checks failed: {len(failed)}.",
        "- PASS_WITH_LIMITATION is required even when all structural checks pass because all three upstream cohorts carry prespecified limitations and k=3 heterogeneity estimates are imprecise.",
        "- This audit verifies structure, estimand boundaries, and required outputs; it does not convert exploratory genes into biomarkers.",
        "",
        "## Failed checks",
        "",
    ]
    if failed.empty:
        report_lines.append("None.")
    else:
        for row in failed.itertuples(index=False):
            report_lines.append(f"- `{row.check}`: {row.detail}")
    report_lines.extend(["", "## Check table", "", "See `M04_AUDIT_CHECKS.csv`.", ""])
    (results / "M04_AUDIT_REPORT.md").write_text("\n".join(report_lines), encoding="utf-8")
    marker = run_dir / ("M04_PASS_WITH_LIMITATION.txt" if failed.empty else "M04_FAIL.txt")
    marker.write_text(
        json.dumps(
            {
                "module": "M04_REPRO_HETEROGENEITY",
                "status": status,
                "checks_passed": int(table["pass"].sum()),
                "checks_total": int(len(table)),
                "failed_checks": failed["check"].tolist(),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"M04 audit status: {status}; checks={int(table['pass'].sum())}/{len(table)}")
    return 0 if failed.empty else 1


if __name__ == "__main__":
    raise SystemExit(main())
