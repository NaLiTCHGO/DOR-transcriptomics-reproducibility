#!/usr/bin/env python3
"""Independent audit for M06 pathway LOCO outputs and no-leakage rules."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


COHORTS = ["GSE274832", "GSE193136", "GSE232306"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    run_dir = parse_args().run_dir.resolve()
    results = run_dir / "results"
    plots = results / "plots"
    checks: list[dict] = []
    required = [
        "M06_LOCO_PATHWAY_EVALUATION.csv",
        "M06_COLLECTION_LOCO_SUMMARY.csv",
        "M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv",
        "M06_ALL_PATHWAY_ROTATION_MATRIX.csv",
        "M06_FROZEN_M05_STRICT_SURVIVAL.csv",
        "M06_ROBUSTNESS_TIER_COUNTS.csv",
        "M06_UNIVERSAL_CORE_LEADING_EDGE.csv",
        "M06_KEY_COUNTS.csv",
        "M06_LEAVE_ONE_COHORT_OUT_SUMMARY.md",
    ]
    for name in required:
        path = results / name
        checks.append({"check": f"required_output::{name}", "pass": path.is_file() and path.stat().st_size > 0, "detail": str(path)})
    plot_names = [
        "M06_LOCO_REPLICATION_FUNNEL.png",
        "M06_STRICT_REPLICATION_RATES.png",
        "M06_UNIVERSAL_CORE_NES_HEATMAP.png",
        "M06_GLOBAL_PAIR_CONTEXT.png",
        "M06_HALLMARK_FROZEN_STRICT_SURVIVAL.png",
        "M06_ROBUSTNESS_TIER_COUNTS.png",
    ]
    for name in plot_names:
        path = plots / name
        checks.append({"check": f"required_plot::{name}", "pass": path.is_file() and path.stat().st_size >= 10000, "detail": f"{path}; bytes={path.stat().st_size if path.is_file() else 0}"})

    try:
        evaluation = pd.read_csv(results / "M06_LOCO_PATHWAY_EVALUATION.csv")
        summary = pd.read_csv(results / "M06_COLLECTION_LOCO_SUMMARY.csv")
        core = pd.read_csv(results / "M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv")
        all_matrix = pd.read_csv(results / "M06_ALL_PATHWAY_ROTATION_MATRIX.csv")
        frozen = pd.read_csv(results / "M06_FROZEN_M05_STRICT_SURVIVAL.csv")
        tier = pd.read_csv(results / "M06_ROBUSTNESS_TIER_COUNTS.csv")
        leading = pd.read_csv(results / "M06_UNIVERSAL_CORE_LEADING_EDGE.csv")
        text = (results / "M06_LEAVE_ONE_COHORT_OUT_SUMMARY.md").read_text(encoding="utf-8")

        selected_recomputed = (
            evaluation["retained_same_direction"].astype(bool)
            & (evaluation["retained_a_padj"] < 0.05)
            & (evaluation["retained_b_padj"] < 0.05)
            & (evaluation["retained_median_abs_nes"] >= 1.0)
        )
        direction_recomputed = selected_recomputed & (
            np.sign(evaluation["held_out_nes"])
            == np.sign(evaluation["retained_a_nes"] + evaluation["retained_b_nes"])
        )
        strict_recomputed = direction_recomputed & (evaluation["held_out_padj"] < 0.05) & (evaluation["held_out_nes"].abs() >= 1.0)
        universal_recomputed = all_matrix[
            [f"pair_selected_without_{c}" for c in COHORTS] + [f"strict_replicated_in_{c}" for c in COHORTS]
        ].astype(bool).all(axis=1)
        checks.extend(
            [
                {"check": "three_rotations_two_collections", "pass": len(summary) == 6 and set(summary["held_out_cohort"]) == set(COHORTS) and set(summary["collection"]) == {"HALLMARK", "REACTOME"}, "detail": summary[["held_out_cohort", "collection"]].to_dict("records")},
                {"check": "all_eligible_pathways_each_rotation", "pass": len(evaluation) == 3 * (49 + 1016) and evaluation.groupby(["held_out_cohort", "collection"])["pathway"].nunique().sum() == len(evaluation), "detail": f"rows={len(evaluation)}"},
                {"check": "pair_selection_rule_recomputed", "pass": np.array_equal(evaluation["pair_selected"].astype(bool).to_numpy(), selected_recomputed.to_numpy()), "detail": f"selected={int(selected_recomputed.sum())}"},
                {"check": "heldout_direction_rule_recomputed", "pass": np.array_equal(evaluation["held_out_direction_replication"].astype(bool).to_numpy(), direction_recomputed.to_numpy()), "detail": f"direction={int(direction_recomputed.sum())}"},
                {"check": "heldout_strict_rule_recomputed", "pass": np.array_equal(evaluation["held_out_strict_replication"].astype(bool).to_numpy(), strict_recomputed.to_numpy()), "detail": f"strict={int(strict_recomputed.sum())}"},
                {"check": "universal_core_recomputed", "pass": np.array_equal(all_matrix["universal_loco_core"].astype(bool).to_numpy(), universal_recomputed.to_numpy()) and len(core) == int(universal_recomputed.sum()), "detail": f"core={len(core)}"},
                {"check": "frozen_strict_survival_complete", "pass": len(frozen) == 55 and frozen["surviving_rotation_n"].between(0, 3).all(), "detail": f"rows={len(frozen)}; tier_rows={len(tier)}"},
                {"check": "core_leading_edges_match", "pass": len(leading) == len(core) and not leading.duplicated(["collection", "pathway"]).any(), "detail": f"leading={len(leading)}; core={len(core)}"},
                {"check": "summary_internal_validation_boundary", "pass": "not external clinical validation" in text and "not pathway activation/inhibition" in text and "No fgsea result or threshold was refit" in text, "detail": "required no-leakage and claim text"},
            ]
        )
    except Exception as exc:
        checks.append({"check": "table_level_audit", "pass": False, "detail": repr(exc)})

    table = pd.DataFrame(checks)
    table.to_csv(results / "M06_AUDIT_CHECKS.csv", index=False)
    failed = table.loc[~table["pass"].astype(bool)]
    status = "PASS_WITH_LIMITATION" if failed.empty else "FAIL"
    lines = [
        "# M06 Independent Output Audit",
        "",
        f"`STATUS: {status}`",
        "",
        f"- Checks passed: {int(table['pass'].sum())}/{len(table)}.",
        f"- Checks failed: {len(failed)}.",
        "- PASS_WITH_LIMITATION remains mandatory because this is conditional internal robustness across the same three cohorts, not external validation, and pathway redundancy/cohort limitations remain.",
        "",
        "## Failed checks",
        "",
    ]
    if failed.empty:
        lines.append("None.")
    else:
        lines.extend(f"- `{row.check}`: {row.detail}" for row in failed.itertuples(index=False))
    lines.extend(["", "## Check table", "", "See `M06_AUDIT_CHECKS.csv`.", ""])
    (results / "M06_AUDIT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    marker = run_dir / ("M06_PASS_WITH_LIMITATION.txt" if failed.empty else "M06_FAIL.txt")
    marker.write_text(json.dumps({"module": "M06_LEAVE_ONE_COHORT_OUT", "status": status, "checks_passed": int(table["pass"].sum()), "checks_total": int(len(table)), "failed_checks": failed["check"].tolist()}, indent=2), encoding="utf-8")
    print(f"M06 audit status: {status}; checks={int(table['pass'].sum())}/{len(table)}")
    return 0 if failed.empty else 1


if __name__ == "__main__":
    raise SystemExit(main())
