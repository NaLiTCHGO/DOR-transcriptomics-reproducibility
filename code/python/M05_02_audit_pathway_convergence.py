#!/usr/bin/env python3
"""Independent structural and interpretation-boundary audit for M05."""

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
    run_dir = parse_args().run_dir.resolve()
    results = run_dir / "results"
    plots = results / "plots"
    checks: list[dict] = []

    required = [
        "M05_RANK_INPUT_SYMBOLS.csv",
        "M05_RANK_INPUT_AUDIT.csv",
        "M05_GENESET_AUDIT.csv",
        "M05_FGSEA_ALL_RESULTS.csv",
        "M05_FGSEA_NUMERICAL_WARNING_AUDIT.csv",
        "M05_PATHWAY_CONVERGENCE.csv",
        "M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv",
        "M05_STRICT_CONVERGENT_PATHWAYS.csv",
        "M05_STRICT_PATHWAY_LEADING_EDGE_OVERLAP.csv",
        "M05_KEY_COUNTS.csv",
        "M05_PATHWAY_CONVERGENCE_SUMMARY.md",
    ]
    for name in required:
        path = results / name
        checks.append({"check": f"required_output::{name}", "pass": path.is_file() and path.stat().st_size > 0, "detail": str(path)})

    required_plots = [
        "M05_HALLMARK_NES_HEATMAP.png",
        "M05_REACTOME_TOP30_NES_HEATMAP.png",
        "M05_PATHWAY_VS_GENE_REPRODUCIBILITY.png",
        "M05_PAIRWISE_NES_SCATTER.png",
        "M05_CONVERGENCE_COUNTS.png",
        "M05_PATHWAY_NES_DISPERSION.png",
    ]
    for name in required_plots:
        path = plots / name
        checks.append({"check": f"required_plot::{name}", "pass": path.is_file() and path.stat().st_size >= 10000, "detail": f"{path}; bytes={path.stat().st_size if path.is_file() else 0}"})

    try:
        rank_symbols = pd.read_csv(results / "M05_RANK_INPUT_SYMBOLS.csv")
        rank_audit = pd.read_csv(results / "M05_RANK_INPUT_AUDIT.csv")
        gene_sets = pd.read_csv(results / "M05_GENESET_AUDIT.csv")
        fgsea = pd.read_csv(results / "M05_FGSEA_ALL_RESULTS.csv")
        numerical = pd.read_csv(results / "M05_FGSEA_NUMERICAL_WARNING_AUDIT.csv")
        convergence = pd.read_csv(results / "M05_PATHWAY_CONVERGENCE.csv")
        pairwise = pd.read_csv(results / "M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv")
        strict = pd.read_csv(results / "M05_STRICT_CONVERGENT_PATHWAYS.csv")
        leading = pd.read_csv(results / "M05_STRICT_PATHWAY_LEADING_EDGE_OVERLAP.csv")
        summary = (results / "M05_PATHWAY_CONVERGENCE_SUMMARY.md").read_text(encoding="utf-8")

        eligible = dict(zip(gene_sets["collection"], gene_sets["eligible_sets_15_to_500"]))
        expected_fgsea = 3 * int(sum(eligible.values()))
        expected_convergence = int(sum(eligible.values()))
        checks.extend(
            [
                {"check": "three_identical_full_rank_vectors", "pass": len(rank_audit) == 3 and rank_audit["unique_rank_symbols"].nunique() == 1 and int(rank_audit["unique_rank_symbols"].iloc[0]) >= 10000 and set(rank_audit["prefilter"]) == {"NONE_FULL_COMMON_MAPPED_UNIVERSE"}, "detail": rank_audit[["cohort", "unique_rank_symbols", "prefilter"]].to_dict("records")},
                {"check": "unique_gene_symbols", "pass": not rank_symbols["gene_symbol"].duplicated().any() and len(rank_symbols) == int(rank_audit["unique_rank_symbols"].iloc[0]), "detail": f"rows={len(rank_symbols)}; duplicates={int(rank_symbols['gene_symbol'].duplicated().sum())}"},
                {"check": "frozen_collection_set_counts", "pass": int(eligible.get("HALLMARK", 0)) >= 45 and int(eligible.get("REACTOME", 0)) >= 900, "detail": eligible},
                {"check": "all_fgsea_results_present", "pass": len(fgsea) == expected_fgsea and fgsea.groupby(["collection", "cohort"])["pathway"].nunique().eq(fgsea.groupby(["collection", "cohort"]).size()).all(), "detail": f"rows={len(fgsea)}; expected={expected_fgsea}"},
                {"check": "fgsea_numerical_warning_audited", "pass": len(numerical) == 6 and (numerical["pvalue_at_eps_1e_10_n"] == numerical["log2err_missing_n"]).all() and "upper-bounded/capped" in " ".join(numerical["interpretation"].astype(str)), "detail": f"eps_floor_total={int(numerical['pvalue_at_eps_1e_10_n'].sum())}"},
                {"check": "pathway_convergence_rows_complete", "pass": len(convergence) == expected_convergence and not convergence.duplicated(["collection", "pathway"]).any(), "detail": f"rows={len(convergence)}; expected={expected_convergence}"},
                {"check": "six_pairwise_collection_results", "pass": len(pairwise) == 6, "detail": f"rows={len(pairwise)}"},
                {"check": "finite_pathway_reproducibility", "pass": np.isfinite(pairwise[["pathway_nes_spearman_rho", "pathway_nes_direction_concordance"]].to_numpy()).all(), "detail": "pathway rho and direction"},
                {"check": "strict_rule_recomputed", "pass": np.array_equal(convergence["strict_pathway_convergence"].astype(bool).to_numpy(), ((convergence["all_three_same_nes_direction"].astype(bool)) & (convergence["fgsea_fdr_lt_0_05_cohort_n"] >= 2) & (convergence["median_abs_nes"] >= 1.0)).to_numpy()), "detail": f"strict_rows={len(strict)}"},
                {"check": "strict_leading_edge_rows_match", "pass": len(leading) == len(strict), "detail": f"leading={len(leading)}; strict={len(strict)}"},
                {"check": "summary_full_rank_and_claim_boundaries", "pass": "No candidate, p-value, meta-FDR, FMNL1, or PGAP1 prefilter" in summary and "does not prove pathway activation or inhibition" in summary and "Signed Stouffer FDR is supportive" in summary and "no exact smaller p-value is claimed" in summary, "detail": "required interpretation text"},
            ]
        )
    except Exception as exc:
        checks.append({"check": "table_level_audit", "pass": False, "detail": repr(exc)})

    table = pd.DataFrame(checks)
    table.to_csv(results / "M05_AUDIT_CHECKS.csv", index=False)
    failed = table.loc[~table["pass"].astype(bool)]
    status = "PASS_WITH_LIMITATION" if failed.empty else "FAIL"
    lines = [
        "# M05 Independent Output Audit",
        "",
        f"`STATUS: {status}`",
        "",
        f"- Checks passed: {int(table['pass'].sum())}/{len(table)}.",
        f"- Checks failed: {len(failed)}.",
        "- A technically complete result remains PASS_WITH_LIMITATION because M05 inherits three cohort-specific limitations, Reactome redundancy, symbol mapping loss, and unresolved GSE232306 dependence pending M06.",
        "- Zero strict pathways would be an allowed outcome and would not fail this audit.",
        "",
        "## Failed checks",
        "",
    ]
    if failed.empty:
        lines.append("None.")
    else:
        lines.extend(f"- `{row.check}`: {row.detail}" for row in failed.itertuples(index=False))
    lines.extend(["", "## Check table", "", "See `M05_AUDIT_CHECKS.csv`.", ""])
    (results / "M05_AUDIT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    marker = run_dir / ("M05_PASS_WITH_LIMITATION.txt" if failed.empty else "M05_FAIL.txt")
    marker.write_text(json.dumps({"module": "M05_PATHWAY_CONVERGENCE", "status": status, "checks_passed": int(table["pass"].sum()), "checks_total": int(len(table)), "failed_checks": failed["check"].tolist()}, indent=2), encoding="utf-8")
    print(f"M05 audit status: {status}; checks={int(table['pass'].sum())}/{len(table)}")
    return 0 if failed.empty else 1


if __name__ == "__main__":
    raise SystemExit(main())
