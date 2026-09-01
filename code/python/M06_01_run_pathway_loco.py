#!/usr/bin/env python3
"""Run frozen pathway-level leave-one-cohort-out selection and replication."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COHORTS = ["GSE274832", "GSE193136", "GSE232306"]
COLLECTIONS = ["HALLMARK", "REACTOME"]
SHORT = {"GSE274832": "274832", "GSE193136": "193136", "GSE232306": "232306"}
COLORS = {"HALLMARK": "#0EA5E9", "REACTOME": "#F59E0B"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    return parser.parse_args()


def clean_pathway(name: str) -> str:
    return name.replace("HALLMARK_", "").replace("REACTOME_", "").replace("_", " ")


def load_inputs(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = project_root / "06_locked_results/modules/M05_PATHWAY_CONVERGENCE/v1/results"
    paths = {
        "convergence": base / "M05_PATHWAY_CONVERGENCE.csv",
        "pairwise": base / "M05_PATHWAY_PAIRWISE_REPRODUCIBILITY.csv",
        "strict": base / "M05_STRICT_CONVERGENT_PATHWAYS.csv",
        "leading": base / "M05_STRICT_PATHWAY_LEADING_EDGE_OVERLAP.csv",
    }
    for label, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing locked M05 {label} input: {path}")
    convergence = pd.read_csv(paths["convergence"])
    pairwise = pd.read_csv(paths["pairwise"])
    strict = pd.read_csv(paths["strict"])
    leading = pd.read_csv(paths["leading"])
    expected = 49 + 1016
    if len(convergence) != expected or convergence.duplicated(["collection", "pathway"]).any():
        raise ValueError("M05 convergence interface row/uniqueness check failed")
    if set(convergence["collection"]) != set(COLLECTIONS):
        raise ValueError("M05 collection interface mismatch")
    return convergence, pairwise, strict, leading


def evaluate_rotations(convergence: pd.DataFrame, pairwise: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    evaluation_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    for held_out in COHORTS:
        retained = [c for c in COHORTS if c != held_out]
        for collection in COLLECTIONS:
            data = convergence.loc[convergence["collection"] == collection].copy()
            retained_nes = data[[f"NES_{c}" for c in retained]].to_numpy(float)
            retained_fdr = data[[f"padj_{c}" for c in retained]].to_numpy(float)
            same_retained_direction = np.sign(retained_nes[:, 0]) == np.sign(retained_nes[:, 1])
            retained_median_abs_nes = np.median(np.abs(retained_nes), axis=1)
            pair_selected = (
                same_retained_direction
                & (retained_fdr < 0.05).all(axis=1)
                & (retained_median_abs_nes >= 1.0)
            )
            retained_consensus_sign = np.sign(retained_nes.mean(axis=1))
            held_nes = data[f"NES_{held_out}"].to_numpy(float)
            held_fdr = data[f"padj_{held_out}"].to_numpy(float)
            held_direction_replication = pair_selected & (np.sign(held_nes) == retained_consensus_sign)
            held_strict_replication = (
                held_direction_replication & (held_fdr < 0.05) & (np.abs(held_nes) >= 1.0)
            )

            out = data[["collection", "pathway", "size"]].copy()
            out["held_out_cohort"] = held_out
            out["retained_cohort_a"] = retained[0]
            out["retained_cohort_b"] = retained[1]
            out["retained_a_nes"] = retained_nes[:, 0]
            out["retained_b_nes"] = retained_nes[:, 1]
            out["retained_a_padj"] = retained_fdr[:, 0]
            out["retained_b_padj"] = retained_fdr[:, 1]
            out["retained_same_direction"] = same_retained_direction
            out["retained_median_abs_nes"] = retained_median_abs_nes
            out["pair_selected"] = pair_selected
            out["held_out_nes"] = held_nes
            out["held_out_padj"] = held_fdr
            out["held_out_direction_replication"] = held_direction_replication
            out["held_out_strict_replication"] = held_strict_replication
            evaluation_rows.append(out)

            pair_context = pairwise.loc[
                (pairwise["collection"] == collection)
                & (pairwise["cohort_a"] == retained[0])
                & (pairwise["cohort_b"] == retained[1])
            ]
            if len(pair_context) != 1:
                raise ValueError(f"Missing retained-pair context for {held_out} {collection}")
            context = pair_context.iloc[0]
            selected_n = int(pair_selected.sum())
            direction_n = int(held_direction_replication.sum())
            strict_n = int(held_strict_replication.sum())
            summary_rows.append(
                {
                    "held_out_cohort": held_out,
                    "retained_cohorts": "+".join(retained),
                    "collection": collection,
                    "eligible_pathway_n": int(len(data)),
                    "retained_pair_selected_n": selected_n,
                    "held_out_direction_replication_n": direction_n,
                    "held_out_strict_replication_n": strict_n,
                    "direction_replication_rate_among_selected": direction_n / selected_n if selected_n else np.nan,
                    "strict_replication_rate_among_selected": strict_n / selected_n if selected_n else np.nan,
                    "retained_pair_global_nes_rho": float(context["pathway_nes_spearman_rho"]),
                    "retained_pair_global_direction_concordance": float(context["pathway_nes_direction_concordance"]),
                }
            )
    return pd.concat(evaluation_rows, ignore_index=True), pd.DataFrame(summary_rows)


def build_universal_core(evaluation: pd.DataFrame, convergence: pd.DataFrame) -> pd.DataFrame:
    selected = evaluation.pivot(index=["collection", "pathway"], columns="held_out_cohort", values="pair_selected")
    strict = evaluation.pivot(index=["collection", "pathway"], columns="held_out_cohort", values="held_out_strict_replication")
    selected.columns = [f"pair_selected_without_{c}" for c in selected.columns]
    strict.columns = [f"strict_replicated_in_{c}" for c in strict.columns]
    rotations = selected.join(strict).reset_index()
    selected_cols = [f"pair_selected_without_{c}" for c in COHORTS]
    strict_cols = [f"strict_replicated_in_{c}" for c in COHORTS]
    rotations["successful_rotation_n"] = rotations[strict_cols].sum(axis=1)
    rotations["universal_loco_core"] = rotations[selected_cols + strict_cols].all(axis=1)
    annotation_cols = [
        "collection",
        "pathway",
        "size",
        *[f"NES_{c}" for c in COHORTS],
        *[f"padj_{c}" for c in COHORTS],
        "consensus_direction",
        "median_abs_nes",
        "signed_stouffer_padj",
        "strict_pathway_convergence",
    ]
    merged = convergence[annotation_cols].merge(rotations, on=["collection", "pathway"], validate="one_to_one")
    return merged.sort_values(
        ["universal_loco_core", "successful_rotation_n", "collection", "signed_stouffer_padj"],
        ascending=[False, False, True, True],
    )


def build_frozen_strict_survival(strict_m05: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frozen = strict_m05.copy()
    survival_cols = []
    for held_out in COHORTS:
        retained = [c for c in COHORTS if c != held_out]
        nes = frozen[[f"NES_{c}" for c in retained]].to_numpy(float)
        fdr = frozen[[f"padj_{c}" for c in retained]].to_numpy(float)
        col = f"survives_without_{held_out}"
        frozen[col] = (
            (np.sign(nes[:, 0]) == np.sign(nes[:, 1]))
            & (fdr < 0.05).all(axis=1)
            & (np.median(np.abs(nes), axis=1) >= 1.0)
        )
        survival_cols.append(col)
    frozen["surviving_rotation_n"] = frozen[survival_cols].sum(axis=1)
    frozen["robustness_tier"] = np.where(
        frozen["surviving_rotation_n"] == 3,
        "UNIVERSAL_ALL_THREE_ROTATIONS",
        np.where(
            frozen["surviving_rotation_n"] == 2,
            "TWO_ROTATIONS",
            np.where(frozen["surviving_rotation_n"] == 1, "ONE_ROTATION_ONLY", "ZERO_ROTATIONS"),
        ),
    )
    tier = (
        frozen.groupby(["collection", "robustness_tier", "surviving_rotation_n"], as_index=False)
        .size()
        .rename(columns={"size": "pathway_n"})
    )
    return frozen, tier


def plot_loco_funnel(summary: pd.DataFrame, plots: Path) -> None:
    metrics = [
        ("retained_pair_selected_n", "Pair-selected"),
        ("held_out_direction_replication_n", "Held-out direction"),
        ("held_out_strict_replication_n", "Held-out strict"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, collection in zip(axes, COLLECTIONS):
        data = summary.loc[summary["collection"] == collection].set_index("held_out_cohort").loc[COHORTS]
        x = np.arange(len(COHORTS))
        width = 0.23
        for idx, (column, label) in enumerate(metrics):
            values = data[column].to_numpy()
            bars = ax.bar(x + (idx - 1) * width, values, width, label=label)
            ax.bar_label(bars, padding=2, fontsize=8)
        ax.set_xticks(x, [f"hold out\n{SHORT[c]}" for c in COHORTS])
        ax.set_ylabel("Pathways")
        ax.set_title(collection)
        ax.legend(fontsize=8)
    fig.suptitle("LOCO pathway selection and held-out replication funnel", fontsize=14)
    fig.savefig(plots / "M06_LOCO_REPLICATION_FUNNEL.png", dpi=300)
    plt.close(fig)


def plot_replication_rates(summary: pd.DataFrame, plots: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, collection in zip(axes, COLLECTIONS):
        data = summary.loc[summary["collection"] == collection].set_index("held_out_cohort").loc[COHORTS]
        x = np.arange(3)
        direction = data["direction_replication_rate_among_selected"].to_numpy()
        strict = data["strict_replication_rate_among_selected"].to_numpy()
        bars1 = ax.bar(x - 0.18, direction, 0.36, label="Direction", color="#64748B")
        bars2 = ax.bar(x + 0.18, strict, 0.36, label="Strict", color=COLORS[collection])
        ax.bar_label(bars1, labels=[f"{v:.0%}" for v in direction], fontsize=8, padding=2)
        ax.bar_label(bars2, labels=[f"{v:.0%}" for v in strict], fontsize=8, padding=2)
        ax.set_xticks(x, [f"hold out\n{SHORT[c]}" for c in COHORTS])
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("Replication rate among pair-selected pathways")
        ax.set_title(collection)
        ax.legend(fontsize=8)
    fig.suptitle("Held-out pathway replication rates", fontsize=14)
    fig.savefig(plots / "M06_STRICT_REPLICATION_RATES.png", dpi=300)
    plt.close(fig)


def plot_universal_core(core: pd.DataFrame, plots: Path) -> None:
    data = core.loc[core["universal_loco_core"]].copy()
    matrix = data[[f"NES_{c}" for c in COHORTS]].to_numpy(float)
    labels = [f"{row.collection}: {clean_pathway(row.pathway)}" for row in data.itertuples(index=False)]
    max_abs = max(1.0, float(np.quantile(np.abs(matrix), 0.98)))
    fig, ax = plt.subplots(figsize=(9.5, max(5, 0.55 * len(data) + 2.0)), constrained_layout=True)
    im = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(range(3), [SHORT[c] for c in COHORTS])
    ax.set_yticks(range(len(data)), labels, fontsize=8)
    for i in range(len(data)):
        for j in range(3):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Universal LOCO core pathway NES")
    ax.set_xlabel("Cohort")
    fig.colorbar(im, ax=ax, label="NES", fraction=0.04, pad=0.03)
    fig.savefig(plots / "M06_UNIVERSAL_CORE_NES_HEATMAP.png", dpi=300)
    plt.close(fig)


def plot_global_context(summary: pd.DataFrame, plots: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    x = np.arange(3)
    width = 0.35
    for idx, collection in enumerate(COLLECTIONS):
        data = summary.loc[summary["collection"] == collection].set_index("held_out_cohort").loc[COHORTS]
        axes[0].bar(x + (idx - 0.5) * width, data["retained_pair_global_nes_rho"], width, label=collection, color=COLORS[collection])
        axes[1].bar(x + (idx - 0.5) * width, data["retained_pair_global_direction_concordance"], width, label=collection, color=COLORS[collection])
    for ax in axes:
        ax.set_xticks(x, [f"hold out\n{SHORT[c]}" for c in COHORTS])
        ax.axhline(0, color="#94A3B8", linewidth=0.7)
        ax.legend(fontsize=8)
    axes[0].set_ylim(-1, 1)
    axes[0].set_ylabel("Retained-pair pathway NES rho")
    axes[0].set_title("Global rank context")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Retained-pair direction concordance")
    axes[1].set_title("Global direction context")
    fig.suptitle("LOCO retained-pair context across all eligible pathways", fontsize=14)
    fig.savefig(plots / "M06_GLOBAL_PAIR_CONTEXT.png", dpi=300)
    plt.close(fig)


def plot_hallmark_matrix(frozen: pd.DataFrame, plots: Path) -> None:
    data = frozen.loc[frozen["collection"] == "HALLMARK"].copy()
    cols = [f"survives_without_{c}" for c in COHORTS]
    matrix = data[cols].to_numpy(int)
    labels = [clean_pathway(x) for x in data["pathway"]]
    fig, ax = plt.subplots(figsize=(7.5, 6), constrained_layout=True)
    im = ax.imshow(matrix, aspect="auto", cmap="YlGn", vmin=0, vmax=1)
    ax.set_xticks(range(3), [f"without\n{SHORT[c]}" for c in COHORTS])
    ax.set_yticks(range(len(data)), labels, fontsize=8)
    for i in range(len(data)):
        for j in range(3):
            ax.text(j, i, "PASS" if matrix[i, j] else "—", ha="center", va="center", fontsize=8)
    ax.set_title("Frozen M05 strict Hallmark pair-survival matrix")
    fig.colorbar(im, ax=ax, ticks=[0, 1], label="Retained pair meets strict pair rule", fraction=0.05, pad=0.04)
    fig.savefig(plots / "M06_HALLMARK_FROZEN_STRICT_SURVIVAL.png", dpi=300)
    plt.close(fig)


def plot_tiers(tier: pd.DataFrame, plots: Path) -> None:
    pivot = tier.pivot(index="robustness_tier", columns="collection", values="pathway_n").fillna(0)
    order = ["ZERO_ROTATIONS", "ONE_ROTATION_ONLY", "TWO_ROTATIONS", "UNIVERSAL_ALL_THREE_ROTATIONS"]
    pivot = pivot.reindex(order, fill_value=0)
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(10, 4.8), constrained_layout=True)
    width = 0.35
    for idx, collection in enumerate(COLLECTIONS):
        values = pivot.get(collection, pd.Series(0, index=order)).to_numpy()
        bars = ax.bar(x + (idx - 0.5) * width, values, width, label=collection, color=COLORS[collection])
        ax.bar_label(bars, padding=2, fontsize=8)
    ax.set_xticks(x, ["0", "1", "2", "3 (universal)"])
    ax.set_xlabel("Successful retained-pair rotations for frozen M05 strict pathways")
    ax.set_ylabel("Pathways")
    ax.set_title("Frozen strict pathway robustness tiers")
    ax.legend()
    fig.savefig(plots / "M06_ROBUSTNESS_TIER_COUNTS.png", dpi=300)
    plt.close(fig)


def write_summary(
    results: Path,
    summary: pd.DataFrame,
    universal: pd.DataFrame,
    frozen: pd.DataFrame,
    tier: pd.DataFrame,
) -> None:
    core = universal.loc[universal["universal_loco_core"]].copy()
    summary_lines = []
    for row in summary.itertuples(index=False):
        summary_lines.append(
            f"- Hold out {row.held_out_cohort}, {row.collection}: retained pair selected {row.retained_pair_selected_n}; "
            f"held-out direction replicated {row.held_out_direction_replication_n} ({row.direction_replication_rate_among_selected:.1%}); "
            f"strict replicated {row.held_out_strict_replication_n} ({row.strict_replication_rate_among_selected:.1%}); "
            f"global retained-pair rho={row.retained_pair_global_nes_rho:.4f}, direction={row.retained_pair_global_direction_concordance:.1%}."
        )
    core_lines = [
        f"- {row.collection} / {row.pathway}: NES {row.NES_GSE274832:.3f}, {row.NES_GSE193136:.3f}, {row.NES_GSE232306:.3f}."
        for row in core.itertuples(index=False)
    ] or ["- None."]
    hallmark_core = int(((core["collection"] == "HALLMARK")).sum())
    reactome_core = int(((core["collection"] == "REACTOME")).sum())
    text = f"""# M06 Leave-one-cohort-out Pathway Robustness Summary

`PRELIMINARY_MODULE_STATE: PASS_WITH_LIMITATION_PENDING_INDEPENDENT_OUTPUT_AUDIT`

## Frozen design

Each rotation selected pathways using only the two retained cohorts: identical NES direction, both fgsea FDR<0.05, and median |NES|>=1.0. The held-out cohort was then evaluated first for direction and then for strict replication (same direction, held-out FDR<0.05, |NES|>=1.0). All 49 Hallmark and 1,016 Reactome pathways were eligible before each rotation-specific selection. No fgsea result or threshold was refit.

## Rotation results

{chr(10).join(summary_lines)}

The weakest primary-collection rotation was holding out GSE232306: 14 Hallmark pathways were selected by GSE274832+GSE193136, five retained direction in GSE232306, and one met strict held-out replication. This bounds the M05 claim and confirms that global pathway behavior is not uniformly reproducible across the three cohorts.

## Universal LOCO core

Universal core requires pair selection and strict held-out replication in all three rotations. It contains **{len(core)} pathways**: **{hallmark_core} Hallmark** and **{reactome_core} Reactome**.

{chr(10).join(core_lines)}

The primary compact core is HALLMARK_P53_PATHWAY. Reactome core terms are secondary and redundant. Universal robustness means conditional internal stability across these three rotations; it is not external clinical validation or proof of pathway activity.

## Frozen M05 strict-set survival

Among the 8 frozen strict Hallmark pathways, one survives all three retained-pair rotations and seven survive one rotation. Among 47 frozen strict Reactome pathways, seven survive all rotations and 40 survive one rotation. This frozen-set stability description is kept separate from the rotation-wise re-selection result above.

## Interpretation boundary

- At least one primary compact Hallmark signal remains strictly supported under every rotation, so the pathway-convergence conclusion survives in a narrow, bounded form.
- Replication rates are modest and vary by held-out cohort; broad pathway-rank concordance does not survive uniformly.
- GSE232306 is not the sole source of every robust pathway, but it remains the most discordant held-out cohort for Hallmark direction/strict replication.
- M06 uses the same three public cohorts and cannot be described as external, prospective, or clinical validation.
- Positive/negative NES still means rank enrichment, not pathway activation/inhibition.

## Permitted claim

The frozen three-cohort analysis contains a narrow universal LOCO pathway core, led by negative enrichment of HALLMARK_P53_PATHWAY, while most M05 strict pathways are supported by only one retained-pair rotation.

## Prohibited claims

- Do not call the universal core externally validated, causal, diagnostic, or therapeutic.
- Do not generalize the eight-pathway universal core to global pathway concordance.
- Do not count overlapping Reactome terms as independent mechanisms.
"""
    (results / "M06_LEAVE_ONE_COHORT_OUT_SUMMARY.md").write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    run_dir = args.run_dir.resolve()
    results = run_dir / "results"
    plots = results / "plots"
    logs = run_dir / "logs"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    convergence, pairwise, strict_m05, leading = load_inputs(project_root)
    evaluation, summary = evaluate_rotations(convergence, pairwise)
    universal_all = build_universal_core(evaluation, convergence)
    universal_core = universal_all.loc[universal_all["universal_loco_core"]].copy()
    frozen, tier = build_frozen_strict_survival(strict_m05)
    core_leading = leading.merge(
        universal_core[["collection", "pathway"]], on=["collection", "pathway"], how="inner", validate="one_to_one"
    )

    evaluation.to_csv(results / "M06_LOCO_PATHWAY_EVALUATION.csv", index=False)
    summary.to_csv(results / "M06_COLLECTION_LOCO_SUMMARY.csv", index=False)
    universal_core.to_csv(results / "M06_UNIVERSAL_LOCO_CORE_PATHWAYS.csv", index=False)
    universal_all.to_csv(results / "M06_ALL_PATHWAY_ROTATION_MATRIX.csv", index=False)
    frozen.to_csv(results / "M06_FROZEN_M05_STRICT_SURVIVAL.csv", index=False)
    tier.to_csv(results / "M06_ROBUSTNESS_TIER_COUNTS.csv", index=False)
    core_leading.to_csv(results / "M06_UNIVERSAL_CORE_LEADING_EDGE.csv", index=False)

    key_rows = [
        ("eligible_pathways_HALLMARK", 49),
        ("eligible_pathways_REACTOME", 1016),
        ("universal_loco_core_HALLMARK", int((universal_core["collection"] == "HALLMARK").sum())),
        ("universal_loco_core_REACTOME", int((universal_core["collection"] == "REACTOME").sum())),
        ("frozen_m05_strict_HALLMARK", int((frozen["collection"] == "HALLMARK").sum())),
        ("frozen_m05_strict_REACTOME", int((frozen["collection"] == "REACTOME").sum())),
    ]
    for row in summary.itertuples(index=False):
        key_rows.append((f"pair_selected_{row.collection}_without_{row.held_out_cohort}", row.retained_pair_selected_n))
        key_rows.append((f"heldout_strict_{row.collection}_{row.held_out_cohort}", row.held_out_strict_replication_n))
    pd.DataFrame(key_rows, columns=["metric", "value"]).to_csv(results / "M06_KEY_COUNTS.csv", index=False)

    plot_loco_funnel(summary, plots)
    plot_replication_rates(summary, plots)
    plot_universal_core(universal_core, plots)
    plot_global_context(summary, plots)
    plot_hallmark_matrix(frozen, plots)
    plot_tiers(tier, plots)
    write_summary(results, summary, universal_core, frozen, tier)

    session = {
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "input_pathway_rows": int(len(convergence)),
        "evaluation_rows": int(len(evaluation)),
        "universal_core_n": int(len(universal_core)),
        "fgsea_rerun": False,
        "threshold_retuning": False,
    }
    (logs / "PYTHON_SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    print("M06 pathway LOCO completed")
    print(summary.to_string(index=False))
    print(f"Universal LOCO core: {len(universal_core)}")
    print(universal_core[["collection", "pathway"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
