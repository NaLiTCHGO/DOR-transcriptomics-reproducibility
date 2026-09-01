#!/usr/bin/env python3
"""M04 cross-cohort DOR effect reproducibility and heterogeneity synthesis.

This script reads only frozen cohort-level M03 effect tables. It never reads,
combines, or writes a sample-level expression matrix.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm, spearmanr


COHORTS = ["GSE274832", "GSE193136", "GSE232306"]
COHORT_SHORT = {
    "GSE274832": "274832",
    "GSE193136": "193136",
    "GSE232306": "232306",
}
COLORS = {
    "GSE274832": "#3B82F6",
    "GSE193136": "#10B981",
    "GSE232306": "#F59E0B",
}
INPUTS = {
    "GSE274832": {
        "path": "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE274832/results/GSE274832_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
        "model": "condition-only DOR-minus-NOR",
        "known_limitation": "n=3/group; metadata template conflict disclosed",
    },
    "GSE193136": {
        "path": "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE193136/results/GSE193136_DESEQ2_AGE_ADJUSTED_DOR_MINUS_NOR_ALL_GENES.csv",
        "model": "age-adjusted DOR-minus-NOR primary",
        "known_limitation": "low transcriptome mapping; age imbalance modeled",
    },
    "GSE232306": {
        "path": "06_locked_results/modules/M03_WITHIN_COHORT_EFFECTS/v1_GSE232306/results/GSE232306_DESEQ2_DOR_MINUS_NOR_ALL_GENES.csv",
        "model": "condition-only DOR-minus-NOR; age unavailable",
        "known_limitation": "phenotype-aligned GC/duplication/PC1; unmeasured confounding possible",
    },
}
TOP_K = [100, 500, 1000]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    return parser.parse_args()


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return adjusted
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    adjusted[valid] = restored
    return adjusted


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = spearmanr(x, y, nan_policy="omit")
    return float(result.statistic), float(result.pvalue)


def load_effect_table(project_root: Path, cohort: str) -> tuple[pd.DataFrame, dict]:
    source = project_root / INPUTS[cohort]["path"]
    if not source.is_file():
        raise FileNotFoundError(f"Missing locked M03 input: {source}")
    raw = pd.read_csv(source)
    required = {
        "gene_id",
        "gene_name",
        "log2FoldChange",
        "lfcSE",
        "pvalue",
        "padj",
    }
    missing = sorted(required.difference(raw.columns))
    if missing:
        raise ValueError(f"{cohort} input missing columns: {missing}")
    if raw["gene_id"].duplicated().any():
        raise ValueError(f"{cohort} contains duplicate gene_id rows")
    if not raw["gene_id"].astype(str).str.match(r"^ENSG[0-9]+$").all():
        raise ValueError(f"{cohort} gene_id namespace is not unversioned Ensembl ENSG")

    numeric = raw[["log2FoldChange", "lfcSE", "pvalue", "padj"]].apply(
        pd.to_numeric, errors="coerce"
    )
    valid = (
        np.isfinite(numeric["log2FoldChange"])
        & np.isfinite(numeric["lfcSE"])
        & (numeric["lfcSE"] > 0)
    )
    table = raw.loc[valid, ["gene_id", "gene_name"]].copy()
    table[f"beta_{cohort}"] = numeric.loc[valid, "log2FoldChange"].to_numpy()
    table[f"se_{cohort}"] = numeric.loc[valid, "lfcSE"].to_numpy()
    table[f"pvalue_{cohort}"] = numeric.loc[valid, "pvalue"].to_numpy()
    table[f"padj_{cohort}"] = numeric.loc[valid, "padj"].to_numpy()
    table[f"z_{cohort}"] = table[f"beta_{cohort}"] / table[f"se_{cohort}"]
    table = table.rename(columns={"gene_name": f"gene_name_{cohort}"})

    audit = {
        "cohort": cohort,
        "input_path": INPUTS[cohort]["path"],
        "effect_model": INPUTS[cohort]["model"],
        "known_limitation": INPUTS[cohort]["known_limitation"],
        "rows_total": int(len(raw)),
        "rows_valid_effect_se": int(valid.sum()),
        "duplicate_gene_ids": int(raw["gene_id"].duplicated().sum()),
        "nonpositive_or_missing_se": int((~np.isfinite(numeric["lfcSE"]) | (numeric["lfcSE"] <= 0)).sum()),
        "comparison": "DOR-minus-NOR",
        "namespace": "Ensembl gene ID without version",
    }
    return table, audit


def resolve_gene_names(common: pd.DataFrame) -> pd.Series:
    name_cols = [f"gene_name_{cohort}" for cohort in COHORTS]
    names = common[name_cols].copy()
    for column in name_cols:
        as_text = names[column].astype("string")
        names[column] = as_text.mask(as_text.str.match(r"^ENSG[0-9]+$", na=False))
    resolved = names.bfill(axis=1).iloc[:, 0]
    return resolved.fillna(common["gene_id"]).astype(str)


def build_pairwise_results(common: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_rows: list[dict] = []
    topk_rows: list[dict] = []
    for cohort_a, cohort_b in combinations(COHORTS, 2):
        beta_a = common[f"beta_{cohort_a}"].to_numpy()
        beta_b = common[f"beta_{cohort_b}"].to_numpy()
        rho, pvalue = safe_spearman(beta_a, beta_b)
        nonzero = (beta_a != 0) & (beta_b != 0)
        same = np.sign(beta_a[nonzero]) == np.sign(beta_b[nonzero])
        cand_a = common[f"candidate_{cohort_a}"].to_numpy(dtype=bool)
        cand_b = common[f"candidate_{cohort_b}"].to_numpy(dtype=bool)
        union = cand_a | cand_b
        overlap = cand_a & cand_b
        overlap_direction = (
            float((np.sign(beta_a[overlap]) == np.sign(beta_b[overlap])).mean())
            if overlap.any()
            else np.nan
        )
        pair_rows.append(
            {
                "cohort_a": cohort_a,
                "cohort_b": cohort_b,
                "n_common_genes": int(len(common)),
                "spearman_rho_log2fc": rho,
                "spearman_pvalue": pvalue,
                "direction_concordance_all_nonzero": float(same.mean()),
                "n_candidate_a": int(cand_a.sum()),
                "n_candidate_b": int(cand_b.sum()),
                "candidate_overlap_n": int(overlap.sum()),
                "candidate_union_n": int(union.sum()),
                "candidate_jaccard": float(overlap.sum() / union.sum()) if union.any() else np.nan,
                "candidate_overlap_direction_concordance": overlap_direction,
            }
        )

        for k in TOP_K:
            top_a = set(
                common.nlargest(k, f"abs_z_{cohort_a}")["gene_id"].astype(str)
            )
            top_b = set(
                common.nlargest(k, f"abs_z_{cohort_b}")["gene_id"].astype(str)
            )
            genes = top_a.intersection(top_b)
            union_top = top_a.union(top_b)
            if genes:
                subset = common.set_index("gene_id").loc[sorted(genes)]
                same_top_direction = float(
                    (
                        np.sign(subset[f"beta_{cohort_a}"])
                        == np.sign(subset[f"beta_{cohort_b}"])
                    ).mean()
                )
            else:
                same_top_direction = np.nan
            topk_rows.append(
                {
                    "cohort_a": cohort_a,
                    "cohort_b": cohort_b,
                    "top_k_by_abs_wald_z": k,
                    "overlap_n": len(genes),
                    "jaccard": len(genes) / len(union_top),
                    "overlap_direction_concordance": same_top_direction,
                }
            )
    return pd.DataFrame(pair_rows), pd.DataFrame(topk_rows)


def add_meta_analysis(common: pd.DataFrame) -> pd.DataFrame:
    betas = common[[f"beta_{c}" for c in COHORTS]].to_numpy(dtype=float)
    ses = common[[f"se_{c}" for c in COHORTS]].to_numpy(dtype=float)
    weights = 1.0 / np.square(ses)
    sum_w = weights.sum(axis=1)
    fixed_beta = (weights * betas).sum(axis=1) / sum_w
    fixed_se = np.sqrt(1.0 / sum_w)
    fixed_z = fixed_beta / fixed_se
    fixed_p = 2.0 * norm.sf(np.abs(fixed_z))

    q = (weights * np.square(betas - fixed_beta[:, None])).sum(axis=1)
    df = len(COHORTS) - 1
    c_term = sum_w - np.square(weights).sum(axis=1) / sum_w
    tau2 = np.maximum(0.0, (q - df) / c_term)
    random_weights = 1.0 / (np.square(ses) + tau2[:, None])
    random_sum_w = random_weights.sum(axis=1)
    random_beta = (random_weights * betas).sum(axis=1) / random_sum_w
    random_se = np.sqrt(1.0 / random_sum_w)
    random_z = random_beta / random_se
    random_p = 2.0 * norm.sf(np.abs(random_z))
    i2 = np.where(q > 0, np.maximum(0.0, (q - df) / q) * 100.0, 0.0)
    q_p = scipy.stats.chi2.sf(q, df)

    candidate_matrix = common[[f"candidate_{c}" for c in COHORTS]].to_numpy(dtype=bool)
    fdr_matrix = common[[f"fdr_{c}" for c in COHORTS]].to_numpy(dtype=bool)
    positive = (betas > 0).sum(axis=1)
    negative = (betas < 0).sum(axis=1)
    candidate_n = candidate_matrix.sum(axis=1)
    candidate_positive = ((betas > 0) & candidate_matrix).sum(axis=1)
    candidate_negative = ((betas < 0) & candidate_matrix).sum(axis=1)
    all_same_direction = (positive == len(COHORTS)) | (negative == len(COHORTS))
    candidate_2plus_same_direction = (candidate_n >= 2) & (
        (candidate_positive == candidate_n) | (candidate_negative == candidate_n)
    )

    result = common.copy()
    result["positive_cohort_n"] = positive
    result["negative_cohort_n"] = negative
    result["all_three_same_direction"] = all_same_direction
    result["consensus_direction"] = np.where(
        positive > negative, "UP_IN_DOR", np.where(negative > positive, "DOWN_IN_DOR", "TIE")
    )
    result["fdr_lt_0_05_cohort_n"] = fdr_matrix.sum(axis=1)
    result["fdr_lfc1_candidate_cohort_n"] = candidate_n
    result["candidate_2plus_same_candidate_direction"] = candidate_2plus_same_direction
    result["fixed_beta"] = fixed_beta
    result["fixed_se"] = fixed_se
    result["fixed_z"] = fixed_z
    result["fixed_pvalue"] = fixed_p
    result["fixed_padj"] = bh_adjust(fixed_p)
    result["cochran_q"] = q
    result["cochran_q_df"] = df
    result["cochran_q_pvalue"] = q_p
    result["tau2_dl"] = tau2
    result["i2_percent"] = i2
    result["random_beta"] = random_beta
    result["random_se"] = random_se
    result["random_z"] = random_z
    result["random_pvalue"] = random_p
    result["random_padj"] = bh_adjust(random_p)
    result["strict_exploratory_consensus"] = (
        (result["random_padj"] < 0.05)
        & (result["random_beta"].abs() >= 0.5)
        & result["all_three_same_direction"]
        & (result["fdr_lfc1_candidate_cohort_n"] >= 2)
        & (result["i2_percent"] < 50.0)
    )
    result["heterogeneity_category"] = pd.cut(
        result["i2_percent"],
        bins=[-np.inf, 25.0, 50.0, 75.0, np.inf],
        labels=["LOW_LT25", "MODERATE_25_TO_LT50", "SUBSTANTIAL_50_TO_LT75", "HIGH_GE75"],
        right=False,
    ).astype(str)
    return result


def build_overlap_summary(meta: pd.DataFrame) -> pd.DataFrame:
    flags = meta[[f"candidate_{c}" for c in COHORTS]].copy()
    flags.columns = COHORTS
    labels = []
    for row in flags.itertuples(index=False, name=None):
        members = [cohort for cohort, present in zip(COHORTS, row) if present]
        labels.append("+".join(members) if members else "NONE")
    summary = pd.Series(labels, name="candidate_membership").value_counts().rename_axis(
        "candidate_membership"
    ).reset_index(name="gene_n")
    summary["cohort_n"] = summary["candidate_membership"].apply(
        lambda x: 0 if x == "NONE" else x.count("+") + 1
    )
    return summary.sort_values(["cohort_n", "gene_n"], ascending=[False, False])


def build_loco(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_level = meta[["gene_id", "gene_name", "random_beta", "strict_exploratory_consensus"]].copy()
    rows: list[dict] = []
    strict = meta["strict_exploratory_consensus"].to_numpy(dtype=bool)
    all_three_same = meta["all_three_same_direction"].to_numpy(dtype=bool)
    full = meta["random_beta"].to_numpy(dtype=float)
    for held_out in COHORTS:
        retained = [c for c in COHORTS if c != held_out]
        betas = meta[[f"beta_{c}" for c in retained]].to_numpy(dtype=float)
        ses = meta[[f"se_{c}" for c in retained]].to_numpy(dtype=float)
        weights = 1.0 / np.square(ses)
        loo_beta = (weights * betas).sum(axis=1) / weights.sum(axis=1)
        loo_se = np.sqrt(1.0 / weights.sum(axis=1))
        rho, pvalue = safe_spearman(full, loo_beta)
        direction = np.sign(full) == np.sign(loo_beta)
        strict_retention = float(direction[strict].mean()) if strict.any() else np.nan
        all_same_retention = (
            float(direction[all_three_same].mean()) if all_three_same.any() else np.nan
        )
        prefix = f"without_{held_out}"
        gene_level[f"{prefix}_fixed_beta"] = loo_beta
        gene_level[f"{prefix}_fixed_se"] = loo_se
        gene_level[f"{prefix}_same_direction_as_full_random"] = direction
        rows.append(
            {
                "held_out_cohort": held_out,
                "retained_cohorts": "+".join(retained),
                "n_genes": int(len(meta)),
                "spearman_vs_full_random_beta": rho,
                "spearman_pvalue": pvalue,
                "direction_agreement_vs_full_random": float(direction.mean()),
                "median_abs_beta_shift": float(np.median(np.abs(loo_beta - full))),
                "all_three_direction_consensus_retention": all_same_retention,
                "strict_consensus_direction_retention": strict_retention,
                "strict_consensus_n": int(strict.sum()),
            }
        )
    return pd.DataFrame(rows), gene_level


def plot_pairwise_scatter(meta: pd.DataFrame, pairwise: pd.DataFrame, plots: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for ax, (cohort_a, cohort_b) in zip(axes, combinations(COHORTS, 2)):
        x = meta[f"beta_{cohort_a}"].to_numpy()
        y = meta[f"beta_{cohort_b}"].to_numpy()
        limit = max(1.0, float(np.quantile(np.abs(np.concatenate([x, y])), 0.99)))
        ax.scatter(x, y, s=5, alpha=0.12, color="#334155", linewidths=0, rasterized=True)
        ax.plot([-limit, limit], [-limit, limit], linestyle="--", color="#94A3B8", linewidth=1)
        ax.axhline(0, color="#CBD5E1", linewidth=0.8)
        ax.axvline(0, color="#CBD5E1", linewidth=0.8)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xlabel(f"{cohort_a} log2FC")
        ax.set_ylabel(f"{cohort_b} log2FC")
        row = pairwise[(pairwise.cohort_a == cohort_a) & (pairwise.cohort_b == cohort_b)].iloc[0]
        ax.set_title(
            f"rho={row.spearman_rho_log2fc:.3f}; direction={row.direction_concordance_all_nonzero:.1%}\n"
            "display clipped at pooled 99th |log2FC| percentile",
            fontsize=10,
        )
    fig.suptitle("Cross-cohort DOR-minus-NOR effect reproducibility", fontsize=14)
    fig.savefig(plots / "M04_PAIRWISE_EFFECT_SCATTER.png", dpi=300)
    plt.close(fig)


def plot_reproducibility_matrices(pairwise: pd.DataFrame, plots: Path) -> None:
    n = len(COHORTS)
    rho = np.eye(n)
    direction = np.ones((n, n))
    for row in pairwise.itertuples(index=False):
        i = COHORTS.index(row.cohort_a)
        j = COHORTS.index(row.cohort_b)
        rho[i, j] = rho[j, i] = row.spearman_rho_log2fc
        direction[i, j] = direction[j, i] = row.direction_concordance_all_nonzero

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), constrained_layout=True)
    configs = [
        (rho, "Spearman rho of log2FC", "coolwarm", -1.0, 1.0, ".3f"),
        (direction, "Effect-direction concordance", "YlGnBu", 0.0, 1.0, ".1%"),
    ]
    labels = [COHORT_SHORT[c] for c in COHORTS]
    for ax, (matrix, title, cmap, vmin, vmax, fmt) in zip(axes, configs):
        image = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(n), labels=labels)
        ax.set_yticks(range(n), labels=labels)
        ax.set_title(title)
        for i in range(n):
            for j in range(n):
                text = format(matrix[i, j], fmt)
                ax.text(j, i, text, ha="center", va="center", fontsize=10)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Three-core-cohort reproducibility matrix", fontsize=14)
    fig.savefig(plots / "M04_REPRODUCIBILITY_MATRICES.png", dpi=300)
    plt.close(fig)


def plot_candidate_overlap(overlap: pd.DataFrame, plots: Path) -> None:
    shown = overlap[overlap["candidate_membership"] != "NONE"].copy()
    shown = shown.sort_values("gene_n", ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.5 * len(shown) + 1.5)), constrained_layout=True)
    y = np.arange(len(shown))
    values = shown["gene_n"].to_numpy(dtype=float)
    ax.hlines(y, 0.8, values, color="#A5B4FC", linewidth=3)
    ax.scatter(values, y, color="#4F46E5", s=55, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.8, values.max() * 1.65)
    ax.set_yticks(y, labels=shown["candidate_membership"])
    for i, value in enumerate(values):
        ax.text(value * 1.08, i, f"{int(value):,}", va="center", fontsize=9)
    ax.grid(axis="x", which="both", linestyle=":", color="#CBD5E1", linewidth=0.8)
    ax.set_xlabel("Genes meeting within-cohort FDR<0.05 and |log2FC|>=1 (log scale)")
    ax.set_ylabel("Candidate membership combination")
    ax.set_title("Thresholded candidate overlap (descriptive; not biomarker validation)")
    fig.savefig(plots / "M04_CANDIDATE_OVERLAP.png", dpi=300)
    plt.close(fig)


def plot_heterogeneity(meta: pd.DataFrame, plots: Path) -> None:
    order = ["LOW_LT25", "MODERATE_25_TO_LT50", "SUBSTANTIAL_50_TO_LT75", "HIGH_GE75"]
    counts = meta["heterogeneity_category"].value_counts().reindex(order, fill_value=0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].hist(meta["i2_percent"], bins=np.linspace(0, 100, 21), color="#0EA5E9", edgecolor="white")
    axes[0].set_xlabel("I2 (%)")
    axes[0].set_ylabel("Genes")
    axes[0].set_title("Gene-level I2 distribution (k=3; descriptive)")
    axes[1].bar(range(len(order)), counts.values, color=["#10B981", "#84CC16", "#F59E0B", "#EF4444"])
    axes[1].set_xticks(range(len(order)), labels=["<25", "25–<50", "50–<75", ">=75"])
    axes[1].set_xlabel("I2 category (%)")
    axes[1].set_ylabel("Genes")
    axes[1].set_title("Heterogeneity categories")
    for i, value in enumerate(counts.values):
        axes[1].text(i, value, f"{int(value):,}", ha="center", va="bottom", fontsize=9)
    fig.savefig(plots / "M04_HETEROGENEITY_DISTRIBUTION.png", dpi=300)
    plt.close(fig)


def plot_effect_heatmap(meta: pd.DataFrame, plots: Path) -> tuple[str, int]:
    strict = meta[meta["strict_exploratory_consensus"]].sort_values(
        ["random_padj", "i2_percent"], ascending=[True, True]
    )
    if len(strict) > 0:
        selected = strict.head(30).copy()
        selection = "strict exploratory consensus"
    else:
        selected = meta.sort_values(
            ["fdr_lfc1_candidate_cohort_n", "all_three_same_direction", "random_padj", "i2_percent"],
            ascending=[False, False, True, True],
        ).head(30).copy()
        selection = "fallback: strongest cross-cohort descriptive genes; strict set was empty"
    matrix = selected[[f"beta_{c}" for c in COHORTS]].to_numpy(dtype=float)
    max_abs = max(0.5, float(np.quantile(np.abs(matrix), 0.95)))
    labels = [
        name if name and name != gid else gid
        for name, gid in zip(selected["gene_name"], selected["gene_id"])
    ]
    fig_height = max(6.0, 0.28 * len(selected) + 2.4)
    fig, ax = plt.subplots(figsize=(7.2, fig_height), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(range(len(COHORTS)), labels=[COHORT_SHORT[c] for c in COHORTS])
    ax.set_yticks(range(len(selected)), labels=labels, fontsize=8)
    ax.set_xlabel("Cohort")
    ax.set_ylabel("Gene")
    ax.set_title(f"Top cross-cohort effects\nSelection: {selection}", fontsize=11)
    for i in range(len(selected)):
        for j in range(len(COHORTS)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
    fig.colorbar(image, ax=ax, label="DOR-minus-NOR log2FC", fraction=0.05, pad=0.03)
    fig.savefig(plots / "M04_TOP_CROSS_COHORT_EFFECT_HEATMAP.png", dpi=300)
    plt.close(fig)
    return selection, int(len(selected))


def plot_loco(loco: pd.DataFrame, plots: Path) -> None:
    labels = [f"without\n{COHORT_SHORT[c]}" for c in loco["held_out_cohort"]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    axes[0].bar(labels, loco["spearman_vs_full_random_beta"], color="#8B5CF6")
    axes[0].set_ylim(-1, 1)
    axes[0].axhline(0, color="#94A3B8", linewidth=0.8)
    axes[0].set_ylabel("Spearman rho")
    axes[0].set_title("LOCO effect rank vs full random-effects estimate")
    axes[1].bar(labels, loco["direction_agreement_vs_full_random"], color="#14B8A6")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Direction agreement")
    axes[1].set_title("LOCO direction vs full random-effects estimate")
    for ax in axes:
        for patch in ax.patches:
            height = patch.get_height()
            ax.text(patch.get_x() + patch.get_width() / 2, height, f"{height:.3f}", ha="center", va="bottom", fontsize=9)
    fig.savefig(plots / "M04_LOCO_STABILITY.png", dpi=300)
    plt.close(fig)


def correlation_label(value: float) -> str:
    magnitude = abs(value)
    if magnitude < 0.1:
        strength = "negligible"
    elif magnitude < 0.3:
        strength = "weak"
    elif magnitude < 0.5:
        strength = "moderate"
    else:
        strength = "strong"
    return f"{strength} {'positive' if value >= 0 else 'negative'}"


def write_summary(
    results: Path,
    meta: pd.DataFrame,
    pairwise: pd.DataFrame,
    topk: pd.DataFrame,
    overlap: pd.DataFrame,
    loco: pd.DataFrame,
    selection_label: str,
) -> pd.DataFrame:
    candidate_counts = {
        cohort: int(meta[f"candidate_{cohort}"].sum()) for cohort in COHORTS
    }
    all_three_candidates = int((meta["fdr_lfc1_candidate_cohort_n"] == 3).sum())
    two_plus_same = int(meta["candidate_2plus_same_candidate_direction"].sum())
    all_three_direction = int(meta["all_three_same_direction"].sum())
    strict_n = int(meta["strict_exploratory_consensus"].sum())
    random_fdr_n = int((meta["random_padj"] < 0.05).sum())
    high_i2_n = int((meta["i2_percent"] >= 75).sum())
    low_i2_n = int((meta["i2_percent"] < 25).sum())
    rho_min = float(pairwise["spearman_rho_log2fc"].min())
    rho_max = float(pairwise["spearman_rho_log2fc"].max())
    dir_min = float(pairwise["direction_concordance_all_nonzero"].min())
    dir_max = float(pairwise["direction_concordance_all_nonzero"].max())
    loco_rho_min = float(loco["spearman_vs_full_random_beta"].min())
    loco_rho_max = float(loco["spearman_vs_full_random_beta"].max())
    loco_dir_min = float(loco["direction_agreement_vs_full_random"].min())
    loco_dir_max = float(loco["direction_agreement_vs_full_random"].max())

    key_rows = [
        ("common_gene_universe", len(meta)),
        ("all_three_same_direction", all_three_direction),
        ("candidate_in_all_three", all_three_candidates),
        ("candidate_in_2plus_same_candidate_direction", two_plus_same),
        ("random_effects_fdr_lt_0_05", random_fdr_n),
        ("strict_exploratory_consensus", strict_n),
        ("i2_lt_25", low_i2_n),
        ("i2_ge_75", high_i2_n),
    ] + [(f"candidate_{cohort}", count) for cohort, count in candidate_counts.items()]
    key_counts = pd.DataFrame(key_rows, columns=["metric", "value"])
    key_counts.to_csv(results / "M04_KEY_COUNTS.csv", index=False)

    pair_lines = []
    for row in pairwise.itertuples(index=False):
        pair_lines.append(
            f"- {row.cohort_a} vs {row.cohort_b}: Spearman rho={row.spearman_rho_log2fc:.4f} "
            f"({correlation_label(row.spearman_rho_log2fc)}), direction concordance={row.direction_concordance_all_nonzero:.2%}, "
            f"candidate overlap={row.candidate_overlap_n:,}, candidate Jaccard={row.candidate_jaccard:.4f}."
        )
    top500 = topk[topk["top_k_by_abs_wald_z"] == 500]
    top_lines = [
        f"- {row.cohort_a} vs {row.cohort_b}: top-500 overlap={row.overlap_n:,}, "
        f"Jaccard={row.jaccard:.4f}, overlap direction concordance={row.overlap_direction_concordance:.2%}."
        for row in top500.itertuples(index=False)
    ]
    loco_lines = [
        f"- Hold out {row.held_out_cohort}: rho versus full random-effects beta={row.spearman_vs_full_random_beta:.4f}; "
        f"direction agreement={row.direction_agreement_vs_full_random:.2%}; median absolute beta shift={row.median_abs_beta_shift:.4f}."
        for row in loco.itertuples(index=False)
    ]

    summary = f"""# M04 Cross-cohort Reproducibility and Heterogeneity Summary

`PRELIMINARY_MODULE_STATE: PASS_WITH_LIMITATION_PENDING_INDEPENDENT_OUTPUT_AUDIT`

## Scope and frozen inputs

Three locked, independent human granulosa-cell RNA-seq cohort effect tables were synthesized. Every input contrast is DOR-minus-NOR. GSE193136 uses its frozen age-adjusted primary model; GSE274832 and GSE232306 use condition-only models. The synthesis key is unversioned Ensembl gene ID. No sample-level expression matrix, cross-cohort normalization, batch correction, or megamatrix was created.

The common valid universe contains **{len(meta):,} genes**.

## Primary reproducibility results

Pairwise all-gene effect-rank correlations ranged from **{rho_min:.4f} to {rho_max:.4f}**, and direction concordance ranged from **{dir_min:.2%} to {dir_max:.2%}**.

{chr(10).join(pair_lines)}

Top-500 absolute Wald-statistic overlap:

{chr(10).join(top_lines)}

These are outcome-neutral reproducibility measurements. Weak or negative concordance is a scientific heterogeneity result, not a technical pipeline failure.

## Thresholded candidate overlap

Within the shared universe, the frozen within-cohort candidate rule (FDR<0.05 and |log2FC|>=1) identified:

- GSE274832: **{candidate_counts['GSE274832']:,} genes**.
- GSE193136: **{candidate_counts['GSE193136']:,} genes**.
- GSE232306: **{candidate_counts['GSE232306']:,} genes**.
- Candidate in all three cohorts: **{all_three_candidates:,} genes**.
- Candidate in at least two cohorts with the same direction among candidate cohorts: **{two_plus_same:,} genes**.
- Same effect direction in all three cohorts regardless of significance: **{all_three_direction:,} genes**.

Thresholded overlaps are descriptive and are not external validation or biomarker confirmation.

## Meta-analysis and heterogeneity

Fixed-effect and DerSimonian-Laird random-effects estimates were computed for all shared genes. Random-effects FDR<0.05 occurred for **{random_fdr_n:,} genes**. I2 was <25% for **{low_i2_n:,} genes** and >=75% for **{high_i2_n:,} genes**. With only three cohorts, Q, tau2, and I2 are descriptive and imprecise.

The prespecified strict exploratory consensus rule returned **{strict_n:,} genes**. The heatmap selection was: **{selection_label}**. Even strict-set genes remain exploratory cross-cohort effects, not clinical biomarkers.

## Leave-one-cohort-out stability

LOCO rank correlations versus the full random-effects estimate ranged from **{loco_rho_min:.4f} to {loco_rho_max:.4f}**; direction agreement ranged from **{loco_dir_min:.2%} to {loco_dir_max:.2%}**.

{chr(10).join(loco_lines)}

## Interpretation boundary

- Direction/rank reproducibility is primary; inverse-variance synthesis is supportive.
- GSE232306 shows phenotype-aligned GC, duplication, and PC1 separation with unavailable individual age and incomplete batch/clinical covariates. Strong effects driven mainly by this cohort may reflect biological composition, technical differences, clinical differences, or a mixture.
- GSE274832 has n=3 per group. GSE193136 has low transcriptome mapping and an age imbalance handled by the frozen adjusted model.
- Gene-level heterogeneity can motivate M05 pathway/module convergence; it does not authorize pathway claims before M05.

## Permitted claims

- The project quantified cross-cohort direction, rank reproducibility, threshold overlap, and gene-level heterogeneity using frozen cohort-specific effects.
- Gene-level concordance and heterogeneity can be reported exactly as measured.
- Selected genes can be described as exploratory consensus effects only when the exact rule is stated.

## Prohibited claims

- No gene is a validated biomarker, diagnostic panel, causal driver, or clinical predictor based on M04.
- M04 is not an individual-level pooled analysis and is not an external clinical validation study.
- Meta-analysis FDR or two-cohort overlap must not be presented without heterogeneity and cohort-limit context.
"""
    (results / "M04_REPRO_HETEROGENEITY_SUMMARY.md").write_text(summary, encoding="utf-8")
    return key_counts


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

    tables: dict[str, pd.DataFrame] = {}
    audits: list[dict] = []
    for cohort in COHORTS:
        tables[cohort], audit = load_effect_table(project_root, cohort)
        audits.append(audit)

    common = tables[COHORTS[0]]
    for cohort in COHORTS[1:]:
        common = common.merge(tables[cohort], on="gene_id", how="inner", validate="one_to_one")
    if len(common) < 10000:
        raise ValueError(f"Common gene universe below gate: {len(common)} < 10000")
    common["gene_name"] = resolve_gene_names(common)
    common = common.drop(columns=[f"gene_name_{c}" for c in COHORTS])
    common = common.sort_values("gene_id").reset_index(drop=True)

    for cohort in COHORTS:
        common[f"abs_z_{cohort}"] = common[f"z_{cohort}"].abs()
        common[f"fdr_{cohort}"] = common[f"padj_{cohort}"] < 0.05
        common[f"candidate_{cohort}"] = (
            common[f"fdr_{cohort}"] & (common[f"beta_{cohort}"].abs() >= 1.0)
        )
    for audit in audits:
        audit["common_valid_genes"] = int(len(common))
    pd.DataFrame(audits).to_csv(results / "M04_INPUT_AUDIT.csv", index=False)

    pairwise, topk = build_pairwise_results(common)
    pairwise.to_csv(results / "M04_PAIRWISE_REPRODUCIBILITY.csv", index=False)
    topk.to_csv(results / "M04_TOPK_OVERLAP.csv", index=False)

    meta = add_meta_analysis(common)
    overlap = build_overlap_summary(meta)
    loco, gene_loco = build_loco(meta)
    meta.to_csv(results / "M04_GENE_LEVEL_META_ANALYSIS.csv", index=False)
    overlap.to_csv(results / "M04_CANDIDATE_OVERLAP_SUMMARY.csv", index=False)
    loco.to_csv(results / "M04_LOCO_STABILITY.csv", index=False)
    gene_loco.to_csv(results / "M04_GENE_LEVEL_LOCO_EFFECTS.csv", index=False)

    plot_pairwise_scatter(meta, pairwise, plots)
    plot_reproducibility_matrices(pairwise, plots)
    plot_candidate_overlap(overlap, plots)
    plot_heterogeneity(meta, plots)
    selection_label, _ = plot_effect_heatmap(meta, plots)
    plot_loco(loco, plots)
    key_counts = write_summary(results, meta, pairwise, topk, overlap, loco, selection_label)

    session = {
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "project_root": str(project_root),
        "run_dir": str(run_dir),
        "input_cohorts": COHORTS,
        "common_gene_universe": int(len(meta)),
        "sample_level_matrix_created": False,
    }
    (logs / "PYTHON_SESSION_INFO.json").write_text(
        json.dumps(session, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("M04 synthesis completed")
    print(f"Common genes: {len(meta):,}")
    for row in pairwise.itertuples(index=False):
        print(
            f"{row.cohort_a} vs {row.cohort_b}: rho={row.spearman_rho_log2fc:.4f}, "
            f"direction={row.direction_concordance_all_nonzero:.2%}"
        )
    print(key_counts.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
