from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B3_GSE193136"
REF = ROOT / "03_data/references/GENCODE_v50_GRCh38p14"
TX2GENE = REF / "gencode.v50.tx2gene.tsv.gz"
CONFIG = ROOT / "04_code/configs/GSE193136_SAMPLES.csv"
SALMON = RUN / "results/salmon"
OUT = RUN / "results/expression_qc"
OUT.mkdir(parents=True, exist_ok=True)

cfg = pd.read_csv(CONFIG, dtype={"srr": str, "gsm": str, "phenotype": str})
cfg["age_years"] = pd.to_numeric(cfg["age_years"])
if len(cfg) != 12 or cfg["srr"].nunique() != 12:
    raise RuntimeError("The frozen GSE193136 sample table must contain 12 unique samples")

tx = pd.read_csv(TX2GENE, sep="\t", dtype=str).drop_duplicates("transcript_id_version")
tpm_columns: dict[str, pd.Series] = {}
count_columns: dict[str, pd.Series] = {}
gene_meta: pd.DataFrame | None = None

for row in cfg.itertuples(index=False):
    quant_path = SALMON / row.srr / "quant.sf"
    if not quant_path.exists():
        raise FileNotFoundError(f"Missing Salmon quantification: {quant_path}")
    quant = pd.read_csv(quant_path, sep="\t")
    quant["transcript_id_version"] = quant["Name"].str.split("|", regex=False).str[0]
    merged = quant.merge(tx, on="transcript_id_version", how="inner")
    if len(merged) < int(len(quant) * 0.95):
        raise RuntimeError(f"Less than 95% of transcripts mapped to GENCODE genes for {row.srr}")
    grouped = (
        merged.groupby(["gene_id_version", "gene_id", "gene_name"], as_index=False)
        .agg(TPM=("TPM", "sum"), NumReads=("NumReads", "sum"))
        .set_index("gene_id_version")
    )
    if gene_meta is None:
        gene_meta = grouped[["gene_id", "gene_name"]]
    tpm_columns[row.srr] = grouped["TPM"]
    count_columns[row.srr] = grouped["NumReads"]

assert gene_meta is not None
sample_order = cfg["srr"].tolist()
tpm = pd.DataFrame(tpm_columns).fillna(0.0)[sample_order]
counts = pd.DataFrame(count_columns).fillna(0.0)[sample_order]
meta = gene_meta.reindex(tpm.index)

pd.concat([meta, tpm], axis=1).reset_index().to_csv(
    OUT / "GSE193136_GENE_TPM.tsv.gz", sep="\t", index=False, compression="gzip"
)
pd.concat([meta, counts], axis=1).reset_index().to_csv(
    OUT / "GSE193136_GENE_ESTIMATED_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip"
)

log_tpm = np.log2(tpm + 1.0)
variable = log_tpm.var(axis=1).sort_values(ascending=False).head(min(2000, len(log_tpm))).index
x = log_tpm.loc[variable].T.to_numpy(dtype=float, copy=True)
x -= x.mean(axis=0, keepdims=True)
u, singular, _ = np.linalg.svd(x, full_matrices=False)
scores = u * singular
variance_pct = singular**2 / np.sum(singular**2) * 100
pca = cfg[["srr", "gsm", "phenotype", "age_years"]].copy()
pca["PC1"] = scores[:, 0]
pca["PC2"] = scores[:, 1]
pca["PC1_variance_percent"] = variance_pct[0]
pca["PC2_variance_percent"] = variance_pct[1]
pca.to_csv(OUT / "GSE193136_PCA_COORDINATES.csv", index=False)

expressed = log_tpm.loc[(tpm >= 1).sum(axis=1) >= 3]
corr = expressed.corr(method="pearson")
corr.to_csv(OUT / "GSE193136_SAMPLE_CORRELATION.csv")
corr_qc = []
for sample in sample_order:
    vals = corr.loc[sample].drop(sample)
    corr_qc.append(
        {
            "srr": sample,
            "phenotype": cfg.set_index("srr").loc[sample, "phenotype"],
            "median_correlation_to_other_samples": float(vals.median()),
            "minimum_correlation_to_other_samples": float(vals.min()),
            "correlation_qc_verdict": "PASS" if vals.median() >= 0.80 else "REVIEW",
        }
    )
pd.DataFrame(corr_qc).to_csv(OUT / "GSE193136_SAMPLE_CORRELATION_QC.csv", index=False)

marker_groups = {
    "granulosa": ["FSHR", "LHCGR", "STAR", "CYP19A1", "CYP11A1", "FOXL2", "INHA", "INHBA", "AMH", "NR5A1", "HSD17B1"],
    "oocyte": ["ZP1", "ZP2", "ZP3", "GDF9", "BMP15", "FIGLA", "NOBOX"],
    "immune": ["PTPRC", "CD3D", "CD3E", "CD14", "CD68"],
    "retina": ["RHO", "GNAT1", "PDE6A", "PDE6B", "RCVRN", "CRX", "RPE65"],
}
marker_rows = []
for group, symbols in marker_groups.items():
    for symbol in symbols:
        hits = meta.index[meta["gene_name"] == symbol]
        for sample in sample_order:
            value = float(tpm.loc[hits, sample].sum()) if len(hits) else 0.0
            marker_rows.append(
                {
                    "marker_group": group,
                    "gene_name": symbol,
                    "srr": sample,
                    "phenotype": cfg.set_index("srr").loc[sample, "phenotype"],
                    "TPM": value,
                    "log2_TPM_plus_1": np.log2(value + 1.0),
                }
            )
markers = pd.DataFrame(marker_rows)
markers.to_csv(OUT / "GSE193136_MARKER_EXPRESSION.csv", index=False)

sample_marker = []
for sample in sample_order:
    sub = markers.loc[markers["srr"] == sample]
    gran = sub.loc[sub["marker_group"] == "granulosa", "TPM"]
    oocyte = sub.loc[sub["marker_group"] == "oocyte", "TPM"]
    immune = sub.loc[sub["marker_group"] == "immune", "TPM"]
    retina = sub.loc[sub["marker_group"] == "retina", "TPM"]
    sample_marker.append(
        {
            "srr": sample,
            "granulosa_markers_tpm_ge_1": int((gran >= 1).sum()),
            "granulosa_marker_median_tpm": float(gran.median()),
            "oocyte_markers_tpm_ge_1": int((oocyte >= 1).sum()),
            "immune_markers_tpm_ge_1": int((immune >= 1).sum()),
            "retina_markers_tpm_ge_1": int((retina >= 1).sum()),
            "granulosa_identity_verdict": "PASS" if int((gran >= 1).sum()) >= 4 else "REVIEW",
        }
    )
pd.DataFrame(sample_marker).to_csv(OUT / "GSE193136_MARKER_QC_BY_SAMPLE.csv", index=False)

pc_age = float(np.corrcoef(pca["PC1"], pca["age_years"])[0, 1])
summary = f"""# GSE193136 Expression QC Summary

- Quantified samples: {len(sample_order)}/12 (6 DOR, 6 NOR)
- Gene rows: {len(tpm):,}
- PCA input: top {len(variable):,} variable genes by log2(TPM+1) variance
- PC1 variance: {variance_pct[0]:.2f}%
- PC2 variance: {variance_pct[1]:.2f}%
- Pearson correlation between PC1 and age: {pc_age:.3f}
- Median age: DOR {cfg.loc[cfg.phenotype == 'DOR', 'age_years'].median():.1f}; NOR {cfg.loc[cfg.phenotype == 'NOR', 'age_years'].median():.1f}
- Interpretation boundary: PCA and marker panels are QC evidence, not differential-expression findings.
"""
(OUT / "GSE193136_EXPRESSION_QC_SUMMARY.md").write_text(summary, encoding="utf-8")
(RUN / "STEP05_EXPRESSION_QC.PASS.txt").write_text("PASS\tgene matrices and numerical QC generated\n", encoding="utf-8")
print(f"PASS: GSE193136 expression matrices and numerical QC generated for {len(sample_order)} samples")

