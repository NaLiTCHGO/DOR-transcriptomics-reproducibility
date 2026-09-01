from __future__ import annotations

import csv
import gzip
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260813_M02_B1_GSE274832"
REF = ROOT / "03_data/references/GENCODE_v50_GRCh38p14"
GTF = REF / "gencode.v50.annotation.gtf.gz"
TX2GENE = REF / "gencode.v50.tx2gene.tsv.gz"
LOCKED = ROOT / "04_code/configs/GSE274832_SAMPLES.csv"
SALMON = RUN / "results/salmon"
OUT = RUN / "results/expression_qc"
OUT.mkdir(parents=True, exist_ok=True)


def attr_value(text: str, key: str) -> str | None:
    hit = re.search(rf'(?:^|;\s*){re.escape(key)} "([^"]+)"', text)
    return hit.group(1) if hit else None


def build_tx2gene() -> pd.DataFrame:
    if not TX2GENE.exists():
        temp = TX2GENE.with_suffix(TX2GENE.suffix + ".building")
        with gzip.open(GTF, "rt", encoding="utf-8") as src, gzip.open(
            temp, "wt", encoding="utf-8", newline=""
        ) as dst:
            writer = csv.writer(dst, delimiter="\t")
            writer.writerow(
                ["transcript_id_version", "transcript_id", "gene_id_version", "gene_id", "gene_name"]
            )
            for line in src:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "transcript":
                    continue
                attrs = fields[8]
                txv = attr_value(attrs, "transcript_id")
                genev = attr_value(attrs, "gene_id")
                if not txv or not genev:
                    continue
                writer.writerow(
                    [txv, txv.split(".")[0], genev, genev.split(".")[0], attr_value(attrs, "gene_name") or ""]
                )
        temp.replace(TX2GENE)
    return pd.read_csv(TX2GENE, sep="\t", dtype=str)


locked = pd.read_csv(LOCKED, dtype=str)
tx = build_tx2gene().drop_duplicates("transcript_id_version")

tpm_columns: dict[str, pd.Series] = {}
count_columns: dict[str, pd.Series] = {}
gene_meta: pd.DataFrame | None = None
for row in locked.itertuples(index=False):
    quant_path = SALMON / row.srr / "quant.sf"
    if not quant_path.exists():
        raise FileNotFoundError(f"Missing Salmon quantification: {quant_path}")
    quant = pd.read_csv(quant_path, sep="\t")
    # GENCODE FASTA headers contain additional pipe-delimited provenance
    # fields after the versioned ENST identifier; Salmon retains that header.
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
tpm = pd.DataFrame(tpm_columns).fillna(0.0)
counts = pd.DataFrame(count_columns).fillna(0.0)
meta = gene_meta.reindex(tpm.index)
sample_order = locked["srr"].tolist()
tpm = tpm[sample_order]
counts = counts[sample_order]

pd.concat([meta, tpm], axis=1).reset_index().to_csv(
    OUT / "GSE274832_GENE_TPM.tsv.gz", sep="\t", index=False, compression="gzip"
)
pd.concat([meta, counts], axis=1).reset_index().to_csv(
    OUT / "GSE274832_GENE_ESTIMATED_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip"
)

log_tpm = np.log2(tpm + 1.0)
variable = log_tpm.var(axis=1).sort_values(ascending=False).head(min(2000, len(log_tpm))).index
x = log_tpm.loc[variable].T.to_numpy(dtype=float, copy=True)
x -= x.mean(axis=0, keepdims=True)
u, singular, _ = np.linalg.svd(x, full_matrices=False)
scores = u * singular
variance = singular**2
variance_pct = variance / variance.sum() * 100
pca = locked[["srr", "gsm", "phenotype"]].copy()
pca["PC1"] = scores[:, 0]
pca["PC2"] = scores[:, 1]
pca["PC1_variance_percent"] = variance_pct[0]
pca["PC2_variance_percent"] = variance_pct[1]
pca.to_csv(OUT / "GSE274832_PCA_COORDINATES.csv", index=False)

expressed = log_tpm.loc[(tpm >= 1).sum(axis=1) >= 2]
corr = expressed.corr(method="pearson")
corr.to_csv(OUT / "GSE274832_SAMPLE_CORRELATION.csv")
corr_qc = []
for sample in sample_order:
    vals = corr.loc[sample].drop(sample)
    corr_qc.append(
        {
            "srr": sample,
            "phenotype": locked.set_index("srr").loc[sample, "phenotype"],
            "median_correlation_to_other_samples": float(vals.median()),
            "minimum_correlation_to_other_samples": float(vals.min()),
            "correlation_qc_verdict": "PASS" if vals.median() >= 0.80 else "REVIEW",
        }
    )
pd.DataFrame(corr_qc).to_csv(OUT / "GSE274832_SAMPLE_CORRELATION_QC.csv", index=False)

marker_groups = {
    "granulosa": ["FSHR", "LHCGR", "STAR", "CYP19A1", "CYP11A1", "FOXL2", "INHA", "INHBA", "AMH", "NR5A1", "HSD17B1"],
    "retina": ["RHO", "GNAT1", "PDE6A", "PDE6B", "RCVRN", "CRX", "OPN1MW", "OPN1LW", "VSX2", "RPE65"],
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
                    "phenotype": locked.set_index("srr").loc[sample, "phenotype"],
                    "TPM": value,
                    "log2_TPM_plus_1": np.log2(value + 1.0),
                }
            )
pd.DataFrame(marker_rows).to_csv(OUT / "GSE274832_MARKER_EXPRESSION.csv", index=False)

max_pairs = []
for i, a in enumerate(sample_order):
    for b in sample_order[i + 1 :]:
        max_pairs.append((float(corr.loc[a, b]), a, b))
max_pairs.sort(reverse=True)

with open(OUT / "GSE274832_EXPRESSION_QC_SUMMARY.md", "w", encoding="utf-8") as handle:
    handle.write("# GSE274832 Expression QC Summary\n\n")
    handle.write(f"- Quantified samples: {len(sample_order)}/6\n")
    handle.write(f"- Gene rows: {len(tpm):,}\n")
    handle.write(f"- PCA input: top {len(variable):,} variable genes by log2(TPM+1) variance\n")
    handle.write(f"- PC1 variance: {variance_pct[0]:.2f}%\n")
    handle.write(f"- PC2 variance: {variance_pct[1]:.2f}%\n")
    handle.write(f"- Highest off-diagonal sample correlation: {max_pairs[0][0]:.4f} ({max_pairs[0][1]} vs {max_pairs[0][2]})\n")
    handle.write("- Interpretation boundary: PCA separation is descriptive only for n=3 vs n=3 and is not a DEG result.\n")

print(f"PASS: expression matrices and QC summaries generated for {len(sample_order)} samples")
