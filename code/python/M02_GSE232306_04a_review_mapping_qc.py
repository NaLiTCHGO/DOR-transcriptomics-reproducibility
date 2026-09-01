from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B4_GSE232306"
SALMON = RUN / "results/salmon"
TX2GENE = ROOT / "03_data/references/GENCODE_v50_GRCh38p14/gencode.v50.tx2gene.tsv.gz"
OUT = RUN / "results/salmon_mapping_review"
OUT.mkdir(parents=True, exist_ok=True)

tx = pd.read_csv(TX2GENE, sep="\t", dtype=str).drop_duplicates("transcript_id_version")
markers = ["FSHR", "LHCGR", "STAR", "CYP19A1", "CYP11A1", "FOXL2", "INHA", "INHBA", "AMH", "NR5A1", "HSD17B1"]
rows = []
marker_rows = []

for quant_path in sorted(SALMON.glob("SRR*/quant.sf")):
    srr = quant_path.parent.name
    meta_path = quant_path.parent / "aux_info/meta_info.json"
    if not meta_path.exists():
        continue
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    quant = pd.read_csv(quant_path, sep="\t")
    quant["transcript_id_version"] = quant["Name"].str.split("|", regex=False).str[0]
    merged = quant.merge(tx, on="transcript_id_version", how="inner")
    gene = merged.groupby("gene_name", as_index=False).agg(TPM=("TPM", "sum"), NumReads=("NumReads", "sum"))
    detected = int((gene["TPM"] >= 1).sum())
    marker_table = gene.set_index("gene_name").reindex(markers).fillna(0.0)
    marker_detected = int((marker_table["TPM"] >= 1).sum())
    rows.append(
        {
            "srr": srr,
            "num_processed": int(meta["num_processed"]),
            "num_mapped": int(meta["num_mapped"]),
            "percent_mapped": float(meta["percent_mapped"]),
            "library_types": ";".join(meta.get("library_types", [])),
            "genes_tpm_ge_1": detected,
            "granulosa_markers_tpm_ge_1": marker_detected,
            "total_gene_numreads": float(gene["NumReads"].sum()),
        }
    )
    for marker, mrow in marker_table.iterrows():
        marker_rows.append({"srr": srr, "gene_name": marker, "TPM": float(mrow["TPM"])})

review = pd.DataFrame(rows).sort_values("srr")
if len(review) != 12:
    raise RuntimeError(f"All 12 Salmon samples are required for mapping review; found {len(review)}")
review.to_csv(OUT / "GSE232306_SALMON_MAPPING_REVIEW.csv", index=False)
pd.DataFrame(marker_rows).to_csv(OUT / "GSE232306_SALMON_MAPPING_GRANULOSA_MARKERS.csv", index=False)

mapping_consistent = review["percent_mapped"].max() - review["percent_mapped"].min() <= 10
identity_supported = (review["granulosa_markers_tpm_ge_1"] >= 4).all()
coverage_supported = (review["genes_tpm_ge_1"] >= 8000).all()
verdict = "PASS" if mapping_consistent and identity_supported and coverage_supported else "HOLD_FOR_TECHNICAL_REVIEW"
summary = f"""# GSE232306 Salmon Mapping and Coverage Review

- Completed samples: {len(review)}/12
- Salmon transcriptome mapping range: {review.percent_mapped.min():.2f}%–{review.percent_mapped.max():.2f}%
- Detected genes (TPM >= 1) range: {review.genes_tpm_ge_1.min()}–{review.genes_tpm_ge_1.max()}
- Prespecified granulosa markers detected (TPM >= 1) range: {review.granulosa_markers_tpm_ge_1.min()}–{review.granulosa_markers_tpm_ge_1.max()} of {len(markers)}
- Inferred library types: {', '.join(sorted(review.library_types.unique()))}
- Joint mapping/coverage verdict: `{verdict}`

This review does not by itself declare M02 PASS. The final M02 decision also requires cohort-wide sample correlation/PCA and granulosa-identity QC.
"""
(OUT / "GSE232306_SALMON_MAPPING_REVIEW.md").write_text(summary, encoding="utf-8")
(RUN / "STEP04A_MAPPING_REVIEW.PASS.txt").write_text(f"{verdict}\t{len(review)}/12 samples\n", encoding="utf-8")
print(summary)
