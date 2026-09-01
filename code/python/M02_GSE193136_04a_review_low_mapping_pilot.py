from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("DOR_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
RUN = ROOT / "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B3_GSE193136"
SALMON = RUN / "results/salmon"
TX2GENE = ROOT / "03_data/references/GENCODE_v50_GRCh38p14/gencode.v50.tx2gene.tsv.gz"
OUT = RUN / "results/salmon_pilot_review"
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
if len(review) < 3:
    raise RuntimeError(f"At least three complete Salmon pilot samples are required; found {len(review)}")
review.to_csv(OUT / "GSE193136_SALMON_LOW_MAPPING_PILOT.csv", index=False)
pd.DataFrame(marker_rows).to_csv(OUT / "GSE193136_SALMON_PILOT_GRANULOSA_MARKERS.csv", index=False)

mapping_consistent = review["percent_mapped"].max() - review["percent_mapped"].min() <= 10
identity_supported = (review["granulosa_markers_tpm_ge_1"] >= 4).all()
coverage_supported = (review["genes_tpm_ge_1"] >= 8000).all()
verdict = "PROVISIONAL_PASS_WITH_LIMITATION" if mapping_consistent and identity_supported and coverage_supported else "HOLD_FOR_TECHNICAL_REVIEW"
summary = f"""# GSE193136 Low Transcriptome-Mapping Pilot Review

- Completed pilot samples: {len(review)}
- Salmon transcriptome mapping range: {review.percent_mapped.min():.2f}%–{review.percent_mapped.max():.2f}%
- Detected genes (TPM >= 1) range: {review.genes_tpm_ge_1.min()}–{review.genes_tpm_ge_1.max()}
- Prespecified granulosa markers detected (TPM >= 1) range: {review.granulosa_markers_tpm_ge_1.min()}–{review.granulosa_markers_tpm_ge_1.max()} of {len(markers)}
- Inferred library types: {', '.join(sorted(review.library_types.unique()))}
- Pilot verdict: `{verdict}`

This pilot does not declare M02 PASS. It only determines whether completing all 12 Salmon quantifications is technically justified. The final M02 decision still requires cohort-wide mapping consistency, correlation/PCA and granulosa-identity QC.
"""
(OUT / "GSE193136_SALMON_LOW_MAPPING_PILOT_REVIEW.md").write_text(summary, encoding="utf-8")
(RUN / "STEP04A_LOW_MAPPING_PILOT.txt").write_text(f"{verdict}\t{len(review)} pilot samples\n", encoding="utf-8")
print(summary)

