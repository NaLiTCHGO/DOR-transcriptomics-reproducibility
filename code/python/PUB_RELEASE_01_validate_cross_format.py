#!/usr/bin/env python3
"""Validate final Figure 2-S6 raster/PDF bundles against frozen candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageStat
from pypdf import PdfReader


SPECS = [
    ("Figure_2", "main", 2008, 2409, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_VISUAL_QC_ROUND6B_HEADER_S4/figures/Figure_2_Cohort_QC_and_Effects_v07_QC.png"),
    ("Figure_3", "main", 2008, 2409, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_BATCH_F3_S6_V04/figures/main/Figure_3_Gene_Reproducibility_and_Heterogeneity_v09_QC.png"),
    ("Figure_4", "main", 2008, 2600, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_BATCH_F3_S6_V04/figures/main/Figure_4_Pathway_Convergence_v05_QC.png"),
    ("Figure_5", "main", 2008, 2200, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_BATCH_F3_S6_V04/figures/main/Figure_5_LOCO_Robustness_v08_QC.png"),
    ("Figure_S1", "supplement", 2008, 2600, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_VISUAL_QC_ROUND5_MICROREPAIR_V02/figures/supplement/Figure_S1_GSE274832_QC_and_Stability_v08_QC.png"),
    ("Figure_S2", "supplement", 2008, 2600, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_VISUAL_QC_ROUND5_MICROREPAIR_V02/figures/supplement/Figure_S2_GSE193136_QC_and_Stability_v08_QC.png"),
    ("Figure_S3", "supplement", 2008, 2600, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_VISUAL_QC_ROUND5_MICROREPAIR_V02/figures/supplement/Figure_S3_GSE232306_QC_and_Stability_v08_QC.png"),
    ("Figure_S4", "supplement", 2008, 1150, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_VISUAL_QC_ROUND6B_HEADER_S4/figures/supplement/Figure_S4_Gene_Synthesis_Sensitivity_v05_QC.png"),
    ("Figure_S5", "supplement", 2008, 2300, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_BATCH_F3_S6_V04/figures/supplement/Figure_S5_Pathway_Context_v07_QC.png"),
    ("Figure_S6", "supplement", 2008, 1150, "05_analysis_steps/PUB_FIGURE_REBUILD/runs/20260824_BATCH_F3_S6_V04/figures/supplement/Figure_S6_LOCO_Context_v04_QC.png"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def normalized_rmse(left: Image.Image, right: Image.Image) -> float:
    left = left.convert("RGB")
    right = right.convert("RGB")
    if right.size != left.size:
        right = right.resize(left.size, Image.Resampling.LANCZOS)
    rms = ImageStat.Stat(ImageChops.difference(left, right)).rms
    return math.sqrt(sum(value * value for value in rms) / len(rms)) / 255.0


def average_hash_similarity(left: Image.Image, right: Image.Image, size: int = 32) -> float:
    """Content-level similarity robust to device-specific antialiasing."""
    def bits(image: Image.Image) -> list[bool]:
        small = image.convert("L").resize((size, size), Image.Resampling.LANCZOS)
        values = list(small.get_flattened_data())
        mean = sum(values) / len(values)
        return [value < mean for value in values]
    left_bits, right_bits = bits(left), bits(right)
    return 1.0 - sum(a != b for a, b in zip(left_bits, right_bits)) / len(left_bits)


def image_dpi(image: Image.Image) -> tuple[float, float]:
    dpi = image.info.get("dpi", (0.0, 0.0))
    if isinstance(dpi, (int, float)):
        return float(dpi), float(dpi)
    return float(dpi[0]), float(dpi[1])


def pdf_font_records(reader: PdfReader) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for page_no, page in enumerate(reader.pages, start=1):
        resources = page.get("/Resources")
        if resources is None:
            continue
        resources = resources.get_object()
        fonts = resources.get("/Font", {})
        fonts = fonts.get_object() if hasattr(fonts, "get_object") else fonts
        for key, ref in fonts.items():
            font = ref.get_object()
            descriptor = font.get("/FontDescriptor")
            if descriptor is None and font.get("/DescendantFonts"):
                descendant = font["/DescendantFonts"][0].get_object()
                descriptor = descendant.get("/FontDescriptor")
            descriptor = descriptor.get_object() if descriptor is not None else None
            embedded = bool(descriptor and any(name in descriptor for name in ("/FontFile", "/FontFile2", "/FontFile3")))
            records.append({
                "page": page_no,
                "resource": str(key),
                "base_font": str(font.get("/BaseFont", "UNKNOWN")),
                "subtype": str(font.get("/Subtype", "UNKNOWN")),
                "embedded": embedded,
            })
    return records


def make_contact_sheet(items: list[tuple[str, Path]], output: Path) -> None:
    thumb_w, thumb_h, gutter, label_h = 720, 800, 35, 34
    rows = math.ceil(len(items) / 2)
    canvas = Image.new("RGB", (2 * thumb_w + 3 * gutter, rows * (thumb_h + label_h) + (rows + 1) * gutter), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    for index, (label, path) in enumerate(items):
        image = Image.open(path).convert("RGB")
        image.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        col, row = index % 2, index // 2
        x = gutter + col * (thumb_w + gutter) + (thumb_w - image.width) // 2
        y = gutter + row * (thumb_h + label_h + gutter) + label_h
        draw.text((gutter + col * (thumb_w + gutter), y - label_h + 8), label, fill="black", font=font)
        canvas.paste(image, (x, y))
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, dpi=(150, 150), optimize=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    run_root = args.run_root.resolve()
    qc_root = run_root / "qc"
    pdf_render_root = qc_root / "pdf_rendered_300dpi"
    qc_root.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    hash_rows: list[dict[str, object]] = []
    preview_items: list[tuple[str, Path]] = []
    pdf_items: list[tuple[str, Path]] = []
    font_rows: list[dict[str, object]] = []

    for figure_id, category, width, height, accepted_rel in SPECS:
        folder = run_root / "figures" / category / figure_id
        preview = folder / f"{figure_id}_preview_300dpi.png"
        png600 = folder / f"{figure_id}_submission_600dpi.png"
        tiff600 = folder / f"{figure_id}_submission_600dpi.tiff"
        pdf = folder / f"{figure_id}_submission.pdf"
        rendered_pdf = pdf_render_root / f"{figure_id}_pdf_300dpi.png"
        accepted = project_root / accepted_rel
        required = [preview, png600, tiff600, pdf, rendered_pdf, accepted]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError("Missing QC input(s): " + "; ".join(missing))

        with Image.open(preview) as im_preview, Image.open(png600) as im_png, Image.open(tiff600) as im_tiff, Image.open(accepted) as im_accepted, Image.open(rendered_pdf) as im_pdf:
            preview_size_ok = im_preview.size == (width, height)
            png_size_ok = im_png.size == (width * 2, height * 2)
            tiff_size_ok = im_tiff.size == (width * 2, height * 2)
            preview_dpi = image_dpi(im_preview)
            png_dpi = image_dpi(im_png)
            tiff_dpi = image_dpi(im_tiff)
            preview_dpi_ok = all(299 <= value <= 301 for value in preview_dpi)
            png_dpi_ok = all(599 <= value <= 601 for value in png_dpi)
            tiff_dpi_ok = all(599 <= value <= 601 for value in tiff_dpi)
            tiff_compression = str(im_tiff.info.get("compression", ""))
            tiff_lzw_ok = "lzw" in tiff_compression.lower()
            accepted_rmse = normalized_rmse(im_preview, im_accepted)
            png_preview_rmse = normalized_rmse(im_preview, im_png)
            png_tiff_rmse = normalized_rmse(im_png, im_tiff)
            pdf_preview_rmse = normalized_rmse(im_preview, im_pdf)
            png_preview_hash = average_hash_similarity(im_preview, im_png)
            png_tiff_hash = average_hash_similarity(im_png, im_tiff)
            pdf_preview_hash = average_hash_similarity(im_preview, im_pdf)

        reader = PdfReader(str(pdf))
        page_count_ok = len(reader.pages) == 1
        page = reader.pages[0]
        width_mm = float(page.mediabox.width) / 72.0 * 25.4
        height_mm = float(page.mediabox.height) / 72.0 * 25.4
        expected_width_mm = round(width / 300.0 * 72.0) / 72.0 * 25.4
        expected_height_mm = round(height / 300.0 * 72.0) / 72.0 * 25.4
        pdf_size_ok = abs(width_mm - expected_width_mm) <= 0.05 and abs(height_mm - expected_height_mm) <= 0.05
        pdf_text = "".join((page.extract_text() or "") for page in reader.pages)
        pdf_text_ok = len(pdf_text.strip()) >= 20
        fonts = pdf_font_records(reader)
        for font in fonts:
            font_rows.append({"figure_id": figure_id, **font})
        fonts_embedded_ok = bool(fonts) and all(bool(font["embedded"]) for font in fonts)
        pdf_under_10mb = pdf.stat().st_size <= 10 * 1024 * 1024

        metrics_ok = (
            preview_size_ok and png_size_ok and tiff_size_ok and preview_dpi_ok and png_dpi_ok
            and tiff_dpi_ok and tiff_lzw_ok and accepted_rmse <= 0.035
            and png_preview_hash >= 0.985 and png_tiff_hash >= 0.995
            and pdf_preview_hash >= 0.900 and page_count_ok and pdf_size_ok
            and pdf_text_ok and fonts_embedded_ok and pdf_under_10mb
        )
        rows.append({
            "figure_id": figure_id,
            "preview_size_ok": preview_size_ok,
            "png600_size_ok": png_size_ok,
            "tiff600_size_ok": tiff_size_ok,
            "preview_dpi": f"{preview_dpi[0]:.3f}x{preview_dpi[1]:.3f}",
            "png_dpi": f"{png_dpi[0]:.3f}x{png_dpi[1]:.3f}",
            "tiff_dpi": f"{tiff_dpi[0]:.3f}x{tiff_dpi[1]:.3f}",
            "tiff_compression": tiff_compression,
            "accepted_preview_rmse": f"{accepted_rmse:.6f}",
            "png600_to_preview_rmse": f"{png_preview_rmse:.6f}",
            "png600_to_tiff600_rmse": f"{png_tiff_rmse:.6f}",
            "pdf_render_to_preview_rmse": f"{pdf_preview_rmse:.6f}",
            "png600_to_preview_hash_similarity": f"{png_preview_hash:.6f}",
            "png600_to_tiff600_hash_similarity": f"{png_tiff_hash:.6f}",
            "pdf_render_to_preview_hash_similarity": f"{pdf_preview_hash:.6f}",
            "pdf_width_mm": f"{width_mm:.3f}",
            "pdf_height_mm": f"{height_mm:.3f}",
            "pdf_one_page": page_count_ok,
            "pdf_text_extractable": pdf_text_ok,
            "pdf_fonts_embedded": fonts_embedded_ok,
            "pdf_under_10mb": pdf_under_10mb,
            "status": "PASS" if metrics_ok else "FAIL",
        })

        for artifact_type, path in (("preview_png", preview), ("submission_png", png600), ("submission_tiff", tiff600), ("submission_pdf", pdf)):
            hash_rows.append({
                "figure_id": figure_id,
                "artifact_type": artifact_type,
                "path": path.relative_to(run_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
        preview_items.append((figure_id, preview))
        pdf_items.append((figure_id, rendered_pdf))

    def write_csv(path: Path, records: list[dict[str, object]]) -> None:
        with path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
            writer.writeheader()
            writer.writerows(records)

    write_csv(qc_root / "CROSS_FORMAT_QC.csv", rows)
    write_csv(qc_root / "PDF_FONT_EMBEDDING_QC.csv", font_rows)
    write_csv(run_root / "manifests" / "FINAL_ARTIFACT_SHA256.csv", hash_rows)
    make_contact_sheet(preview_items, qc_root / "CONTACT_SHEET_FINAL_PREVIEWS.png")
    make_contact_sheet(pdf_items, qc_root / "CONTACT_SHEET_PDF_RENDERS.png")

    failures = [row["figure_id"] for row in rows if row["status"] != "PASS"]
    report = [
        "# Cross-format figure QC",
        "",
        "- Scope: Figure 2-Figure S6 final release artifacts.",
        "- Candidate identity: the 300-ppi release preview was compared with the user-approved candidate.",
        "- Independent render identity: 600-ppi PNG, LZW TIFF and rendered PDF were compared with the release preview using both pixel RMSE (reported) and content-level average-hash similarity (gated).",
        "- PDF checks: one page, exact physical dimensions, extractable text, embedded fonts and <10 MB.",
        "- Runtime-warning checks: stored in the renderer-specific warning reconciliation cards.",
        f"- Figures checked: {len(rows)}.",
        f"- Final status: {'PASS' if not failures else 'FAIL'}.",
    ]
    if failures:
        report.append("- Failed figures: " + ", ".join(failures) + ".")
    (qc_root / "CROSS_FORMAT_QC.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"{'PASS' if not failures else 'FAIL'}: cross-format QC for {len(rows)} figures")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
