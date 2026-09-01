#!/usr/bin/env python3
"""Validate the final DOR Figure 1 package against the user-approved v03 preview."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from PIL import Image, ImageChops
from pypdf import PdfReader


def resolve(value):
    return value.get_object() if hasattr(value, "get_object") else value


def pdf_font_audit(reader: PdfReader) -> list[dict]:
    page = reader.pages[0]
    resources = resolve(page.get("/Resources", {}))
    fonts = resolve(resources.get("/Font", {})) if resources else {}
    records: list[dict] = []
    for resource_name, reference in fonts.items():
        font = resolve(reference)
        descriptors = []
        direct_descriptor = font.get("/FontDescriptor")
        if direct_descriptor is not None:
            descriptors.append(resolve(direct_descriptor))
        for descendant_reference in font.get("/DescendantFonts", []):
            descendant = resolve(descendant_reference)
            descriptor = descendant.get("/FontDescriptor")
            if descriptor is not None:
                descriptors.append(resolve(descriptor))
        subtype = str(font.get("/Subtype", ""))
        embedded = subtype == "/Type3" or any(
            any(key in descriptor for key in ("/FontFile", "/FontFile2", "/FontFile3"))
            for descriptor in descriptors
        )
        records.append(
            {
                "resource": str(resource_name),
                "base_font": str(font.get("/BaseFont", "")),
                "subtype": subtype,
                "embedded": embedded,
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--approved-preview", required=True, type=Path)
    args = parser.parse_args()

    root = args.package_root.resolve()
    approved = args.approved_preview.resolve()
    outputs = root / "outputs"
    qc = root / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    svg_path = outputs / "DOR_Figure1_final_v01.svg"
    pdf_path = outputs / "DOR_Figure1_final_v01.pdf"
    preview_path = outputs / "DOR_Figure1_final_v01_preview.png"
    png_path = outputs / "DOR_Figure1_final_v01_600dpi.png"
    tiff_path = outputs / "DOR_Figure1_final_v01_600dpi.tiff"
    required = [svg_path, pdf_path, preview_path, png_path, tiff_path, approved]
    failures = [f"Missing required file: {path}" for path in required if not path.is_file()]

    svg_text = svg_path.read_text(encoding="utf-8") if svg_path.is_file() else ""
    width_match = re.search(r'width="([0-9.]+)mm"', svg_text)
    height_match = re.search(r'height="([0-9.]+)mm"', svg_text)
    svg_width_mm = float(width_match.group(1)) if width_match else 0.0
    svg_height_mm = float(height_match.group(1)) if height_match else 0.0
    if abs(svg_width_mm - 180.0) > 0.01:
        failures.append(f"SVG width is not 180 mm: {svg_width_mm}")
    if abs(svg_height_mm - 180.9474) > 0.02:
        failures.append(f"SVG height mismatch: {svg_height_mm}")

    preview_pixel_identical = False
    preview_render_equivalent = False
    comparison_metrics = {}
    if preview_path.is_file() and approved.is_file():
        with Image.open(preview_path) as generated_image, Image.open(approved) as approved_image:
            generated_rgb = generated_image.convert("RGB")
            approved_rgb = approved_image.convert("RGB")
            same_size = generated_rgb.size == approved_rgb.size
            if same_size:
                difference = ImageChops.difference(generated_rgb, approved_rgb)
                pixels = list(difference.get_flattened_data())
                changed_pixels = sum(pixel != (0, 0, 0) for pixel in pixels)
                channel_abs_sum = sum(sum(pixel) for pixel in pixels)
                max_channel_difference = max(max(pixel) for pixel in pixels)
                total_pixels = generated_rgb.width * generated_rgb.height
                mean_absolute_channel_difference = channel_abs_sum / (total_pixels * 3)
                changed_pixel_percent = changed_pixels / total_pixels * 100
                preview_pixel_identical = changed_pixels == 0
                preview_render_equivalent = (
                    mean_absolute_channel_difference <= 0.01
                    and changed_pixel_percent <= 0.5
                    and max_channel_difference <= 20
                )
                comparison_metrics = {
                    "same_size": True,
                    "mean_absolute_channel_difference": mean_absolute_channel_difference,
                    "changed_pixels": changed_pixels,
                    "changed_pixel_percent": changed_pixel_percent,
                    "max_channel_difference": max_channel_difference,
                    "thresholds": {
                        "mean_absolute_channel_difference_max": 0.01,
                        "changed_pixel_percent_max": 0.5,
                        "max_channel_difference_max": 20,
                    },
                }
            else:
                comparison_metrics = {
                    "same_size": False,
                    "generated_size": list(generated_rgb.size),
                    "approved_size": list(approved_rgb.size),
                }
    if not preview_render_equivalent:
        failures.append("Final preview differs materially from the user-approved v03 candidate")

    raster_metrics = {}
    if png_path.is_file():
        with Image.open(png_path) as image:
            png_dpi = image.info.get("dpi", (0, 0))
            raster_metrics["png"] = {
                "size_px": list(image.size),
                "mode": image.mode,
                "dpi": [float(value) for value in png_dpi],
            }
            if image.size != (4252, 4274):
                failures.append(f"Final PNG dimensions mismatch: {image.size}")
            if min(png_dpi) < 599.0:
                failures.append(f"Final PNG DPI is below 600: {png_dpi}")
    if tiff_path.is_file():
        with Image.open(tiff_path) as image:
            tiff_dpi = image.info.get("dpi", (0, 0))
            compression = image.info.get("compression", "unknown")
            raster_metrics["tiff"] = {
                "size_px": list(image.size),
                "mode": image.mode,
                "dpi": [float(value) for value in tiff_dpi],
                "compression": compression,
            }
            if image.size != (4252, 4274):
                failures.append(f"Final TIFF dimensions mismatch: {image.size}")
            if min(tiff_dpi) < 599.0:
                failures.append(f"Final TIFF DPI is below 600: {tiff_dpi}")
            if compression != "tiff_lzw":
                failures.append(f"Final TIFF is not LZW-compressed: {compression}")

    pdf_metrics = {}
    font_records = []
    if pdf_path.is_file():
        reader = PdfReader(str(pdf_path))
        if len(reader.pages) != 1:
            failures.append(f"Final PDF page count is not 1: {len(reader.pages)}")
        page = reader.pages[0]
        width_pt = float(page.mediabox.width)
        height_pt = float(page.mediabox.height)
        width_mm = width_pt / 72 * 25.4
        height_mm = height_pt / 72 * 25.4
        pdf_metrics = {
            "pages": len(reader.pages),
            "width_pt": width_pt,
            "height_pt": height_pt,
            "width_mm": width_mm,
            "height_mm": height_mm,
        }
        if abs(width_mm - 180.0) > 0.05:
            failures.append(f"Final PDF width is not 180 mm: {width_mm:.4f}")
        if abs(height_mm - 180.9474) > 0.05:
            failures.append(f"Final PDF height mismatch: {height_mm:.4f}")
        font_records = pdf_font_audit(reader)
        if not font_records:
            failures.append("Final PDF contains no auditable font resources")
        elif any(not record["embedded"] for record in font_records):
            failures.append("Final PDF contains an unembedded font resource")

    report = {
        "hard_pass": not failures,
        "failures": failures,
        "approved_preview": str(approved),
        "preview_pixel_identical": preview_pixel_identical,
        "preview_render_equivalent": preview_render_equivalent,
        "approved_preview_comparison": comparison_metrics,
        "svg": {"width_mm": svg_width_mm, "height_mm": svg_height_mm},
        "pdf": pdf_metrics,
        "pdf_fonts": font_records,
        "rasters": raster_metrics,
        "manual_post_edit_required": False,
    }
    report_path = qc / "Figure1_final_release_validation.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
