#!/usr/bin/env python3
"""Generate the locked DOR Figure 1 submission artwork from one layout spec.

The script is the project-specific entry point required by the Figure 1 Suite.
It uses only the frozen helpers and inputs inside the submission package.
No GUI post-edit is required or permitted in the reproducibility chain.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image


OUTPUT_STEM = "DOR_Figure1_final_v01"
TARGET_WIDTH_MM = 180.0
FINAL_DPI = 600


def run_checked(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        check=False,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    if result.returncode != 0:
        detail = (result.stdout or "") + ("\n" + result.stderr if result.stderr else "")
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n{detail}")
    return result


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def find_font(candidates: list[Path]) -> Path:
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("Required publication font was not found: " + ", ".join(map(str, candidates)))


def freeze_physical_svg_size(svg_path: Path, canvas_width: int, canvas_height: int) -> float:
    height_mm = TARGET_WIDTH_MM * canvas_height / canvas_width
    svg_text = svg_path.read_text(encoding="utf-8")
    old_width = f'width="{canvas_width}"'
    old_height = f'height="{canvas_height}"'
    if old_width not in svg_text or old_height not in svg_text:
        raise ValueError("Rendered SVG does not contain the expected canvas dimensions")
    svg_text = svg_text.replace(old_width, f'width="{TARGET_WIDTH_MM:.4f}mm"', 1)
    svg_text = svg_text.replace(old_height, f'height="{height_mm:.4f}mm"', 1)
    svg_path.write_text(svg_text, encoding="utf-8")
    return height_mm


def validate_frozen_provenance(package_root: Path, spec: dict) -> dict:
    core_path = package_root / "source" / "frozen_inputs" / "CORE_RNASEQ_COHORT_SUMMARY.csv"
    etab_path = package_root / "source" / "frozen_inputs" / "E-MTAB-391_PROVENANCE_REPORT.md"

    with core_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = {row["series"]: row for row in csv.DictReader(handle)}

    expected = {
        "GSE274832": {"context": "147", "omics": "6", "independent": "6", "dor": "3", "nor": "3"},
        "GSE193136": {"context": "24", "omics": "12", "independent": "12", "dor": "6", "nor": "6"},
        "GSE232306": {"context": "60", "omics": "12", "independent": "12", "dor": "6", "nor": "6"},
    }
    failures: list[str] = []
    for series, values in expected.items():
        row = rows.get(series)
        if row is None:
            failures.append(f"Missing locked cohort row: {series}")
            continue
        if values["context"] not in row["N_recruited_context"]:
            failures.append(f"{series} recruited/context count mismatch")
        for key, column in (
            ("omics", "N_omics"),
            ("independent", "N_independent"),
            ("dor", "DOR_omics"),
            ("nor", "NOR_omics"),
        ):
            if row[column] != values[key]:
                failures.append(f"{series} {column} mismatch: {row[column]} != {values[key]}")

    etab_text = etab_path.read_text(encoding="utf-8")
    for required in ("28", "26", "15 DOR + 13 NOR"):
        if required not in etab_text:
            failures.append(f"E-MTAB-391 provenance text is missing: {required}")

    spec_text = json.dumps(spec, ensure_ascii=False)
    for required in (
        "GSE274832",
        "GSE193136",
        "GSE232306",
        "147 clinical context",
        "30 independent transcriptomes",
        "15 DOR + 15 NOR",
        "28 cycle-samples / 26 subjects",
        "not pooled into core synthesis",
    ):
        if required not in spec_text:
            failures.append(f"Locked layout is missing required provenance wording: {required}")

    report = {
        "hard_pass": not failures,
        "failures": failures,
        "core_cohorts_checked": sorted(expected),
        "etab391_cycle_samples": 28,
        "etab391_independent_subjects": 26,
        "stochastic": False,
    }
    if failures:
        raise ValueError("Frozen provenance validation failed: " + "; ".join(failures))
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--package-root",
        type=Path,
        help="Submission package root; defaults to the parent of the code directory.",
    )
    parser.add_argument(
        "--inkscape",
        type=Path,
        default=Path(r"C:\Program Files\Inkscape\bin\inkscape.exe"),
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    package_root = (args.package_root or script_path.parents[1]).resolve()
    inkscape = args.inkscape.resolve()

    spec_path = package_root / "config" / "Figure1_layout_spec.json"
    legend_source = package_root / "config" / "Figure1_legend.txt"
    vendor = package_root / "code" / "vendor"
    renderer = vendor / "render_layout_spec.py"
    layout_validator = vendor / "validate_layout_spec.py"
    text_validator = vendor / "validate_rendered_output.py"
    grammar_validator = vendor / "validate_visual_grammar.py"
    outputs = package_root / "outputs"
    qc = package_root / "qc"
    environment = package_root / "environment"
    for folder in (outputs, qc, environment):
        folder.mkdir(parents=True, exist_ok=True)

    required = [
        spec_path,
        legend_source,
        renderer,
        layout_validator,
        text_validator,
        grammar_validator,
        inkscape,
        package_root / "source" / "Figure1_source_manifest.csv",
        package_root / "source" / "frozen_inputs" / "CORE_RNASEQ_COHORT_SUMMARY.csv",
        package_root / "source" / "frozen_inputs" / "E-MTAB-391_PROVENANCE_REPORT.md",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing final Figure 1 input/helper:\n" + "\n".join(missing))

    layout = json.loads(spec_path.read_text(encoding="utf-8"))
    canvas_width = int(layout["canvas"]["width_px"])
    canvas_height = int(layout["canvas"]["height_px"])
    if float(layout["publication"]["target_width_mm"]) != TARGET_WIDTH_MM:
        raise ValueError("Layout target width is not the locked 180 mm")
    if int(layout["publication"]["final_dpi"]) != FINAL_DPI:
        raise ValueError("Layout final DPI is not the locked 600 dpi")

    provenance_report = validate_frozen_provenance(package_root, layout)
    (qc / "Figure1_final_provenance_validation.json").write_text(
        json.dumps(provenance_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    layout_result = run_checked(
        [sys.executable, str(layout_validator), str(spec_path), "--min-pt", "7.5"]
    )
    (qc / "Figure1_final_layout_preflight.json").write_text(layout_result.stdout, encoding="utf-8")

    svg_path = outputs / f"{OUTPUT_STEM}.svg"
    preview_path = outputs / f"{OUTPUT_STEM}_preview.png"
    png_path = outputs / f"{OUTPUT_STEM}_600dpi.png"
    tiff_path = outputs / f"{OUTPUT_STEM}_600dpi.tiff"
    pdf_path = outputs / f"{OUTPUT_STEM}.pdf"
    legend_path = outputs / "Figure1_legend_final.txt"

    run_checked([sys.executable, str(renderer), str(spec_path), "--svg", str(svg_path)])
    physical_height_mm = freeze_physical_svg_size(svg_path, canvas_width, canvas_height)

    run_checked(
        [
            str(inkscape),
            str(svg_path),
            "--export-area-page",
            "--export-type=png",
            f"--export-filename={preview_path}",
            f"--export-width={canvas_width}",
            f"--export-height={canvas_height}",
        ]
    )

    final_width_px = round(TARGET_WIDTH_MM / 25.4 * FINAL_DPI)
    final_height_px = round(final_width_px * canvas_height / canvas_width)
    run_checked(
        [
            str(inkscape),
            str(svg_path),
            "--export-area-page",
            "--export-type=png",
            f"--export-filename={png_path}",
            f"--export-width={final_width_px}",
            f"--export-height={final_height_px}",
        ]
    )
    run_checked(
        [
            str(inkscape),
            str(svg_path),
            "--export-area-page",
            "--export-type=pdf",
            f"--export-filename={pdf_path}",
        ]
    )

    with Image.open(png_path) as image:
        rgb = image.convert("RGB")
        if rgb.size != (final_width_px, final_height_px):
            raise ValueError(f"Unexpected final PNG size: {rgb.size}")
        temp_png = png_path.with_suffix(".tmp.png")
        rgb.save(temp_png, format="PNG", dpi=(FINAL_DPI, FINAL_DPI), compress_level=6)
        temp_png.replace(png_path)
        rgb.save(
            tiff_path,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(FINAL_DPI, FINAL_DPI),
        )

    shutil.copy2(legend_source, legend_path)

    regular_font = find_font(
        [
            Path(r"C:\Windows\Fonts\arial.ttf"),
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        ]
    )
    bold_font = find_font(
        [
            Path(r"C:\Windows\Fonts\arialbd.ttf"),
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ]
    )
    run_checked(
        [
            sys.executable,
            str(text_validator),
            str(spec_path),
            "--regular-font",
            str(regular_font),
            "--bold-font",
            str(bold_font),
            "--report",
            str(qc / "Figure1_final_rendered_text_fit.json"),
        ]
    )
    grammar_result = run_checked(
        [sys.executable, str(grammar_validator), str(spec_path), "--profile", "minimal_journal_workflow"]
    )
    (qc / "Figure1_final_visual_grammar.json").write_text(grammar_result.stdout, encoding="utf-8")

    with Image.open(preview_path) as preview_image:
        preview_metrics = {"size_px": list(preview_image.size), "mode": preview_image.mode}
    with Image.open(png_path) as final_image:
        final_metrics = {
            "size_px": list(final_image.size),
            "mode": final_image.mode,
            "dpi": [round(value, 3) for value in final_image.info.get("dpi", (0, 0))],
        }
    with Image.open(tiff_path) as tiff_image:
        tiff_metrics = {
            "size_px": list(tiff_image.size),
            "mode": tiff_image.mode,
            "dpi": [float(value) for value in tiff_image.info.get("dpi", (0, 0))],
            "compression": tiff_image.info.get("compression", "unknown"),
        }

    inkscape_version = run_checked([str(inkscape), "--version"]).stdout.strip()
    environment_path = environment / "Figure1_environment.txt"
    environment_path.write_text(
        "\n".join(
            [
                f"generated_at={datetime.now().astimezone().isoformat(timespec='seconds')}",
                f"python={platform.python_version()}",
                f"platform={platform.platform()}",
                f"pillow={package_version('Pillow')}",
                f"inkscape={inkscape_version}",
                f"regular_font={regular_font}",
                f"bold_font={bold_font}",
                "stochastic=false",
                "manual_post_edit_required=false",
                "render_chain=Figure1_layout_spec.json -> frozen SVG renderer -> physical-size SVG -> Inkscape PDF/PNG -> Pillow TIFF/DPI normalization",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    output_report = {
        "hard_pass": True,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "target_width_mm": TARGET_WIDTH_MM,
        "target_height_mm": round(physical_height_mm, 4),
        "final_dpi": FINAL_DPI,
        "preview": preview_metrics,
        "png_600dpi": final_metrics,
        "tiff_600dpi": tiff_metrics,
        "outputs": {
            "svg": str(svg_path.relative_to(package_root)),
            "pdf": str(pdf_path.relative_to(package_root)),
            "png_600dpi": str(png_path.relative_to(package_root)),
            "tiff_600dpi": str(tiff_path.relative_to(package_root)),
            "preview_png": str(preview_path.relative_to(package_root)),
            "legend": str(legend_path.relative_to(package_root)),
        },
    }
    (qc / "Figure1_final_output_validation.json").write_text(
        json.dumps(output_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(output_report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
