#!/usr/bin/env python3
"""Run structural and real-regeneration QC for the locked DOR Figure 1 package."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image, ImageChops


IMAGE_OUTPUTS = (
    "outputs/DOR_Figure1_final_v01_preview.png",
    "outputs/DOR_Figure1_final_v01_600dpi.png",
    "outputs/DOR_Figure1_final_v01_600dpi.tiff",
)
TEXT_OUTPUTS = (
    "outputs/DOR_Figure1_final_v01.svg",
    "outputs/Figure1_legend_final.txt",
)


def run(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )


def image_identical(left: Path, right: Path) -> tuple[bool, dict[str, object]]:
    with Image.open(left) as a0, Image.open(right) as b0:
        a = a0.convert("RGBA")
        b = b0.convert("RGBA")
        same_size = a.size == b.size
        identical = same_size and ImageChops.difference(a, b).getbbox() is None
        return identical, {
            "locked_size_px": list(a.size),
            "regenerated_size_px": list(b.size),
            "pixel_identical": identical,
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", required=True)
    args = parser.parse_args()

    package = Path(args.package_root).resolve()
    qc_dir = package / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    structural_cmd = [
        sys.executable,
        str(package / "code/vendor/validate_reproducibility_package.py"),
        str(package),
    ]
    structural = run(structural_cmd, package)
    try:
        structural_report = json.loads(structural.stdout)
    except json.JSONDecodeError:
        structural_report = {
            "hard_pass": False,
            "failures": ["Structural validator did not return valid JSON"],
            "stdout": structural.stdout,
            "stderr": structural.stderr,
        }
    structural_report["validator_exit_code"] = structural.returncode
    (qc_dir / "Figure1_reproducibility_structural_qc.json").write_text(
        json.dumps(structural_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    replay_report: dict[str, object] = {
        "package": str(package),
        "test": "independent temporary-copy one-command regeneration",
        "declared_command": "python code/generate_Figure1.py --package-root .",
        "structural_hard_pass": bool(structural_report.get("hard_pass", False)),
        "comparisons": {},
        "failures": [],
    }
    failures: list[str] = replay_report["failures"]  # type: ignore[assignment]

    if not replay_report["structural_hard_pass"]:
        failures.append("Structural reproducibility gate failed")
    else:
        with tempfile.TemporaryDirectory(prefix="dor_figure1_replay_") as tmp:
            replay_package = Path(tmp) / "DOR_FIGURE1_REPLAY"
            shutil.copytree(package, replay_package)
            generator = run(
                [sys.executable, "code/generate_Figure1.py", "--package-root", "."],
                replay_package,
            )
            replay_report["generator_exit_code"] = generator.returncode
            replay_report["generator_stdout_tail"] = generator.stdout[-4000:]
            replay_report["generator_stderr_tail"] = generator.stderr[-4000:]
            if generator.returncode != 0:
                failures.append("Generator failed in the independent temporary package")
            else:
                comparisons: dict[str, object] = replay_report["comparisons"]  # type: ignore[assignment]
                for rel in IMAGE_OUTPUTS:
                    identical, metrics = image_identical(package / rel, replay_package / rel)
                    comparisons[rel] = metrics
                    if not identical:
                        failures.append(f"Regenerated image is not pixel-identical: {rel}")
                for rel in TEXT_OUTPUTS:
                    identical = (package / rel).read_bytes() == (replay_package / rel).read_bytes()
                    comparisons[rel] = {"byte_identical": identical}
                    if not identical:
                        failures.append(f"Regenerated text/vector file differs: {rel}")

    replay_report["temporary_copy_removed_after_test"] = True
    replay_report["manual_post_edit_required"] = False
    replay_report["hard_pass"] = not failures
    replay_report["submission_ready"] = not failures
    (qc_dir / "Figure1_real_regeneration_test.json").write_text(
        json.dumps(replay_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(replay_report, indent=2, ensure_ascii=False))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
