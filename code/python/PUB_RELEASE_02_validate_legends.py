#!/usr/bin/env python3
"""Check journal title/legend limits for the locked Figure 2-S6 legends."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def count_words(text: str) -> int:
    return len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*|[α-ωΑ-Ω²≥≤−]+", text))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legend", required=True, type=Path)
    parser.add_argument("--qc-root", required=True, type=Path)
    args = parser.parse_args()
    text = args.legend.read_text(encoding="utf-8")
    matches = list(re.finditer(r"^## (Figure (?:[2-5]|S[1-6])\. .+)$", text, re.MULTILINE))
    if len(matches) != 10:
        raise ValueError(f"Expected 10 figure legend headings, found {len(matches)}")

    rows = []
    for index, match in enumerate(matches):
        heading = match.group(1).strip()
        body_start = match.end()
        body_end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        body = text[body_start:body_end].strip()
        title = re.sub(r"^Figure (?:[2-5]|S[1-6])\.\s*", "", heading)
        title_words = count_words(title)
        legend_words = count_words(body)
        status = "PASS" if title_words <= 15 and legend_words <= 300 and legend_words > 0 else "FAIL"
        rows.append({
            "figure": heading.split(".", 1)[0],
            "title": title,
            "title_words": title_words,
            "legend_words": legend_words,
            "title_limit": 15,
            "legend_limit": 300,
            "status": status,
        })

    args.qc_root.mkdir(parents=True, exist_ok=True)
    csv_path = args.qc_root / "FIGURE_LEGEND_WORD_COUNT_QC.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    failures = [row["figure"] for row in rows if row["status"] != "PASS"]
    report = [
        "# Figure legend journal-limit QC",
        "",
        f"- Legends checked: {len(rows)}.",
        "- Maximum title length: 15 words.",
        "- Maximum legend length: 300 words.",
        f"- Final status: {'PASS' if not failures else 'FAIL'}.",
    ]
    if failures:
        report.append("- Failed figures: " + ", ".join(failures) + ".")
    (args.qc_root / "FIGURE_LEGEND_WORD_COUNT_QC.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"{'PASS' if not failures else 'FAIL'}: title and legend word limits")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
