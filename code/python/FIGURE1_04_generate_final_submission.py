#!/usr/bin/env python3
"""Compatibility entry point for the self-contained accepted Figure 1 generator."""

from pathlib import Path
import runpy
import sys


def main() -> None:
    repository_root = Path(__file__).resolve().parents[2]
    package_root = repository_root / "figure1_reproducibility"
    target = package_root / "code" / "generate_Figure1.py"
    args = sys.argv[1:]
    if not any(arg == "--package-root" or arg.startswith("--package-root=") for arg in args):
        args = ["--package-root", str(package_root), *args]
    sys.argv = [str(target), *args]
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()

