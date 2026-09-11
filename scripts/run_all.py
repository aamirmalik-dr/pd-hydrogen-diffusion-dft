"""Run the full post-processing chain on the committed calculations."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [
    "extract_profile.py",
    "analyze_diffusion.py",
    "plot_electronic_structure.py",
    "make_path_figures.py",
]


def main() -> None:
    for step in STEPS:
        print(f"\n=== {step} ===")
        sys.argv = [step]
        runpy.run_path(str(HERE / step), run_name="__main__")


if __name__ == "__main__":
    main()
