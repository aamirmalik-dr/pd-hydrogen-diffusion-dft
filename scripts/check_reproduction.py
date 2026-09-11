"""Compare freshly regenerated result files with the committed versions.

Run after ``scripts/run_all.py``. Every numeric value in ``results/metrics.json``,
``results/sites.json`` and ``results/electronic_structure.json`` and every row of
``results/profile.csv`` and ``results/diffusion_vs_T.csv`` must agree with the
version at ``HEAD`` to a relative tolerance of 1e-6. Figures are not compared:
their bytes depend on the Matplotlib and font versions.
"""

from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-6
ATOL = 1e-12


def committed(path: str) -> str:
    return subprocess.run(
        ["git", "show", f"HEAD:{path}"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout


def compare(a, b, where: str, problems: list[str]) -> None:
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            problems.append(f"{where}: keys differ {sorted(set(a) ^ set(b))}")
        for k in a.keys() & b.keys():
            compare(a[k], b[k], f"{where}.{k}", problems)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            problems.append(f"{where}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{where}[{i}]", problems)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if abs(a - b) > ATOL + RTOL * max(abs(a), abs(b)):
            problems.append(f"{where}: {a!r} vs committed {b!r}")
    elif a != b:
        problems.append(f"{where}: {a!r} vs committed {b!r}")


def compare_csv(new: str, old: str, where: str, problems: list[str]) -> None:
    rows_new = list(csv.DictReader(io.StringIO(new)))
    rows_old = list(csv.DictReader(io.StringIO(old)))
    if len(rows_new) != len(rows_old):
        problems.append(f"{where}: {len(rows_new)} rows vs committed {len(rows_old)}")
        return
    for i, (rn, ro) in enumerate(zip(rows_new, rows_old)):
        for k in rn:
            try:
                compare(float(rn[k]), float(ro[k]), f"{where}[{i}].{k}", problems)
            except (ValueError, TypeError):
                compare(rn[k], ro.get(k), f"{where}[{i}].{k}", problems)


def main() -> int:
    problems: list[str] = []
    for name in (
        "metrics.json",
        "sites.json",
        "electronic_structure.json",
        "relaxation_diagnostics.json",
    ):
        rel = f"results/{name}"
        new = json.loads((ROOT / rel).read_text(encoding="utf-8"))
        old = json.loads(committed(rel))
        compare(new, old, name, problems)
    for name in ("profile.csv", "diffusion_vs_T.csv"):
        rel = f"results/{name}"
        compare_csv((ROOT / rel).read_text(encoding="utf-8"), committed(rel), name, problems)
    if problems:
        print("reproduction check FAILED:")
        for p in problems[:50]:
            print("  " + p)
        return 1
    print("reproduction check passed: regenerated results match the committed files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
