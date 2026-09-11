"""Generate CP-PAW structure files for a constrained O to T path.

Takes the structure file of the relaxed octahedral calculation as a template,
moves the hydrogen atom to ``R_H(x) = O + x (T - O)`` for every requested ``x``,
appends the translation-invariant linear constraint block and writes one
``.strc`` per point. The control files for the two stages are copied alongside.
This reproduces the inputs under ``calculations/04_path_octa_to_tetra``.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np

from pdhdiff.cppaw_io import parse_strc
from pdhdiff.structure import (
    TRIANGLE_LABELS,
    constraint_block,
    h_position_along_path,
    reaction_coordinate,
)

ROOT = Path(__file__).resolve().parents[1]
_H_LINE = re.compile(r"(!ATOM\s+NAME='H_01'\s+R=)\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+(\s*!END)")


def strip_constraints(text: str) -> str:
    """Remove an existing ``!CONSTRAINTS ... !END !END`` block from a structure file."""
    out, depth, skipping = [], 0, False
    for line in text.splitlines():
        if not skipping and line.strip().upper().startswith("!CONSTRAINTS"):
            skipping, depth = True, 1
            continue
        if skipping:
            depth += len(re.findall(r"!(?!END)[A-Z]+", line.upper()))
            depth -= line.upper().count("!END")
            if depth <= 0:
                skipping = False
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def make_structure(template: str, x: float) -> str:
    """Return the structure file text for reaction coordinate ``x``."""
    r = h_position_along_path(x)
    body = strip_constraints(template)
    body = _H_LINE.sub(lambda m: f"{m.group(1)} {r[0]:.4f} {r[1]:.4f} {r[2]:.4f}{m.group(2)}", body)
    block = constraint_block("H_01", TRIANGLE_LABELS, move=False)
    # insert the block before the closing !END of !STRUCTURE (the line before !EOB)
    lines = body.rstrip().splitlines()
    idx = max(i for i, ln in enumerate(lines) if ln.strip().upper() == "!END")
    lines[idx:idx] = ["", *block.splitlines(), ""]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--template",
        type=Path,
        default=ROOT
        / "calculations"
        / "03_h_octahedral"
        / "stage2_relaxation"
        / "pd_super_octa.strc",
    )
    ap.add_argument(
        "--stage1-cntl",
        type=Path,
        default=ROOT
        / "calculations"
        / "03_h_octahedral"
        / "stage1_wavefunctions"
        / "pd_super_octa.cntl",
    )
    ap.add_argument(
        "--stage2-cntl",
        type=Path,
        default=ROOT
        / "calculations"
        / "03_h_octahedral"
        / "stage2_relaxation"
        / "pd_super_octa.cntl",
    )
    ap.add_argument("--out", type=Path, default=ROOT / "scratch" / "path_inputs")
    ap.add_argument("--points", type=int, default=11, help="number of equally spaced values of x")
    args = ap.parse_args()

    template = args.template.read_text(encoding="utf-8")
    for x in np.linspace(0, 1, args.points):
        d = args.out / f"g_{x:.1f}"
        (d / "stage1_wavefunctions").mkdir(parents=True, exist_ok=True)
        (d / "stage2_relaxation").mkdir(parents=True, exist_ok=True)
        text = make_structure(template, float(x))
        for stage, cntl in (
            ("stage1_wavefunctions", args.stage1_cntl),
            ("stage2_relaxation", args.stage2_cntl),
        ):
            (d / stage / "pd_super_octa.strc").write_text(text, encoding="utf-8")
            shutil.copy(cntl, d / stage / "pd_super_octa.cntl")
        strc = parse_strc(d / "stage2_relaxation" / "pd_super_octa.strc")
        tri = np.array([strc.positions[strc.index(n)] for n in TRIANGLE_LABELS])
        g = reaction_coordinate(strc.positions[strc.index("H_01")], tri)
        print(
            f"x = {x:.2f}: wrote {d}, g(R) = {g:.4f}, constraints = {len(strc.linear_constraints)}"
        )


if __name__ == "__main__":
    main()
