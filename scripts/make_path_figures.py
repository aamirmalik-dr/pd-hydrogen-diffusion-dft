"""Relaxed structures along the path: trajectory file, geometry and convergence figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from pdhdiff.constants import A_PD_ANG
from pdhdiff.cppaw_io import parse_prot
from pdhdiff.plotting import plot_convergence, plot_path_geometry, plot_path_structures
from pdhdiff.profile import STAGE2
from pdhdiff.structure import O_SITE, T_SITE, write_extxyz_frames

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calc", type=Path, default=ROOT / "calculations" / "04_path_octa_to_tetra")
    ap.add_argument("--figures", type=Path, default=ROOT / "figures")
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)

    dirs = sorted(d for d in args.calc.glob("g_*") if d.is_dir())
    frames, traces, rows = [], {}, []
    a = A_PD_ANG
    line_dir = (T_SITE - O_SITE) / np.linalg.norm(T_SITE - O_SITE)
    with open(args.results / "profile.csv", encoding="utf-8") as fh:
        prof = {row["label"]: row for row in csv.DictReader(fh)}
    for d in dirs:
        prot = parse_prot(next((d / STAGE2).glob("*.prot*")))
        symbols = [at.element for at in prot.atoms]
        pos = prot.positions()
        cell = float(prot.lattice[0, 0])
        g = float(prof[d.name]["g_nominal"])
        e_rel = float(prof[d.name]["energy_rel_mev"])
        frames.append(
            (
                symbols,
                pos,
                f'g={g:.1f} E_rel_meV={e_rel:.2f} Lattice="{cell} 0 0 0 {cell} 0 0 0 {cell}" '
                'pbc="T T T" Properties=species:S:1:pos:R:3',
            )
        )
        h = next(at.position for at in prot.atoms if at.element == "H")
        pd = np.array([at.position for at in prot.atoms if at.element == "Pd"])
        # offset of H from the ideal straight line O -> T
        rel = h / a - O_SITE
        offset = np.linalg.norm(rel - (rel @ line_dir) * line_dir) * a
        rows.append(
            (
                g,
                offset,
                float(prof[d.name]["max_pd_displacement_ang"]),
                float(prof[d.name]["triangle_pd_displacement_ang"]),
            )
        )
        if d.name in ("g_0.0", "g_0.6", "g_1.0"):
            traces[f"g = {g:.1f} ({len(prot.last_trace)} steps)"] = prot.last_trace
    write_extxyz_frames(args.results / "path_trajectory.xyz", frames)

    rows = np.array(rows)
    fig = plot_path_geometry(rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3])
    fig.tight_layout()
    fig.savefig(args.figures / "path_geometry.png")

    fig = plot_convergence(traces)
    fig.tight_layout()
    fig.savefig(args.figures / "relaxation_convergence.png")

    cage_frames = []
    for symbols, pos, _ in frames:
        pd = np.array([p for s, p in zip(symbols, pos) if s == "Pd"])
        h = next(p for s, p in zip(symbols, pos) if s == "H")
        cage_frames.append((pd, h))
    fig = plot_path_structures(cage_frames, a)
    fig.savefig(args.figures / "path_structures.png", bbox_inches="tight")
    print(f"wrote {len(frames)} frames to results/path_trajectory.xyz and three figures")
    for g, off, mx, tri in rows:
        print(
            f"g = {g:.1f}: H offset from line {off:.3f} A, "
            f"max Pd displacement {mx:.3f} A, triangle {tri:.3f} A"
        )


if __name__ == "__main__":
    main()
