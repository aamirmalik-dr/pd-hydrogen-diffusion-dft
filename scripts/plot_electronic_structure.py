"""Plot the projected DOS and the band structure of bulk fcc Pd from the CP-PAW output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pdhdiff.cppaw_io import parse_bands, parse_dos, parse_prot
from pdhdiff.plotting import plot_bands, plot_dos

ROOT = Path(__file__).resolve().parents[1]
K_LABELS = ["$\\Gamma$", "X", "W", "L", "$\\Gamma$", "K", "X"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calc", type=Path, default=ROOT / "calculations" / "01_pd_primitive")
    ap.add_argument("--figures", type=Path, default=ROOT / "figures")
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)

    prot = parse_prot(args.calc / "pd.prot.gz")
    ef = prot.fermi_energy_ev
    curves, occ_curves = {}, {}
    for key in ("s", "p", "d"):
        dos = parse_dos(args.calc / "dos" / f"{key}.dos.gz")
        curves[key] = dos.total
        occ_curves[key] = dos.occupied
        energy = dos.energy_ev
    total = parse_dos(args.calc / "dos" / "total.dos.gz")
    fig = plot_dos(energy, total.total, curves, ef)
    fig.tight_layout()
    fig.savefig(args.figures / "pd_dos.png")

    segments = parse_bands(args.calc / "pd_bands.dat")
    fig = plot_bands(segments, ef, K_LABELS)
    fig.tight_layout()
    fig.savefig(args.figures / "pd_bands.png")

    # numbers for the results file: d-band edges and integrated occupations
    de = float(np.median(np.diff(energy)))
    e_rel = energy - ef
    occ_total = float(np.trapezoid(total.occupied, energy))
    occ = {k: 2 * float(np.trapezoid(c, energy)) for k, c in occ_curves.items()}
    d_occ = occ_curves["d"]
    nonzero = e_rel[d_occ > 0.05 * d_occ.max()]
    summary = {
        "note": "DOS files are per spin; electron counts below are multiplied by 2 for the "
        "non-spin-polarised run. Projections are sphere weights and do not add up to the total.",
        "fermi_energy_ev": ef,
        "dos_at_fermi_states_per_ev_both_spins": 2 * float(np.interp(0.0, e_rel, total.total)),
        "occupied_electrons_total_both_spins": 2 * occ_total,
        "occupied_electrons_by_projection_both_spins": occ,
        "d_band_lower_edge_rel_ef_ev": float(nonzero.min()),
        "d_band_upper_edge_rel_ef_ev": float(nonzero.max()),
        "energy_grid_spacing_ev": float(de),
        "n_kpoints": prot.n_kpoints,
        "bands_n_segments": len(segments),
        "bands_n_bands": int(segments[0].energies_ev.shape[1]),
        "band_energies_at_gamma_rel_ef_ev": [float(v) for v in segments[0].energies_ev[0] - ef],
    }
    with open(args.results / "electronic_structure.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
