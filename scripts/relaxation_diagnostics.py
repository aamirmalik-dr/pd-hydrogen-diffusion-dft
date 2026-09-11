"""How well relaxed are the final structures? Diagnostics from the protocol files.

CP-PAW's autopilot stops a damped relaxation on the energy and the ionic kinetic
energy, not on the forces. The periodic atom lists in the protocol carry the
forces that propagated the atoms, so the last one shows what was left when the
run stopped, and the sequence shows the force on a cage atom against its
distance to the hydrogen while the cage opens. This script collects, for every
stage-2 run:

* the largest residual force on the constrained atoms and on the free atoms,
* the energy drift over the last 30 steps and the ionic temperature at the end,
* how far the final energy lies above the lowest energy seen with cold electrons,
* a harmonic upper bound on the energy still to be gained, using the radial
  stiffness of the first-shell atoms read off the tetrahedral run itself,
* for the path points, the Lagrange multiplier printed by CP-PAW against the
  slope of the fitted profile (they agree only at a constrained minimum).

Outputs ``results/relaxation_diagnostics.json`` and
``figures/relaxation_diagnostics.png``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from pdhdiff.constants import A_PD_ANG, BOHR_ANG, HARTREE_EV
from pdhdiff.cppaw_io import parse_prot, parse_reports
from pdhdiff.plotting import plot_relaxation_diagnostics
from pdhdiff.profile import STAGE2, read_profile_csv
from pdhdiff.structure import TRIANGLE_LABELS

ROOT = Path(__file__).resolve().parents[1]
MH_BOHR_TO_EV_ANG = 1e-3 * HARTREE_EV / BOHR_ANG
COLD_ELECTRONS_H = 2e-4  # fictitious wave-function kinetic energy below which E(RHO) is trusted
FIRST_SHELL_ANG = 2.2  # Pd atoms closer than this to H form the cage


def _minimum_image(d: np.ndarray, cell: float) -> np.ndarray:
    return d - np.round(d / cell) * cell


def radial_force_history(prot_path: Path, atom: str, h_name: str = "H_01") -> list[dict]:
    """Distance to H and radial force on ``atom`` at every printed report."""
    rows = []
    for rep in parse_reports(prot_path):
        p, f = rep.atom(atom).position, rep.atom(atom).force
        cell = parse_prot(prot_path).lattice[0, 0]
        d = _minimum_image(p - rep.atom(h_name).position, cell)
        r = float(np.linalg.norm(d))
        rows.append(
            {
                "nfi": rep.nfi,
                "total_energy_h": rep.total_energy_h,
                "distance_to_h_ang": r,
                "radial_force_outward_mh_bohr": float(f @ (d / r)),
                "force_norm_mh_bohr": float(np.linalg.norm(f)),
            }
        )
    return rows


def run_summary(prot_path: Path, constrained: tuple[str, ...], h_name: str = "H_01") -> dict:
    """Residual forces, energy drift and shell forces of one relaxation run."""
    prot = parse_prot(prot_path)
    cell = prot.lattice[0, 0]
    h = prot.atom(h_name).position
    free, shell, cons = [], [], []
    for a in prot.atoms:
        fn = float(np.linalg.norm(a.force))
        if a.name in constrained:
            cons.append(fn)
        else:
            free.append(fn)
        if a.element == "Pd":
            r = float(np.linalg.norm(_minimum_image(a.position - h, cell)))
            if r < FIRST_SHELL_ANG:
                shell.append((a.name, r, fn, a.name in constrained))
    tr = prot.last_trace
    cold = tr[tr[:, 3] < COLD_ELECTRONS_H]
    e_min_cold = float(cold[:, 4].min()) if len(cold) else float("nan")
    return {
        "n_steps_last_run": int(len(tr)),
        "final_energy_h": prot.final_energy,
        "energy_drift_last_30_steps_mev": float((tr[-1, 4] - tr[-30, 4]) * HARTREE_EV * 1e3),
        "final_energy_above_lowest_cold_trace_energy_mev": float(
            (tr[-1, 4] - e_min_cold) * HARTREE_EV * 1e3
        ),
        "ionic_temperature_end_k": float(tr[-1, 2]),
        "max_force_free_atoms_mh_bohr": max(free) if free else None,
        "max_force_constrained_atoms_mh_bohr": max(cons) if cons else None,
        "first_shell": [
            {"atom": n, "distance_to_h_ang": r, "force_mh_bohr": fn, "constrained": c}
            for n, r, fn, c in sorted(shell, key=lambda x: x[1])
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calc", type=Path, default=ROOT / "calculations")
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--figures", type=Path, default=ROOT / "figures")
    args = ap.parse_args()

    runs = {
        "octahedral": (args.calc / "03_h_octahedral" / STAGE2 / "pd_super_octa.prot.gz", ()),
        "tetrahedral": (args.calc / "03_h_tetrahedral" / STAGE2 / "pd_super_tetra.prot.gz", ()),
    }
    for d in sorted((args.calc / "04_path_octa_to_tetra").glob("g_*")):
        runs[d.name] = (d / STAGE2 / "pd_super_octa.prot.gz", ("H_01", *TRIANGLE_LABELS))

    out: dict = {"units": "forces in milli-Hartree per bohr (1 mH/bohr = 0.0514 eV/A)"}
    out["runs"] = {name: run_summary(p, c) for name, (p, c) in runs.items()}

    # force against distance for one cage atom of the tetrahedral run: the stiffness
    hist_t = radial_force_history(runs["tetrahedral"][0], "PD03")
    hist_o = radial_force_history(runs["octahedral"][0], "PD04")
    r = np.array([h["distance_to_h_ang"] for h in hist_t])
    f = np.array([h["radial_force_outward_mh_bohr"] for h in hist_t])
    slope, intercept = np.polyfit(r, f, 1)  # mH/bohr per angstrom, negative
    k_radial_ev_ang2 = -slope * MH_BOHR_TO_EV_ANG  # eV/A^2
    r_zero = -intercept / slope
    out["cage_stiffness"] = {
        "atom": "PD03 of the tetrahedral run",
        "history": hist_t,
        "history_octahedral_PD04": hist_o,
        "k_radial_ev_per_ang2": float(k_radial_ev_ang2),
        "zero_force_distance_ang": float(r_zero),
        "note": "Linear fit of the outward radial force against the Pd-H distance over the "
        "printed reports. The surroundings relax at the same time, so this is the soft, "
        "relaxed response and the harmonic bound below is an upper estimate.",
    }

    # harmonic upper bound on the unrelaxed energy: sum over first-shell atoms of F^2 / 2k
    for name, summary in out["runs"].items():
        e_bound = sum(
            (s["force_mh_bohr"] * MH_BOHR_TO_EV_ANG) ** 2 / (2 * k_radial_ev_ang2)
            for s in summary["first_shell"]
        )
        summary["harmonic_upper_bound_unrelaxed_energy_mev"] = float(e_bound * 1e3)

    # Lagrange multiplier against the slope of the profile
    g, e = read_profile_csv(args.results / "profile.csv")
    spline = CubicSpline(g, e)
    a_bohr = A_PD_ANG / BOHR_ANG
    lam_rows = []
    for gi in g:
        prot = parse_prot(runs[f"g_{gi:.1f}"][0])
        lam = -prot.constraint_values[-1][1]  # multiplier mu, sign as dE/dVAL
        slope_h_per_bohr = float(spline(gi, 1)) / HARTREE_EV / a_bohr
        lam_rows.append(
            {
                "g": float(gi),
                "multiplier_h_per_bohr": lam,
                "profile_slope_h_per_bohr": slope_h_per_bohr,
                "ratio": lam / slope_h_per_bohr if abs(slope_h_per_bohr) > 1e-5 else None,
            }
        )
    ti = -a_bohr * np.concatenate(
        [
            [0.0],
            np.cumsum(
                0.5
                * (np.diff(g))
                * (
                    np.array([r["multiplier_h_per_bohr"] for r in lam_rows])[1:]
                    + np.array([r["multiplier_h_per_bohr"] for r in lam_rows])[:-1]
                )
            ),
        ]
    )
    out["generalized_force_check"] = {
        "note": "At a constrained minimum the Lagrange multiplier equals dE/dVAL. The shortfall "
        "of the multiplier against the slope of the energy profile measures how far the "
        "unconstrained atoms are from their minimum along the path.",
        "rows": lam_rows,
        "thermodynamic_integration_of_multiplier_mev": [float(v * HARTREE_EV * 1e3) for v in -ti],
        "direct_profile_mev": [float(v * 1e3) for v in e],
    }

    args.results.mkdir(parents=True, exist_ok=True)
    with open(args.results / "relaxation_diagnostics.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    fig = plot_relaxation_diagnostics(hist_t, hist_o, k_radial_ev_ang2, r_zero, lam_rows)
    fig.savefig(args.figures / "relaxation_diagnostics.png")

    print(
        f"cage stiffness (PD03, tetrahedral run): {k_radial_ev_ang2:.1f} eV/A^2, "
        f"zero force at {r_zero:.3f} A"
    )
    head = (
        f"{'run':12s} {'steps':>5} {'drift30[meV]':>12} {'T_end[K]':>8} "
        f"{'Fmax free':>9} {'Fmax cons':>9} {'E-Emin_cold':>11} {'bound[meV]':>10}"
    )
    print(head)
    for name, s in out["runs"].items():
        fc = s["max_force_constrained_atoms_mh_bohr"]
        fc = fc if fc is not None else float("nan")
        print(
            f"{name:12s} {s['n_steps_last_run']:5d} {s['energy_drift_last_30_steps_mev']:12.3f}"
            f" {s['ionic_temperature_end_k']:8.1f} {s['max_force_free_atoms_mh_bohr']:9.2f}"
            f" {fc:9.2f} {s['final_energy_above_lowest_cold_trace_energy_mev']:11.2f}"
            f" {s['harmonic_upper_bound_unrelaxed_energy_mev']:10.1f}"
        )
    ratios = [f"{r['ratio']:.2f}" if r["ratio"] else "-" for r in lam_rows]
    print("multiplier / profile slope:", ratios)


if __name__ == "__main__":
    main()
