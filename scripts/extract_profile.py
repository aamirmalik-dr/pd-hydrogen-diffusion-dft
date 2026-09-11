"""Extract the O to T energy profile and the site energies from the CP-PAW runs.

Reads every ``calculations/04_path_octa_to_tetra/g_*`` directory and the two
unconstrained site calculations, checks convergence flags and writes
``results/profile.csv`` and ``results/sites.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pdhdiff.constants import HARTREE_EV
from pdhdiff.cppaw_io import parse_prot
from pdhdiff.profile import STAGE2, collect_path, write_profile_csv

ROOT = Path(__file__).resolve().parents[1]


def site_summary(site_dir: Path) -> dict:
    """Final energy, iteration count and H position of an unconstrained site run."""
    prot = parse_prot(next((site_dir / STAGE2).glob("*.prot*")))
    h = next(a for a in prot.atoms if a.element == "H")
    pd = np.array([a.position for a in prot.atoms if a.element == "Pd"])
    d = pd - h.position
    d -= np.round(d / prot.lattice[0, 0]) * prot.lattice[0, 0]
    dist = np.sort(np.linalg.norm(d, axis=1))
    stage1 = parse_prot(next((site_dir / "stage1_wavefunctions").glob("*.prot*")))
    return {
        "energy_h": prot.final_energy,
        "energy_stage1_h": stage1.final_energy,
        "n_kpoints": prot.n_kpoints,
        "n_bands": prot.n_bands,
        "cutoff_ry": prot.cutoff_ry,
        "n_iterations": prot.n_iterations,
        "wallclock": prot.wallclock,
        "finished": prot.finished,
        "h_position_ang": [float(v) for v in h.position],
        "nearest_pd_distances_ang": [float(v) for v in dist[:6]],
        "max_force_mh_per_bohr": float(np.abs(prot.forces()).max()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calc", type=Path, default=ROOT / "calculations")
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    points = collect_path(args.calc / "04_path_octa_to_tetra")
    write_profile_csv(points, args.out / "profile.csv")

    sites = {
        "octahedral": site_summary(args.calc / "03_h_octahedral"),
        "tetrahedral": site_summary(args.calc / "03_h_tetrahedral"),
    }
    bulk = parse_prot(args.calc / "01_pd_primitive" / "pd.prot.gz")
    super_ = parse_prot(args.calc / "02_pd_supercell" / "pd_super.prot.gz")
    sites["bulk_primitive"] = {
        "energy_h": bulk.final_energy,
        "n_kpoints": bulk.n_kpoints,
        "n_bands": bulk.n_bands,
        "fermi_energy_ev": bulk.fermi_energy_ev,
        "n_iterations": bulk.n_iterations,
    }
    sites["bulk_supercell"] = {
        "energy_h": super_.final_energy,
        "energy_per_atom_h": super_.final_energy / 32,
        "n_kpoints": super_.n_kpoints,
        "n_bands": super_.n_bands,
        "fermi_energy_ev": super_.fermi_energy_ev,
        "n_iterations": super_.n_iterations,
    }
    e_o, e_t = sites["octahedral"]["energy_h"], sites["tetrahedral"]["energy_h"]
    sites["derived"] = {
        "e_t_minus_e_o_mev": (e_t - e_o) * HARTREE_EV * 1e3,
        "path_start_minus_octahedral_mev": (points[0].energy_h - e_o) * HARTREE_EV * 1e3,
        "path_end_minus_tetrahedral_mev": (points[-1].energy_h - e_t) * HARTREE_EV * 1e3,
        "supercell_minus_primitive_per_atom_mev": (super_.final_energy / 32 - bulk.final_energy)
        * HARTREE_EV
        * 1e3,
        "h_solution_energy_relative_to_free_atom_not_computed": True,
    }
    with open(args.out / "sites.json", "w", encoding="utf-8") as fh:
        json.dump(sites, fh, indent=2)

    print(
        f"{'g':>5} {'g(constraint)':>14} {'g(relaxed)':>11} {'E [H]':>14} "
        f"{'E-E_O [meV]':>12} {'steps':>6} done"
    )
    e0 = points[0].energy_h
    for p in points:
        print(
            f"{p.g_nominal:5.2f} {p.g_constraint:14.4f} {p.g_relaxed:11.4f} {p.energy_h:14.7f}"
            f" {(p.energy_h - e0) * HARTREE_EV * 1e3:12.1f} {p.n_iterations!s:>6} {p.finished}"
        )
    print(f"E_T - E_O (unconstrained) = {sites['derived']['e_t_minus_e_o_mev']:.2f} meV")
    print(
        "path end-point checks: "
        f"{sites['derived']['path_start_minus_octahedral_mev']:.3f} meV (O), "
        f"{sites['derived']['path_end_minus_tetrahedral_mev']:.3f} meV (T)"
    )
    print(f"wrote {args.out / 'profile.csv'} and {args.out / 'sites.json'}")


if __name__ == "__main__":
    main()
