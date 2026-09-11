"""Fit the energy profile, evaluate harmonic TST rates and the diffusion constant.

Inputs: ``results/profile.csv`` and ``results/sites.json`` from ``extract_profile.py``.
Outputs: ``results/metrics.json``, ``results/diffusion_vs_T.csv``,
``figures/energy_profile.png``, ``figures/arrhenius.png`` and the README hero.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from pdhdiff.constants import (
    A_PD_ANG,
    D_EXP_298K_M2_S,
    EA_EXP_EV,
    MASS_D_AMU,
    MASS_H_AMU,
    MASS_T_AMU,
)
from pdhdiff.plotting import hero_figure, plot_arrhenius, plot_energy_profile
from pdhdiff.profile import (
    fit_octahedral_well,
    fit_tetrahedral_well,
    locate_barrier,
    read_profile_csv,
)
from pdhdiff.structure import hop_geometry, path_length_ang
from pdhdiff.tst import (
    SiteModel,
    arrhenius_fit,
    attempt_frequency,
    diffusion_constant,
    force_constant_si,
    vibrational_quantum_mev,
    with_mass,
)

ROOT = Path(__file__).resolve().parents[1]


def d_at(model: SiteModel, temperature_k: float, geo) -> float:
    """Diffusion constant of ``model`` at one temperature with the primary conventions."""
    return float(diffusion_constant(model, [temperature_k], geo).d_m2_s[0])


def sensitivity(
    model: SiteModel, barrier, well_o, well_t, e_t_ev: float, length: float, geo
) -> dict:
    """Spread of D(298 K) over the alternative estimators of the fitted quantities.

    Barrier: local cubic (primary), cubic spline, highest calculated point.
    Curvatures: least-squares fits (primary) and the two-point estimates from the
    calculated point nearest to each minimum. The extreme combinations bound the
    range that the choice of estimator alone can produce.
    """
    barriers = {
        "local_cubic": barrier.e_ts_ev,
        "spline": barrier.e_ts_spline_ev,
        "raw_maximum": barrier.e_raw_max_ev,
    }
    k_o = {"least_squares": well_o.k_g_ev, "two_point": well_o.k_g_two_point_ev}
    k_t = {"least_squares": well_t.k_g_ev, "two_point": well_t.k_g_two_point_ev}

    def build(eb: float, ko: float, kt: float) -> SiteModel:
        return SiteModel(
            de_t_minus_o_ev=e_t_ev,
            ea_o_to_t_ev=eb,
            ea_t_to_o_ev=eb - e_t_ev,
            omega_o=attempt_frequency(force_constant_si(ko, length)),
            omega_t=attempt_frequency(force_constant_si(kt, length)),
        )

    single = {}
    for name, eb in barriers.items():
        single[f"barrier={name}"] = d_at(build(eb, well_o.k_g_ev, well_t.k_g_ev), 298.0, geo)
    for name, ko in k_o.items():
        single[f"k_O={name}"] = d_at(build(barrier.e_ts_ev, ko, well_t.k_g_ev), 298.0, geo)
    for name, kt in k_t.items():
        single[f"k_T={name}"] = d_at(build(barrier.e_ts_ev, well_o.k_g_ev, kt), 298.0, geo)
    grid = [
        d_at(build(eb, ko, kt), 298.0, geo)
        for eb in barriers.values()
        for ko in k_o.values()
        for kt in k_t.values()
    ]
    return {
        "note": "D(298 K) for each alternative estimator of one fitted quantity, the others at "
        "their primary values, and the min/max over all combinations.",
        "barrier_estimates_meV": {k: v * 1e3 for k, v in barriers.items()},
        "k_O_estimates_eV_per_g2": k_o,
        "k_T_estimates_eV_per_g2": k_t,
        "d_298K_single_factor_m2_s": single,
        "d_298K_min_m2_s": min(grid),
        "d_298K_max_m2_s": max(grid),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--figures", type=Path, default=ROOT / "figures")
    ap.add_argument("--tmin", type=float, default=200.0)
    ap.add_argument("--tmax", type=float, default=1000.0)
    args = ap.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)

    g, e = read_profile_csv(args.results / "profile.csv")
    with open(args.results / "sites.json", encoding="utf-8") as fh:
        sites = json.load(fh)
    de_uncon = sites["derived"]["e_t_minus_e_o_mev"] / 1e3

    barrier = locate_barrier(g, e)
    well_o = fit_octahedral_well(g, e)
    well_t = fit_tetrahedral_well(g, e)
    length = path_length_ang(A_PD_ANG)

    k_o = force_constant_si(well_o.k_g_ev, length)
    k_t = force_constant_si(well_t.k_g_ev, length)
    w_o = attempt_frequency(k_o, MASS_H_AMU)
    w_t = attempt_frequency(k_t, MASS_H_AMU)

    ea1 = barrier.e_ts_ev
    ea2 = barrier.e_ts_ev - e[-1]
    model = SiteModel(
        de_t_minus_o_ev=float(e[-1]),
        ea_o_to_t_ev=ea1,
        ea_t_to_o_ev=ea2,
        omega_o=w_o,
        omega_t=w_t,
    )
    geo = hop_geometry(A_PD_ANG)
    temps = np.linspace(args.tmin, args.tmax, 161)
    primary = diffusion_constant(
        model, temps, geo, population_mode="harmonic1d", einstein_half=True
    )
    alternatives = {
        "populations with 3D harmonic prefactor": diffusion_constant(
            model, temps, geo, "harmonic3d", True
        ),
        "populations from site energies only": diffusion_constant(
            model, temps, geo, "boltzmann", True
        ),
        "3D prefactor, no Einstein 1/2 (seminar convention)": diffusion_constant(
            model, temps, geo, "harmonic3d", False
        ),
    }
    ea_fit, d0_fit = arrhenius_fit(temps, primary.d_m2_s)

    def at(res, t):
        return float(np.interp(t, res.temperature_k, res.d_m2_s))

    table_t = [200, 300, 400, 600, 800, 1000]
    metrics = {
        "inputs": {
            "lattice_constant_ang": A_PD_ANG,
            "path_length_ang": length,
            "hydrogen_mass_amu": MASS_H_AMU,
            "n_path_points": int(len(g)),
        },
        "profile": {
            "e_t_minus_e_o_path_mev": float(e[-1] * 1e3),
            "e_t_minus_e_o_unconstrained_mev": float(de_uncon * 1e3),
            "barrier": asdict(barrier),
            "e_act_o_to_t_mev": ea1 * 1e3,
            "e_act_t_to_o_mev": ea2 * 1e3,
            "well_octahedral": asdict(well_o),
            "well_tetrahedral": asdict(well_t),
        },
        "vibrations": {
            "k_octahedral_N_per_m": k_o,
            "k_tetrahedral_N_per_m": k_t,
            "omega_octahedral_rad_s": w_o,
            "omega_tetrahedral_rad_s": w_t,
            "nu_octahedral_Hz": w_o / (2 * np.pi),
            "nu_tetrahedral_Hz": w_t / (2 * np.pi),
            "hbar_omega_octahedral_meV": vibrational_quantum_mev(w_o),
            "hbar_omega_tetrahedral_meV": vibrational_quantum_mev(w_t),
        },
        "geometry": {
            "t_neighbours_of_o": geo.n_t_around_o,
            "o_neighbours_of_t": geo.n_o_around_t,
            "geometric_factor_o_ang2_isotropic": float(np.trace(geo.factor_o_ang2) / 3),
            "geometric_factor_t_ang2_isotropic": float(np.trace(geo.factor_t_ang2) / 3),
            "geometric_factor_o_over_a2": float(np.trace(geo.factor_o_ang2) / 3 / A_PD_ANG**2),
            "geometric_factor_t_over_a2": float(np.trace(geo.factor_t_ang2) / 3 / A_PD_ANG**2),
        },
        "diffusion": {
            "convention": "D = 1/2 sum_i P_i sum_j Gamma_ji |r_j - r_i|^2 / 3 (Einstein relation); "
            "sublattice populations from the 1D harmonic wells, consistent with the 1D rates",
            "arrhenius_fit_ea_eV": ea_fit,
            "arrhenius_fit_d0_m2_s": d0_fit,
            "d_298K_m2_s": at(primary, 298.0),
            "d_298K_experiment_m2_s": D_EXP_298K_M2_S,
            "ratio_calc_over_experiment_298K": at(primary, 298.0) / D_EXP_298K_M2_S,
            "ea_experiment_eV_quoted": EA_EXP_EV,
            "detailed_balance_ratio_300K": float(
                np.interp(300.0, primary.temperature_k, primary.detailed_balance_ratio)
            ),
            "table": [
                {
                    "T_K": t,
                    "gamma_t_from_o_per_s": float(
                        np.interp(t, primary.temperature_k, primary.gamma_t_from_o)
                    ),
                    "gamma_o_from_t_per_s": float(
                        np.interp(t, primary.temperature_k, primary.gamma_o_from_t)
                    ),
                    "p_tet": float(np.interp(t, primary.temperature_k, primary.p_tet)),
                    "d_m2_s": at(primary, t),
                }
                for t in table_t
            ],
            "alternatives_d_298K_m2_s": {k: at(v, 298.0) for k, v in alternatives.items()},
        },
        "sensitivity": sensitivity(model, barrier, well_o, well_t, float(e[-1]), length, geo),
        "isotopes_classical": {
            "note": "Classical harmonic TST with the same profile: only the attempt "
            "frequencies change, by sqrt(m_H / m). Zero-point energy and tunnelling, "
            "which dominate the measured isotope effect, are not included.",
            **{
                label: {
                    "mass_amu": mass,
                    "d_298K_m2_s": d_at(with_mass(model, mass), 298.0, geo),
                    "d_298K_over_h": d_at(with_mass(model, mass), 298.0, geo)
                    / d_at(model, 298.0, geo),
                }
                for label, mass in (("H", MASS_H_AMU), ("D", MASS_D_AMU), ("T", MASS_T_AMU))
            },
        },
    }
    with open(args.results / "metrics.json", "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, indent=2)

    with open(args.results / "diffusion_vs_T.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["T_K", "gamma_t_from_o_per_s", "gamma_o_from_t_per_s", "p_oct", "p_tet", "d_m2_s"]
        )
        for i, t in enumerate(temps):
            w.writerow(
                [
                    f"{t:.1f}",
                    f"{primary.gamma_t_from_o[i]:.4e}",
                    f"{primary.gamma_o_from_t[i]:.4e}",
                    f"{primary.p_oct[i]:.6f}",
                    f"{primary.p_tet[i]:.6f}",
                    f"{primary.d_m2_s[i]:.4e}",
                ]
            )

    prof_kwargs = dict(
        g=g, e_ev=e, barrier=barrier, well_o=well_o, well_t=well_t, e_t_unconstrained_ev=de_uncon
    )
    arr_kwargs = dict(result=primary, ea_fit_ev=ea_fit, d0_fit=d0_fit, alt=alternatives)
    fig = plot_energy_profile(**prof_kwargs)
    fig.tight_layout()
    fig.savefig(args.figures / "energy_profile.png")
    fig = plot_arrhenius(**arr_kwargs)
    fig.tight_layout()
    fig.savefig(args.figures / "arrhenius.png")
    fig = hero_figure(prof_kwargs, arr_kwargs)
    fig.savefig(args.figures / "hero_profile_arrhenius.png")

    print(
        f"barrier: g = {barrier.g_ts:.3f}, E = {barrier.e_ts_ev*1e3:.1f} meV "
        f"(spline {barrier.e_ts_spline_ev*1e3:.1f} meV at g = {barrier.g_ts_spline:.3f}; "
        f"raw max {barrier.e_raw_max_ev*1e3:.1f})"
    )
    print(
        f"E_act O->T = {ea1*1e3:.1f} meV, T->O = {ea2*1e3:.1f} meV, E_T - E_O = {e[-1]*1e3:.1f} meV"
    )
    print(
        f"k_O = {well_o.k_g_ev:.3f} eV/g^2 ({k_o:.2f} N/m), "
        f"hbar w_O = {vibrational_quantum_mev(w_o):.1f} meV"
    )
    print(
        f"k_T = {well_t.k_g_ev:.3f} eV/g^2 ({k_t:.2f} N/m), "
        f"hbar w_T = {vibrational_quantum_mev(w_t):.1f} meV"
    )
    print(
        f"D(298 K) = {at(primary, 298):.3e} m^2/s  (experiment {D_EXP_298K_M2_S:.3e}), "
        f"Arrhenius fit E_a = {ea_fit:.4f} eV, D0 = {d0_fit:.3e} m^2/s"
    )
    for name, res in alternatives.items():
        print(f"  alternative {name}: D(298 K) = {at(res, 298):.3e} m^2/s")
    print(
        "detailed-balance ratio at 300 K: "
        f"{metrics['diffusion']['detailed_balance_ratio_300K']:.3f}"
    )


if __name__ == "__main__":
    main()
