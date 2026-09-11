"""Profile fitting and transition-state-theory tests, including regression against results/."""

import json

import numpy as np
import pytest

from pdhdiff.constants import A_PD_ANG, HARTREE_EV, KB_EV, MASS_H_AMU
from pdhdiff.profile import (
    collect_path,
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
    jump_rate,
    site_populations,
)


@pytest.fixture(scope="module")
def profile(results_root):
    return read_profile_csv(results_root / "profile.csv")


@pytest.fixture(scope="module")
def metrics(results_root):
    with open(results_root / "metrics.json", encoding="utf-8") as fh:
        return json.load(fh)


def test_collect_path_matches_committed_csv(path_root, profile):
    points = collect_path(path_root)
    g, e = profile
    assert len(points) == 11
    assert np.allclose([p.g_nominal for p in points], g)
    e0 = points[0].energy_h
    assert np.allclose([(p.energy_h - e0) * HARTREE_EV for p in points], e, atol=1e-6)
    assert all(p.finished for p in points)
    assert all(abs(p.g_constraint - p.g_nominal) < 2e-3 for p in points)


def test_synthetic_profile_fits():
    g = np.linspace(0, 1, 11)
    e = 0.5 * 1.0 * g**2 - 0.3 * g**4  # even well at O, known curvature 1.0
    fit = fit_octahedral_well(g, e, g_max=0.3)
    assert abs(fit.k_g_ev - 1.0) < 1e-9 and abs(fit.anharmonic_coeff + 0.3) < 1e-9
    d = g - 1
    e2 = 0.08 + 0.5 * 8.0 * d**2 + 10.0 * d**3
    fit_t = fit_tetrahedral_well(g, e2, g_min=0.7)
    assert abs(fit_t.k_g_ev - 8.0) < 1e-9 and abs(fit_t.anharmonic_coeff - 10.0) < 1e-9
    e3 = -((g - 0.63) ** 2) + 0.2
    bar = locate_barrier(g, e3)
    assert abs(bar.g_ts - 0.63) < 2e-3 and abs(bar.e_ts_ev - 0.2) < 1e-4


def test_barrier_and_wells_regression(profile, metrics):
    g, e = profile
    bar = locate_barrier(g, e)
    assert abs(bar.e_ts_ev * 1e3 - metrics["profile"]["barrier"]["e_ts_ev"] * 1e3) < 1e-6
    assert 0.6 < bar.g_ts < 0.7
    assert abs(bar.e_ts_ev - bar.e_ts_spline_ev) < 0.003  # two estimates agree within 3 meV
    wo = fit_octahedral_well(g, e)
    wt = fit_tetrahedral_well(g, e)
    assert abs(wo.k_g_ev - metrics["profile"]["well_octahedral"]["k_g_ev"]) < 1e-9
    assert abs(wt.k_g_ev - metrics["profile"]["well_tetrahedral"]["k_g_ev"]) < 1e-9
    assert wo.k_g_ev > 0 and wt.k_g_ev > wo.k_g_ev


def test_attempt_frequency_units():
    # k = 1 eV/A^2 on a 1 A path with the proton mass: omega = sqrt(16.02 N/m / 1.673e-27 kg)
    k = force_constant_si(1.0, 1.0)
    assert abs(k - 16.0218) < 1e-3
    w = attempt_frequency(k, MASS_H_AMU)
    assert abs(w / 9.78e13 - 1) < 0.01


def test_rates_and_detailed_balance():
    length = path_length_ang(A_PD_ANG)
    w_o = attempt_frequency(force_constant_si(1.0, length))
    w_t = attempt_frequency(force_constant_si(8.0, length))
    model = SiteModel(0.08, 0.20, 0.12, w_o, w_t)
    t = np.array([200.0, 300.0, 600.0])
    assert np.allclose(jump_rate(w_o, 0.2, t), w_o / (2 * np.pi) * np.exp(-0.2 / (KB_EV * t)))
    res = diffusion_constant(model, t, hop_geometry(), population_mode="harmonic1d")
    assert np.allclose(res.detailed_balance_ratio, 1.0)
    assert np.allclose(res.term_o_m2_s, res.term_t_m2_s)  # both fluxes equal at detailed balance
    p_o, p_t = site_populations(model, t, "boltzmann")
    assert np.all(p_o + p_t == pytest.approx(1.0))
    assert np.all(np.diff(res.d_m2_s) > 0)
    grid = np.linspace(200, 1000, 50)
    ea, d0 = arrhenius_fit(grid, diffusion_constant(model, grid).d_m2_s)
    assert 0.19 < ea < 0.21 and d0 > 0


def test_diffusion_regression(metrics):
    d = metrics["diffusion"]
    assert 1e-10 < d["d_298K_m2_s"] < 1e-9
    assert 0.19 < d["arrhenius_fit_ea_eV"] < 0.21
    assert abs(d["detailed_balance_ratio_300K"] - 1.0) < 1e-6
    table = {row["T_K"]: row["d_m2_s"] for row in d["table"]}
    assert table[200] < table[300] < table[600] < table[1000]
