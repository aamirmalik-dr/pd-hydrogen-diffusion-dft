"""Tests of the periodic atom-list parser and of the committed relaxation diagnostics."""

import json

import numpy as np

from pdhdiff.cppaw_io import parse_prot, parse_reports


def test_reports_of_the_tetrahedral_run(calc_root):
    path = calc_root / "03_h_tetrahedral" / "stage2_relaxation" / "pd_super_tetra.prot.gz"
    reports = parse_reports(path)
    assert len(reports) == 6
    assert [r.nfi for r in reports] == [645, 700, 800, 900, 1000, 1022]
    energies = [r.total_energy_h for r in reports]
    assert energies == sorted(energies, reverse=True)[: len(energies)] or energies[0] > energies[-1]
    final = parse_prot(path)
    assert abs(reports[-1].total_energy_h - final.final_energy) < 1e-9
    assert np.allclose(reports[-1].atom("PD03").position, final.atom("PD03").position)
    # the cage opens during the run and the radial force changes sign
    h = reports[0].atom("H_01").position
    d_first = np.linalg.norm(reports[0].atom("PD03").position - h)
    d_last = np.linalg.norm(reports[-1].atom("PD03").position - h)
    assert d_first < 1.69 < 1.75 < d_last
    f_first = reports[0].atom("PD03").force @ (reports[0].atom("PD03").position - h)
    f_last = reports[-1].atom("PD03").force @ (reports[-1].atom("PD03").position - h)
    assert f_first > 0 > f_last


def test_reports_skip_frozen_atom_runs(calc_root):
    path = calc_root / "03_h_tetrahedral" / "stage1_wavefunctions" / "pd_super_tetra.prot.gz"
    assert parse_reports(path) == []


def test_committed_diagnostics(results_root):
    with open(results_root / "relaxation_diagnostics.json", encoding="utf-8") as fh:
        d = json.load(fh)
    k = d["cage_stiffness"]["k_radial_ev_per_ang2"]
    assert 20 < k < 45
    assert 1.70 < d["cage_stiffness"]["zero_force_distance_ang"] < 1.75
    runs = d["runs"]
    assert set(runs) >= {"octahedral", "tetrahedral", "g_0.0", "g_1.0"}
    assert runs["octahedral"]["max_force_free_atoms_mh_bohr"] < 2.0
    assert runs["tetrahedral"]["max_force_free_atoms_mh_bohr"] > 10.0
    for name, r in runs.items():
        assert abs(r["energy_drift_last_30_steps_mev"]) < 1.5, name
        assert r["ionic_temperature_end_k"] <= 5.0, name
        assert r["harmonic_upper_bound_unrelaxed_energy_mev"] < 60.0, name
    rows = d["generalized_force_check"]["rows"]
    ratios = [r["ratio"] for r in rows if r["ratio"] is not None and 0.15 < r["g"] < 0.55]
    assert all(0.6 < x < 0.9 for x in ratios)  # multiplier undershoots the slope on the O side
    ti = d["generalized_force_check"]["thermodynamic_integration_of_multiplier_mev"]
    direct = d["generalized_force_check"]["direct_profile_mev"]
    assert ti[-1] < 0.7 * direct[-1]  # integration of the multiplier does not recover E_T - E_O
