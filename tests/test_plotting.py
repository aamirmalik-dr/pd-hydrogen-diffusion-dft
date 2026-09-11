"""Smoke tests: every figure function renders to a file from small synthetic inputs."""

import numpy as np

from pdhdiff.cppaw_io import BandSegment
from pdhdiff.plotting import (
    hero_figure,
    plot_arrhenius,
    plot_bands,
    plot_convergence,
    plot_dos,
    plot_energy_profile,
    plot_path_geometry,
    plot_path_structures,
)
from pdhdiff.profile import fit_octahedral_well, fit_tetrahedral_well, locate_barrier
from pdhdiff.structure import supercell_pd_positions
from pdhdiff.tst import SiteModel, arrhenius_fit, diffusion_constant


def _profile():
    g = np.linspace(0, 1, 11)
    e = 0.2 * np.sin(np.pi * g) ** 2 + 0.08 * g**2
    return g, e


def test_profile_and_arrhenius_figures(tmp_path):
    g, e = _profile()
    bar = locate_barrier(g, e)
    wo, wt = fit_octahedral_well(g, e), fit_tetrahedral_well(g, e)
    fig = plot_energy_profile(g, e, bar, wo, wt, e_t_unconstrained_ev=e[-1])
    fig.savefig(tmp_path / "profile.png")
    model = SiteModel(e[-1], bar.e_ts_ev, bar.e_ts_ev - e[-1], 6e13, 1.5e14)
    t = np.linspace(200, 1000, 21)
    res = diffusion_constant(model, t)
    ea, d0 = arrhenius_fit(t, res.d_m2_s)
    fig = plot_arrhenius(
        res, ea, d0, alt={"boltzmann": diffusion_constant(model, t, None, "boltzmann")}
    )
    fig.savefig(tmp_path / "arrhenius.png")
    fig = hero_figure(
        dict(g=g, e_ev=e, barrier=bar, well_o=wo, well_t=wt),
        dict(result=res, ea_fit_ev=ea, d0_fit=d0),
    )
    fig.savefig(tmp_path / "hero.png")
    assert (tmp_path / "hero.png").stat().st_size > 10_000


def test_electronic_structure_figures(tmp_path):
    e = np.linspace(10, 30, 500)
    total = np.exp(-((e - 20) ** 2) / 2)
    curves = {k: total * f for k, f in (("s", 0.1), ("p", 0.1), ("d", 0.8))}
    fig = plot_dos(e, total, curves, e_fermi=21.0)
    fig.savefig(tmp_path / "dos.png")
    x = np.linspace(0, 1, 20)
    segs = [
        BandSegment(np.zeros(3), np.ones(3), x + i, np.outer(np.ones(20), [15.0, 18.0, 25.0]))
        for i in range(3)
    ]
    fig = plot_bands(segs, 21.0, ["A", "B", "C", "D"])
    fig.savefig(tmp_path / "bands.png")
    assert (tmp_path / "bands.png").exists()


def test_path_figures(tmp_path):
    a = 3.89
    pd = supercell_pd_positions(2) * a
    frames = [
        (pd, np.array([0.0, 0.0, 0.5 * a]) + x * np.array([0.25, 0.25, 0.25]) * a)
        for x in np.linspace(0, 1, 5)
    ]
    fig = plot_path_structures(frames, a)
    fig.savefig(tmp_path / "cage.png")
    g = np.linspace(0, 1, 5)
    fig = plot_path_geometry(g, g * 0.01, g * 0.1, g * 0.1)
    fig.savefig(tmp_path / "geom.png")
    tr = np.column_stack(
        [
            np.arange(50),
            np.zeros(50),
            np.zeros(50),
            np.zeros(50),
            -1 - np.exp(-np.arange(50) / 10),
            np.zeros(50),
        ]
    )
    fig = plot_convergence({"a": tr, "b": tr})
    fig.savefig(tmp_path / "conv.png")
    assert (tmp_path / "conv.png").exists()
