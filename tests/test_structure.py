"""Geometry and reaction-coordinate tests."""

import numpy as np

from pdhdiff.constants import A_PD_ANG, BOHR_ANG
from pdhdiff.cppaw_io import parse_prot, parse_strc
from pdhdiff.structure import (
    O_SITE,
    T_SITE,
    TRIANGLE_LABELS,
    TRIANGLE_PD,
    constraint_block,
    constraint_value_to_g,
    h_position_along_path,
    hop_geometry,
    reaction_coordinate,
    supercell_pd_positions,
)


def test_reaction_coordinate_endpoints():
    assert abs(reaction_coordinate(O_SITE, TRIANGLE_PD)) < 1e-12
    assert abs(reaction_coordinate(T_SITE, TRIANGLE_PD) - 1.0) < 1e-12
    assert abs(reaction_coordinate(h_position_along_path(0.35), TRIANGLE_PD) - 0.35) < 1e-12


def test_reaction_coordinate_is_translation_invariant():
    shift = np.array([0.123, -0.4, 0.77])
    g0 = reaction_coordinate(h_position_along_path(0.6), TRIANGLE_PD)
    g1 = reaction_coordinate(h_position_along_path(0.6) + shift, TRIANGLE_PD + shift)
    assert abs(g0 - g1) < 1e-12


def test_supercell_matches_project_structure_file(calc_root):
    strc = parse_strc(calc_root / "02_pd_supercell" / "pd_super.strc")
    assert np.allclose(strc.positions, supercell_pd_positions(2))
    # the triangle atoms carry the labels used in the constraint block
    octa = parse_strc(calc_root / "03_h_octahedral" / "stage2_relaxation" / "pd_super_octa.strc")
    tri = np.array([octa.positions[octa.index(n)] for n in TRIANGLE_LABELS])
    assert np.allclose(tri, TRIANGLE_PD)


def test_cppaw_constraint_value_convention(path_root):
    """The value printed by CP-PAW equals a (g - 2/3) in bohr for every path point."""
    for d in sorted(path_root.glob("g_*")):
        strc = parse_strc(d / "stage2_relaxation" / "pd_super_octa.strc")
        prot = parse_prot(d / "stage2_relaxation" / "pd_super_octa.prot.gz")
        tri = np.array([strc.positions[strc.index(n)] for n in TRIANGLE_LABELS])
        g_in = reaction_coordinate(strc.positions[strc.index("H_01")], tri)
        g_out = constraint_value_to_g(prot.final_constraint_value, A_PD_ANG)
        assert abs(g_in - g_out) < 2e-3, d.name
        # and recomputing g from the relaxed Cartesian positions gives the same value
        tri_rel = np.array([prot.atom(n).position for n in TRIANGLE_LABELS]) / A_PD_ANG
        g_rel = reaction_coordinate(prot.atom("H_01").position / A_PD_ANG, tri_rel)
        assert abs(g_rel - g_in) < 1e-4, d.name


def test_constraint_block_text_matches_project_input(path_root):
    block = constraint_block()
    strc_text = (path_root / "g_0.0" / "stage2_relaxation" / "pd_super_octa.strc").read_text()
    normalised = " ".join(strc_text.split())
    for line in block.splitlines():
        if "!ATOM" in line:
            assert " ".join(line.split()) in normalised


def test_hop_geometry():
    geo = hop_geometry(A_PD_ANG)
    assert geo.n_t_around_o == 8
    assert geo.n_o_around_t == 4
    assert np.allclose(geo.factor_o_ang2, A_PD_ANG**2 / 2 * np.eye(3))
    assert np.allclose(geo.factor_t_ang2, A_PD_ANG**2 / 4 * np.eye(3))
    assert abs(geo.hop_length_ang - np.sqrt(3) / 4 * A_PD_ANG) < 1e-12


def test_bohr_constant():
    assert abs(A_PD_ANG / BOHR_ANG - 7.35103) < 1e-4
