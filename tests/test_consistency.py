"""Cross-checks between independently produced files of the calculations tree."""

import json
from pathlib import Path

import numpy as np
import pytest

from pdhdiff.constants import A_PD_ANG
from pdhdiff.cppaw_io import parse_prot, parse_strc
from pdhdiff.structure import TRIANGLE_LABELS, reaction_coordinate

ROOT = Path(__file__).resolve().parents[1]


def _read_xyz(path: Path):
    lines = [ln for ln in path.read_text().splitlines()[2:] if ln.strip()]
    symbols = [ln.split()[0] for ln in lines]
    pos = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines])
    return symbols, pos


@pytest.mark.parametrize(
    "run_dir",
    sorted((ROOT / "calculations" / "04_path_octa_to_tetra").glob("g_*"))
    + [ROOT / "calculations" / "03_h_octahedral", ROOT / "calculations" / "03_h_tetrahedral"],
    ids=lambda p: p.name,
)
def test_xyz_written_by_cppaw_matches_protocol_atom_list(run_dir):
    """The xyz files keep input order (H first); the protocol sorts atoms by species."""
    stage = run_dir / "stage2_relaxation"
    prot = parse_prot(next(stage.glob("*.prot.gz")))
    strc = parse_strc(next(stage.glob("*.strc")))
    symbols, pos = _read_xyz(next(stage.glob("*.xyz")))
    assert symbols[0] == "H" and symbols.count("PD") == 32
    cell = prot.lattice[0, 0]
    for name, p in zip(strc.names, pos):  # xyz rows follow the structure-file order
        d = p - prot.atom(name).position
        d -= np.round(d / cell) * cell
        assert np.linalg.norm(d) < 1e-4, name


def test_protocol_atom_order_is_species_sorted(path_root):
    prot = parse_prot(path_root / "g_0.5" / "stage2_relaxation" / "pd_super_octa.prot.gz")
    assert [a.name for a in prot.atoms][:2] == ["PD01", "PD02"]
    assert prot.atoms[-1].name == "H_01"


def test_energy_is_stationary_at_the_end_of_every_relaxation(path_root):
    """Convergence evidence used in the write-up: energy drift and ionic temperature."""
    for d in sorted(path_root.glob("g_*")):
        prot = parse_prot(d / "stage2_relaxation" / "pd_super_octa.prot.gz")
        tr = prot.last_trace
        drift_mev = abs(tr[-1, 4] - tr[-30, 4]) * 27.211386 * 1e3
        assert drift_mev < 1.5, d.name
        assert tr[-1, 2] <= 5.0, d.name  # ionic temperature in K
        assert prot.finished and prot.n_runs == 2


def test_input_generator_round_trip(tmp_path):
    import runpy
    import sys

    script = ROOT / "scripts" / "build_path_inputs.py"
    sys.argv = ["build_path_inputs.py", "--out", str(tmp_path), "--points", "3"]
    runpy.run_path(str(script), run_name="__main__")
    for x in (0.0, 0.5, 1.0):
        strc = parse_strc(tmp_path / f"g_{x:.1f}" / "stage2_relaxation" / "pd_super_octa.strc")
        tri = np.array([strc.positions[strc.index(n)] for n in TRIANGLE_LABELS])
        assert abs(reaction_coordinate(strc.positions[strc.index("H_01")], tri) - x) < 1e-9
        assert len(strc.linear_constraints) == 1 and strc.constraint_move is False
    ref = parse_strc(
        ROOT / "calculations/04_path_octa_to_tetra/g_0.5/stage2_relaxation/pd_super_octa.strc"
    )
    gen = parse_strc(tmp_path / "g_0.5" / "stage2_relaxation" / "pd_super_octa.strc")
    assert np.allclose(ref.positions, gen.positions)
    assert ref.linear_constraints[0].keys() == gen.linear_constraints[0].keys()


def test_metrics_sensitivity_block_brackets_the_headline_value(results_root):
    with open(results_root / "metrics.json", encoding="utf-8") as fh:
        m = json.load(fh)
    s = m["sensitivity"]
    d = m["diffusion"]["d_298K_m2_s"]
    assert s["d_298K_min_m2_s"] <= d <= s["d_298K_max_m2_s"]
    assert s["d_298K_max_m2_s"] / s["d_298K_min_m2_s"] < 3.0
    iso = m["isotopes_classical"]
    assert abs(iso["D"]["d_298K_over_h"] - 1 / np.sqrt(2.01410 / 1.00794)) < 1e-6


def test_lattice_constant_matches_protocol(path_root):
    prot = parse_prot(path_root / "g_0.0" / "stage2_relaxation" / "pd_super_octa.prot.gz")
    assert abs(prot.lattice[0, 0] - 2 * A_PD_ANG) < 1e-9
