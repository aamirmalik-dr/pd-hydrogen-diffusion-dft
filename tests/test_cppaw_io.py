"""Parser tests against the committed CP-PAW files."""

import numpy as np

from pdhdiff.cppaw_io import parse_bands, parse_blocks, parse_dos, parse_prot, parse_strc


def test_prot_primitive_cell(calc_root):
    prot = parse_prot(calc_root / "01_pd_primitive" / "pd.prot.gz")
    assert prot.n_kpoints == 63
    assert prot.n_bands == 10
    assert prot.cutoff_ry == 30.0
    assert prot.finished
    assert abs(prot.final_energy - (-30.0020008)) < 1e-6
    assert abs(prot.fermi_energy_ev - 21.43526) < 1e-4
    assert len(prot.atoms) == 1 and prot.atoms[0].element == "Pd"
    assert np.allclose(prot.lattice, 1.945 * (np.ones((3, 3)) - np.eye(3)))
    assert prot.trace.shape[1] == 6 and len(prot.trace) > 100


def test_prot_relaxation_has_forces_and_constraint(path_root):
    prot = parse_prot(path_root / "g_0.5" / "stage2_relaxation" / "pd_super_octa.prot.gz")
    assert prot.n_kpoints == 8
    assert len(prot.atoms) == 33
    assert prot.n_iterations == 364  # from the run-time report, not the NSTEP setting
    h = prot.atom("H_01")
    assert h.element == "H" and h.force is not None
    # residual forces: the constrained atoms carry the constraint force (order 10 mH/bohr),
    # all other atoms are relaxed below 3 mH/bohr
    free = [a for a in prot.atoms if a.name not in ("H_01", "PD17", "PD03", "PD04")]
    assert np.abs(np.array([a.force for a in free])).max() < 3.0
    assert np.abs(prot.forces()).max() < 30.0
    val, force = prot.constraint_values[-1]
    assert abs(val - (-1.23)) < 0.01
    assert abs(force) < 5e-3


def test_strc_supercell_with_constraint(path_root):
    strc = parse_strc(path_root / "g_0.3" / "stage2_relaxation" / "pd_super_octa.strc")
    assert strc.lunit_ang == 3.89
    assert np.allclose(strc.lattice, 2 * np.eye(3))
    assert len(strc.names) == 33 and strc.species.count("PD") == 32
    assert strc.kpoint_density == 20.0 and strc.empty_bands == 5
    assert np.allclose(strc.positions[strc.index("H_01")], [0.075, 0.075, 0.575])
    assert len(strc.linear_constraints) == 1
    coeff = strc.linear_constraints[0]
    assert set(coeff) == {"H_01", "PD17", "PD03", "PD04"}
    assert abs(sum(v[0] for v in coeff.values())) < 1e-6  # translation invariant
    assert strc.constraint_move is False


def test_strc_out_bohr_units(tmp_path):
    text = """!STRUCTURE
!GENERIC
LUNIT=
   7.3510341111090964
!END
!LATTICE
T=
 0.0 0.5 0.5
 0.5 0.0 0.5
 0.5 0.5 0.0
!END
!ATOM
NAME=
 'PD1'
R=
 0.0 0.0 0.0
SP=
 'PD'
!END
!END
!EOB
"""
    p = tmp_path / "x.strc_out"
    p.write_text(text)
    strc = parse_strc(p)
    assert abs(strc.lunit_ang - 3.89) < 1e-4
    assert strc.species == ["PD"]


def test_block_parser_handles_fortran_exponents():
    root = parse_blocks("!A X=1.D-6 Y='q' Z=T !B W=2 !END !END")
    a = root.first("A")
    assert a.scalar("X") == 1e-6 and a.scalar("Y") == "q" and a.scalar("Z") == "T"
    assert a.first("B").scalar("W") == 2.0


def test_dos_is_positive_monotonic_and_integrates_to_valence(calc_root):
    dos = parse_dos(calc_root / "01_pd_primitive" / "dos" / "total.dos.gz")
    assert np.all(np.diff(dos.energy_ev) > 0)
    assert dos.occupied.min() >= 0 and dos.unoccupied.min() >= 0
    n_occ = 2 * np.trapezoid(dos.occupied, dos.energy_ev)  # both spins
    assert abs(n_occ - 10.0) < 0.1  # Pd has 10 valence electrons in this setup


def test_bands_segments(calc_root):
    segs = parse_bands(calc_root / "01_pd_primitive" / "pd_bands.dat")
    assert len(segs) == 6
    assert segs[0].energies_ev.shape == (40, 8)
    x = np.concatenate([s.x for s in segs])
    assert np.all(np.diff(x) >= 0)
    # X point in reduced coordinates of the primitive reciprocal lattice
    assert np.allclose(segs[0].k_start, 0) and np.allclose(segs[0].k_end, [0, 0.5, 0.5], atol=1e-4)
