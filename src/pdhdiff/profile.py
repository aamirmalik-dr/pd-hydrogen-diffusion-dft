"""Assemble and fit the energy profile of the O to T hop from the CP-PAW runs.

The profile is the relaxed total energy as a function of the linear reaction
coordinate ``g``. Every point of the path was produced by two CP-PAW runs: a
wave-function optimisation with frozen atoms (stage 1) and a damped atomic
relaxation with the constraint held (stage 2). Only stage 2 enters the profile.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from pdhdiff.constants import A_PD_ANG, BOHR_ANG, HARTREE_EV
from pdhdiff.cppaw_io import ProtocolSummary, parse_prot, parse_strc
from pdhdiff.structure import TRIANGLE_LABELS, constraint_value_to_g, reaction_coordinate

STAGE1 = "stage1_wavefunctions"
STAGE2 = "stage2_relaxation"


@dataclass
class PathPoint:
    """One relaxed point of the constrained path.

    Attributes:
        label: Directory name, for example ``g_0.5``.
        g_nominal: Reaction coordinate set by the initial hydrogen position.
        g_constraint: Reaction coordinate from the constraint value printed by CP-PAW.
        g_relaxed: Reaction coordinate recomputed from the relaxed positions.
        energy_h: Final total energy of the relaxation in Hartree.
        energy_stage1_h: Final energy of the frozen-atom wave-function run in Hartree.
        h_position_ang: Relaxed hydrogen position in angstrom.
        h_nominal_ang: Initial hydrogen position in angstrom.
        max_pd_displacement_ang: Largest Pd displacement from the ideal lattice.
        triangle_pd_displacement_ang: Mean displacement of the three triangle Pd atoms.
        constraint_force_h_per_bohr: Force on the constraint in Hartree per bohr.
        max_force_mh_per_bohr: Largest residual force component in milli-Hartree per bohr.
        n_iterations: Number of steps of the relaxation run.
        finished: Whether the relaxation protocol ends with a finished banner.
    """

    label: str
    g_nominal: float
    g_constraint: float
    g_relaxed: float
    energy_h: float
    energy_stage1_h: float | None
    h_position_ang: tuple[float, float, float]
    h_nominal_ang: tuple[float, float, float]
    max_pd_displacement_ang: float
    triangle_pd_displacement_ang: float
    constraint_force_h_per_bohr: float | None
    max_force_mh_per_bohr: float
    n_iterations: int | None
    finished: bool


def _ideal_pd_positions_ang(strc, a: float) -> np.ndarray:
    return strc.cartesian()[[i for i, s in enumerate(strc.species) if s.upper() == "PD"]]


def analyse_point(point_dir: Path, a: float = A_PD_ANG) -> PathPoint:
    """Extract the energy and geometry of one path point from its two-stage run.

    Args:
        point_dir: Directory holding ``stage1_wavefunctions`` and ``stage2_relaxation``.
        a: Lattice constant in angstrom.

    Returns:
        A :class:`PathPoint`.
    """
    stage2 = point_dir / STAGE2
    stage1 = point_dir / STAGE1
    prot_path = next(stage2.glob("*.prot*"))
    strc_path = next(stage2.glob("*.strc"))
    prot: ProtocolSummary = parse_prot(prot_path)
    strc = parse_strc(strc_path)
    e1 = None
    if stage1.exists():
        p1 = list(stage1.glob("*.prot*"))
        if p1:
            e1 = parse_prot(p1[0]).final_energy

    h_name = next(n for n, s in zip(strc.names, strc.species) if s.upper() == "H")
    h_nominal = strc.positions[strc.index(h_name)]
    g_nominal = reaction_coordinate(
        h_nominal, np.array([strc.positions[strc.index(n)] for n in TRIANGLE_LABELS])
    )

    h_relaxed = prot.atom(h_name).position
    tri = np.array([prot.atom(n).position for n in TRIANGLE_LABELS])
    g_relaxed = reaction_coordinate(h_relaxed / a, tri / a)

    pd_ideal = _ideal_pd_positions_ang(strc, a)
    pd_relaxed = np.array([at.position for at in prot.atoms if at.element == "Pd"])
    disp = pd_relaxed - pd_ideal
    disp -= np.round(disp / (2 * a)) * 2 * a  # minimum image in the 2x2x2 cell
    tri_ideal = np.array([strc.cartesian()[strc.index(n)] for n in TRIANGLE_LABELS])
    tri_disp = np.linalg.norm(tri - tri_ideal, axis=1).mean()

    forces = prot.forces()
    val = prot.final_constraint_value
    return PathPoint(
        label=point_dir.name,
        g_nominal=round(g_nominal, 6),
        g_constraint=constraint_value_to_g(val, a) if val is not None else float("nan"),
        g_relaxed=g_relaxed,
        energy_h=prot.final_energy,
        energy_stage1_h=e1,
        h_position_ang=tuple(float(v) for v in h_relaxed),
        h_nominal_ang=tuple(float(v) for v in h_nominal * a),
        max_pd_displacement_ang=float(np.linalg.norm(disp, axis=1).max()),
        triangle_pd_displacement_ang=float(tri_disp),
        constraint_force_h_per_bohr=(
            prot.constraint_values[-1][1] if prot.constraint_values else None
        ),
        max_force_mh_per_bohr=float(np.abs(forces).max()) if forces.size else float("nan"),
        n_iterations=prot.n_iterations,
        finished=prot.finished,
    )


def collect_path(path_root: Path, a: float = A_PD_ANG) -> list[PathPoint]:
    """Analyse every ``g_*`` directory below ``path_root``, sorted by ``g``."""
    points = [analyse_point(d, a) for d in sorted(path_root.glob("g_*")) if d.is_dir()]
    return sorted(points, key=lambda p: p.g_nominal)


def write_profile_csv(points: list[PathPoint], path: Path) -> None:
    """Write the path points to a CSV file, one row per point."""
    rows = []
    e0 = points[0].energy_h
    for p in points:
        d = asdict(p)
        d["energy_rel_mev"] = (p.energy_h - e0) * HARTREE_EV * 1e3
        d["h_position_ang"] = " ".join(f"{v:.5f}" for v in p.h_position_ang)
        d["h_nominal_ang"] = " ".join(f"{v:.5f}" for v in p.h_nominal_ang)
        rows.append(d)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_profile_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(g, E_rel_eV)`` from a profile CSV written by :func:`write_profile_csv`."""
    g, e = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            g.append(float(row["g_nominal"]))
            e.append(float(row["energy_rel_mev"]) / 1e3)
    return np.array(g), np.array(e)


# --------------------------------------------------------------------------------------
# fits
# --------------------------------------------------------------------------------------


@dataclass
class BarrierFit:
    """Location and height of the transition state.

    Attributes:
        g_ts: Reaction coordinate of the maximum.
        e_ts_ev: Energy of the maximum relative to the octahedral minimum.
        g_ts_spline: Maximum position from a cubic spline through all points.
        e_ts_spline_ev: Maximum energy from the spline.
        e_raw_max_ev: Highest calculated point.
        g_raw_max: Reaction coordinate of the highest calculated point.
    """

    g_ts: float
    e_ts_ev: float
    g_ts_spline: float
    e_ts_spline_ev: float
    e_raw_max_ev: float
    g_raw_max: float


def locate_barrier(g: np.ndarray, e_ev: np.ndarray, n_local: int = 4) -> BarrierFit:
    """Locate the transition state from the discrete profile.

    A cubic polynomial is fitted through the ``n_local`` points around the highest
    calculated point and maximised analytically on the fitted interval. A
    not-a-knot cubic spline through all points gives an independent estimate.

    Args:
        g: Reaction coordinates, increasing.
        e_ev: Energies relative to the octahedral minimum in eV.
        n_local: Number of points of the local polynomial fit.

    Returns:
        A :class:`BarrierFit`.
    """
    i_max = int(np.argmax(e_ev))
    lo = max(0, i_max - n_local // 2)
    hi = min(len(g), lo + n_local)
    lo = hi - n_local
    coeff = np.polyfit(g[lo:hi], e_ev[lo:hi], 3)
    grid = np.linspace(g[lo], g[hi - 1], 2001)
    vals = np.polyval(coeff, grid)
    j = int(np.argmax(vals))
    spline = CubicSpline(g, e_ev)
    grid_all = np.linspace(g[0], g[-1], 5001)
    sv = spline(grid_all)
    k = int(np.argmax(sv))
    return BarrierFit(
        g_ts=float(grid[j]),
        e_ts_ev=float(vals[j]),
        g_ts_spline=float(grid_all[k]),
        e_ts_spline_ev=float(sv[k]),
        e_raw_max_ev=float(e_ev[i_max]),
        g_raw_max=float(g[i_max]),
    )


@dataclass
class WellFit:
    """Harmonic curvature of a minimum of the profile along ``g``.

    Attributes:
        site: ``"O"`` or ``"T"``.
        k_g_ev: Force constant in eV per unit ``g^2``.
        anharmonic_coeff: Coefficient of the anharmonic term of the fit (quartic for O,
            cubic for T), in eV per unit ``g^4`` or ``g^3``.
        n_points: Number of points used.
        rms_residual_mev: Root-mean-square residual of the fit in meV.
        k_g_two_point_ev: Two-point estimate ``2 dE / dg^2`` from the nearest point, for
            comparison with the least-squares value.
    """

    site: str
    k_g_ev: float
    anharmonic_coeff: float
    n_points: int
    rms_residual_mev: float
    k_g_two_point_ev: float


def fit_octahedral_well(g: np.ndarray, e_ev: np.ndarray, g_max: float = 0.3) -> WellFit:
    """Fit ``E = E_O + k/2 g^2 + c g^4`` to the points with ``g <= g_max``.

    The octahedral site is a centre of inversion, so the profile is even in ``g``
    and the fit has no odd terms. ``E_O`` is fixed to the calculated energy at ``g = 0``.
    """
    mask = g <= g_max + 1e-9
    gg, ee = g[mask], e_ev[mask] - e_ev[0]
    design = np.column_stack([0.5 * gg**2, gg**4])
    coeff, *_ = np.linalg.lstsq(design, ee, rcond=None)
    resid = ee - design @ coeff
    two_point = 2.0 * (e_ev[1] - e_ev[0]) / g[1] ** 2
    return WellFit(
        site="O",
        k_g_ev=float(coeff[0]),
        anharmonic_coeff=float(coeff[1]),
        n_points=int(mask.sum()),
        rms_residual_mev=float(np.sqrt(np.mean(resid**2)) * 1e3),
        k_g_two_point_ev=float(two_point),
    )


def fit_tetrahedral_well(g: np.ndarray, e_ev: np.ndarray, g_min: float = 0.7) -> WellFit:
    """Fit ``E = E_T + k/2 d^2 + c d^3`` with ``d = g - 1`` to the points with ``g >= g_min``.

    The tetrahedral site has no inversion symmetry along the path (beyond ``g = 1``
    the H atom runs into a Pd atom), so a cubic term is required. ``E_T`` is fixed
    to the calculated energy at ``g = 1``.
    """
    mask = g >= g_min - 1e-9
    dd, ee = g[mask] - 1.0, e_ev[mask] - e_ev[-1]
    design = np.column_stack([0.5 * dd**2, dd**3])
    coeff, *_ = np.linalg.lstsq(design, ee, rcond=None)
    resid = ee - design @ coeff
    two_point = 2.0 * (e_ev[-2] - e_ev[-1]) / (g[-2] - 1.0) ** 2
    return WellFit(
        site="T",
        k_g_ev=float(coeff[0]),
        anharmonic_coeff=float(coeff[1]),
        n_points=int(mask.sum()),
        rms_residual_mev=float(np.sqrt(np.mean(resid**2)) * 1e3),
        k_g_two_point_ev=float(two_point),
    )


def spline_curve(g: np.ndarray, e_ev: np.ndarray, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """Dense cubic-spline interpolation of the profile for plotting."""
    spline = CubicSpline(g, e_ev)
    grid = np.linspace(g[0], g[-1], n)
    return grid, spline(grid)


def constraint_value_bohr(g: float, a: float = A_PD_ANG) -> float:
    """Inverse of :func:`pdhdiff.structure.constraint_value_to_g`."""
    return (g - 2.0 / 3.0) * a / BOHR_ANG
