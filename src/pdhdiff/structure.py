"""Geometry of fcc palladium, its interstitial sites and the O to T reaction coordinate.

All positions in this module are in units of the cubic lattice constant unless a
function name or argument says otherwise. The conventions follow the CP-PAW
project: the octahedral start site is ``O = (0, 0, 1/2)``, the tetrahedral end
site is ``T = (1/4, 1/4, 3/4)``, and the three Pd atoms of the triangle the H
atom crosses are at ``(0, 0, 1)``, ``(1/2, 0, 1/2)`` and ``(0, 1/2, 1/2)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from pdhdiff.constants import A_PD_ANG, BOHR_ANG

FCC_BASIS = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
"""Pd positions in the cubic cell."""

OCTAHEDRAL_SITES = np.array([[0.5, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])
"""Octahedral interstitial sites of the cubic cell (one per Pd atom)."""

TETRAHEDRAL_SITES = np.array(list(product([0.25, 0.75], repeat=3)))
"""Tetrahedral interstitial sites of the cubic cell (two per Pd atom)."""

O_SITE = np.array([0.0, 0.0, 0.5])
"""Octahedral start site of the diffusion path."""

T_SITE = np.array([0.25, 0.25, 0.75])
"""Tetrahedral end site of the diffusion path."""

TRIANGLE_PD = np.array([[0.0, 0.0, 1.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
"""The three Pd atoms forming the face shared by the O and T cages."""

TRIANGLE_LABELS = ("PD17", "PD03", "PD04")
"""Labels of the triangle atoms in the 32-atom supercell used in the project."""

PATH_LENGTH = float(np.linalg.norm(T_SITE - O_SITE))
"""Distance from O to T in units of the lattice constant (sqrt(3)/4)."""


def path_length_ang(a: float = A_PD_ANG) -> float:
    """Length of the O to T path in angstrom."""
    return PATH_LENGTH * a


def supercell_pd_positions(n: int = 2) -> np.ndarray:
    """Pd positions of an ``n x n x n`` repetition of the cubic cell.

    The atom order matches the structure files of the project: the four basis
    atoms of each cube, cubes ordered with x fastest, then y, then z.

    Args:
        n: Number of cubic cells along each axis.

    Returns:
        Array of shape ``(4 n^3, 3)`` in units of the lattice constant.
    """
    pos = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                shift = np.array([i, j, k], dtype=float)
                for b in FCC_BASIS:
                    pos.append(b + shift)
    return np.array(pos)


def h_position_along_path(x: float) -> np.ndarray:
    """Ideal hydrogen position for reaction coordinate ``x`` between O (0) and T (1)."""
    return O_SITE + x * (T_SITE - O_SITE)


def reaction_coordinate(r_h: np.ndarray, r_pd_triangle: np.ndarray) -> float:
    """Translation-invariant linear reaction coordinate of the project.

    The coordinate is the projection of the H position, measured relative to the
    centre of the Pd triangle, onto the O to T direction, normalised so that it is
    0 at O and 1 at T for the unrelaxed crystal::

        g = (1, 1, 1) . (4/3 R_H - 4/9 sum_i R_Pd,i) + 2/3

    Args:
        r_h: Hydrogen position in units of the lattice constant.
        r_pd_triangle: Positions of the three triangle Pd atoms, shape ``(3, 3)``.

    Returns:
        The reaction coordinate ``g``.
    """
    ones = np.ones(3)
    return float(ones @ (4.0 / 3.0 * r_h - 4.0 / 9.0 * np.sum(r_pd_triangle, axis=0)) + 2.0 / 3.0)


def constraint_value_to_g(value_bohr: float, a: float = A_PD_ANG) -> float:
    """Convert the constraint value printed by CP-PAW (in bohr) to ``g``.

    CP-PAW evaluates the linear form with the coefficients of the ``!LINEAR`` block
    on Cartesian positions in bohr, which gives ``a (g - 2/3)``.
    """
    return value_bohr * BOHR_ANG / a + 2.0 / 3.0


def constraint_coefficients(h_label: str = "H_01", pd_labels=TRIANGLE_LABELS) -> dict:
    """Coefficient vectors of the ``!CONSTRAINTS!LINEAR`` block.

    The H coefficient is written as ``3 x 0.44444`` so that the coefficients sum to
    zero exactly after rounding, which keeps the constraint independent of a global
    translation.
    """
    c = 0.44444
    coeffs = {h_label: np.full(3, 3 * c)}
    for lab in pd_labels:
        coeffs[lab] = np.full(3, -c)
    return coeffs


def constraint_block(h_label: str = "H_01", pd_labels=TRIANGLE_LABELS, move: bool = False) -> str:
    """Text of the ``!CONSTRAINTS`` block for the structure file."""
    lines = ["\t!CONSTRAINTS", f"\t\t!LINEAR SHOW=T MOVE={'T' if move else 'F'}"]
    for name, vec in constraint_coefficients(h_label, pd_labels).items():
        v = " ".join(f"{x:.5f}" for x in vec)
        lines.append(f"\t\t\t!ATOM NAME='{name}' R= {v} !END")
    lines += ["\t\t!END", "\t!END"]
    return "\n".join(lines)


# --------------------------------------------------------------------------------------
# neighbour topology of the interstitial network
# --------------------------------------------------------------------------------------


def _images(points: np.ndarray, n: int = 1) -> np.ndarray:
    shifts = np.array(list(product(range(-n, n + 1), repeat=3)), dtype=float)
    return (points[:, None, :] + shifts[None, :, :]).reshape(-1, 3)


def interstitial_neighbours(site: np.ndarray, targets: np.ndarray, tol: float = 1e-6) -> np.ndarray:
    """Displacement vectors from ``site`` to the nearest periodic images of ``targets``.

    Only the shell at the O to T hop distance (``sqrt(3)/4`` lattice constants) is
    returned, which is the only hop the project considers.
    """
    imgs = _images(targets)
    d = imgs - site
    dist = np.linalg.norm(d, axis=1)
    mask = np.abs(dist - PATH_LENGTH) < tol
    return d[mask]


def geometric_factor(site: np.ndarray, targets: np.ndarray, a: float = A_PD_ANG) -> np.ndarray:
    """Sum of outer products ``(r_j - r_i) (x) (r_j - r_i)`` over hop targets, in angstrom^2.

    Args:
        site: Starting site in units of the lattice constant.
        targets: Candidate end sites in units of the lattice constant.
        a: Lattice constant in angstrom.

    Returns:
        A ``3 x 3`` tensor in angstrom^2.
    """
    d = interstitial_neighbours(site, targets) * a
    return d.T @ d


@dataclass
class HopGeometry:
    """Neighbour counts and geometric factors of the O and T sublattices.

    Attributes:
        n_t_around_o: Number of tetrahedral neighbours of an octahedral site.
        n_o_around_t: Number of octahedral neighbours of a tetrahedral site.
        factor_o_ang2: Geometric factor for hops out of an O site in angstrom^2.
        factor_t_ang2: Geometric factor for hops out of a T site in angstrom^2.
        hop_length_ang: O to T distance in angstrom.
    """

    n_t_around_o: int
    n_o_around_t: int
    factor_o_ang2: np.ndarray
    factor_t_ang2: np.ndarray
    hop_length_ang: float


def hop_geometry(a: float = A_PD_ANG) -> HopGeometry:
    """Evaluate the interstitial hop topology of fcc Pd numerically."""
    d_ot = interstitial_neighbours(O_SITE, TETRAHEDRAL_SITES)
    d_to = interstitial_neighbours(T_SITE, OCTAHEDRAL_SITES)
    return HopGeometry(
        n_t_around_o=len(d_ot),
        n_o_around_t=len(d_to),
        factor_o_ang2=geometric_factor(O_SITE, TETRAHEDRAL_SITES, a),
        factor_t_ang2=geometric_factor(T_SITE, OCTAHEDRAL_SITES, a),
        hop_length_ang=path_length_ang(a),
    )


# --------------------------------------------------------------------------------------
# writers
# --------------------------------------------------------------------------------------


def write_xyz(path, symbols, positions_ang, comment: str = "") -> None:
    """Write a single-frame xyz file.

    Args:
        path: Output path.
        symbols: Element symbols, one per atom.
        positions_ang: Cartesian positions in angstrom, shape ``(n_atoms, 3)``.
        comment: Second-line comment.
    """
    lines = [str(len(symbols)), comment]
    for s, p in zip(symbols, positions_ang):
        lines.append(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def write_extxyz_frames(path, frames: list[tuple[list[str], np.ndarray, str]]) -> None:
    """Write a multi-frame xyz trajectory (one block per frame)."""
    chunks = []
    for symbols, pos, comment in frames:
        lines = [str(len(symbols)), comment]
        for s, p in zip(symbols, pos):
            lines.append(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
        chunks.append("\n".join(lines))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(chunks) + "\n")
