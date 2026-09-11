"""Figures for the H in Pd diffusion study.

Every function takes data arrays and returns a Matplotlib figure. Colours are
assigned by role: blue for the primary series, orange for the second, aqua for
the third, with text always in neutral ink. The palette was validated for
colour-vision deficiency with the adjacent-pair test.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pdhdiff.constants import D_EXP_298K_M2_S, KB_EV  # noqa: E402
from pdhdiff.profile import BarrierFit, WellFit, spline_curve  # noqa: E402
from pdhdiff.tst import DiffusionResult  # noqa: E402

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e1"
FILL_BLUE = "#cde2fb"

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.facecolor": "white",
}


def _style() -> None:
    plt.rcParams.update(STYLE)


# --------------------------------------------------------------------------------------
# energy profile
# --------------------------------------------------------------------------------------


def plot_energy_profile(
    g: np.ndarray,
    e_ev: np.ndarray,
    barrier: BarrierFit,
    well_o: WellFit,
    well_t: WellFit,
    e_t_unconstrained_ev: float | None = None,
    ax=None,
):
    """Relaxed energy along the reaction coordinate with the fits used for TST.

    Args:
        g: Reaction coordinates.
        e_ev: Energies relative to the octahedral minimum in eV.
        barrier: Transition-state fit.
        well_o: Harmonic fit of the octahedral well.
        well_t: Harmonic fit of the tetrahedral well.
        e_t_unconstrained_ev: Energy of the unconstrained tetrahedral relaxation
            relative to the unconstrained octahedral one, drawn as a reference marker.
        ax: Existing axes, or ``None`` to create a figure.

    Returns:
        The Matplotlib figure.
    """
    _style()
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
    else:
        fig = ax.figure
    mev = 1e3
    xs, ys = spline_curve(g, e_ev)
    ax.fill_between(xs, 0, ys * mev, color=FILL_BLUE, alpha=0.45, linewidth=0)
    ax.plot(xs, ys * mev, color=BLUE, lw=2, label="cubic spline through the 11 points")
    ax.plot(
        g,
        e_ev * mev,
        "o",
        color=BLUE,
        ms=6,
        mec="white",
        mew=1.2,
        label="constrained relaxations (PBE)",
    )

    go = np.linspace(0, 0.42, 100)
    ax.plot(
        go,
        (0.5 * well_o.k_g_ev * go**2 + well_o.anharmonic_coeff * go**4) * mev,
        color=AQUA,
        lw=2,
        ls="--",
        label=f"O well fit, k = {well_o.k_g_ev:.2f} eV/g$^2$",
    )
    gt = np.linspace(0.62, 1.05, 100)
    d = gt - 1
    ax.plot(
        gt,
        (e_ev[-1] + 0.5 * well_t.k_g_ev * d**2 + well_t.anharmonic_coeff * d**3) * mev,
        color=ORANGE,
        lw=2,
        ls="--",
        label=f"T well fit, k = {well_t.k_g_ev:.1f} eV/g$^2$",
    )
    ax.plot(barrier.g_ts, barrier.e_ts_ev * mev, marker="x", color=INK, ms=9, mew=2, ls="none")
    ax.annotate(
        f"transition state\ng = {barrier.g_ts:.2f}, {barrier.e_ts_ev * mev:.0f} meV",
        (barrier.g_ts, barrier.e_ts_ev * mev),
        xytext=(barrier.g_ts + 0.08, barrier.e_ts_ev * mev + 28),
        color=INK,
        fontsize=9,
        ha="left",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8),
    )
    ax.annotate(
        f"T site\n+{e_ev[-1] * mev:.0f} meV",
        (1.0, e_ev[-1] * mev),
        xytext=(0.90, e_ev[-1] * mev - 40),
        color=INK2,
        fontsize=9,
        ha="center",
        va="top",
    )
    ax.annotate("O site, 0 meV", (0.0, 0.0), xytext=(0.03, -14), color=INK2, fontsize=9, va="top")
    if e_t_unconstrained_ev is not None:
        ax.plot(
            1.0,
            e_t_unconstrained_ev * mev,
            marker="s",
            color=ORANGE,
            ms=7,
            mfc="none",
            mew=1.5,
            ls="none",
            label="unconstrained T relaxation",
        )
    ax.axvline(2 / 3, color=INK2, lw=0.8, ls=":")
    ax.text(2 / 3 + 0.01, 8, "Pd triangle\n(g = 2/3)", color=INK2, fontsize=8, va="bottom")
    ax.set_xlabel("reaction coordinate g (0 = octahedral, 1 = tetrahedral)")
    ax.set_ylabel("E $-$ E$_\\mathrm{O}$ (meV)")
    ax.set_xlim(-0.03, 1.07)
    ax.set_ylim(-32, max(e_ev) * mev + 75)
    ax.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 1.0))
    ax.set_title(
        "H hop in fcc Pd: relaxed energy profile, 32-atom supercell", loc="left", fontsize=11
    )
    return fig


# --------------------------------------------------------------------------------------
# Arrhenius plot
# --------------------------------------------------------------------------------------


def plot_arrhenius(
    result: DiffusionResult,
    ea_fit_ev: float,
    d0_fit: float,
    alt: dict[str, DiffusionResult] | None = None,
    ax=None,
):
    """``ln D`` against ``1 / k_B T`` with the experimental room-temperature value.

    Args:
        result: Primary diffusion result.
        ea_fit_ev: Fitted activation energy of the primary curve.
        d0_fit: Fitted prefactor of the primary curve.
        alt: Optional named alternative results drawn as thin lines.
        ax: Existing axes, or ``None``.

    Returns:
        The Matplotlib figure.
    """
    _style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
    else:
        fig = ax.figure
    x = 1 / (KB_EV * result.temperature_k)
    ax.plot(x, np.log(result.d_m2_s), color=BLUE, lw=2.2, label="this work, harmonic TST")
    if alt:
        styles = [(ORANGE, "--"), (AQUA, ":"), (VIOLET, "-.")]
        for (name, res), (col, ls) in zip(alt.items(), styles):
            ax.plot(
                1 / (KB_EV * res.temperature_k),
                np.log(res.d_m2_s),
                color=col,
                lw=1.4,
                ls=ls,
                label=name,
            )
    x298 = 1 / (KB_EV * 298.0)
    ax.plot(
        x298,
        np.log(D_EXP_298K_M2_S),
        marker="D",
        color=INK,
        ms=7,
        ls="none",
        label="experiment, 298 K (Powell 1991)",
    )
    d298 = float(np.interp(298.0, result.temperature_k, result.d_m2_s))
    ax.text(
        0.03,
        0.06,
        f"this work: D(298 K) = {d298:.1e} m$^2$/s, E$_a$ = {ea_fit_ev:.3f} eV, "
        f"D$_0$ = {d0_fit:.1e} m$^2$/s\nexperiment: D(298 K) = {D_EXP_298K_M2_S:.1e} m$^2$/s",
        transform=ax.transAxes,
        fontsize=8.5,
        color=INK,
        ha="left",
        va="bottom",
    )
    ax.set_xlabel("1 / (k$_B$T)  (eV$^{-1}$)")
    ax.set_ylabel("ln D  (D in m$^2$/s)")
    t_ticks = np.array([1000, 600, 400, 300, 200])
    sec = ax.secondary_xaxis("top", functions=(lambda v: v, lambda v: v))
    sec.set_xticks(1 / (KB_EV * t_ticks))
    sec.set_xticklabels([f"{t:d} K" for t in t_ticks])
    sec.tick_params(colors=INK2)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Arrhenius plot of the diffusion constant", loc="left", fontsize=11)
    return fig


# --------------------------------------------------------------------------------------
# electronic structure
# --------------------------------------------------------------------------------------


def plot_dos(
    energy_ev,
    total: np.ndarray,
    curves: dict[str, np.ndarray],
    e_fermi: float,
    spin_factor: float = 2.0,
    ax=None,
):
    """Total and projected density of states of bulk Pd with E_F at zero.

    The projected s, p and d weights are not additive with the total (they are
    projections onto the augmentation sphere), so they are drawn as lines over
    the filled total.

    Args:
        energy_ev: Absolute energy grid in eV.
        total: Total DOS per spin.
        curves: Mapping ``"s"``, ``"p"``, ``"d"`` to the projected DOS per spin.
        e_fermi: Chemical potential in eV.
        spin_factor: Multiplier for both spins of a non-spin-polarised run.
        ax: Existing axes, or ``None``.
    """
    _style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 3.6))
    else:
        fig = ax.figure
    e = energy_ev - e_fermi
    occ = e <= 0
    ax.fill_between(e, 0, total * spin_factor, color="#d9d8d2", lw=0, label="total")
    ax.fill_between(
        e[occ], 0, (total * spin_factor)[occ], color="#b9b8b0", lw=0, label="total, occupied"
    )
    colors = {"d": BLUE, "p": ORANGE, "s": AQUA}
    for key in ("d", "p", "s"):
        if key in curves:
            ax.plot(
                e,
                curves[key] * spin_factor,
                color=colors[key],
                lw=1.6,
                label=f"Pd {key} projection",
            )
    ax.axvline(0, color=INK, lw=1.0, ls="--")
    ymax = float((total * spin_factor).max()) * 1.08
    ax.text(0.35, ymax * 0.93, "E$_F$", color=INK, fontsize=9)
    ax.set_xlim(-8.5, 8.0)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("E $-$ E$_F$ (eV)")
    ax.set_ylabel("DOS (states / eV / atom, both spins)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Bulk Pd, primitive cell, PBE: density of states", loc="left", fontsize=11)
    return fig


def plot_bands(segments, e_fermi: float, labels: list[str], ax=None):
    """Band structure along a k-path with ``E_F`` at zero.

    Args:
        segments: Output of :func:`pdhdiff.cppaw_io.parse_bands`.
        e_fermi: Chemical potential in eV.
        labels: Names of the segment boundaries (``len(segments) + 1`` entries).
        ax: Existing axes, or ``None``.
    """
    _style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.2, 4.6))
    else:
        fig = ax.figure
    ticks = [segments[0].x[0]]
    for seg in segments:
        ax.plot(seg.x, seg.energies_ev - e_fermi, color=BLUE, lw=1.3)
        ticks.append(seg.x[-1])
    for t in ticks[1:-1]:
        ax.axvline(t, color=INK2, lw=0.6)
    ax.axhline(0, color=INK, lw=1.0, ls="--")
    ax.axhspan(-9, 0, color=GRID, alpha=0.35, lw=0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_ylim(-8.5, 6)
    ax.set_ylabel("E $-$ E$_F$ (eV)")
    ax.grid(False, axis="x")
    ax.set_title("Bulk Pd, PBE: band structure", loc="left", fontsize=11)
    return fig


# --------------------------------------------------------------------------------------
# convergence and geometry along the path
# --------------------------------------------------------------------------------------


def plot_convergence(traces: dict[str, np.ndarray], ax=None):
    """Energy against step for a few relaxations, relative to each final energy.

    Args:
        traces: Mapping ``label -> array`` with the ``!>`` trace columns
            (``nfi, t_ps, T_K, ekin_psi, e_rho, econs``).
    """
    _style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
    else:
        fig = ax.figure
    cols = [BLUE, ORANGE, AQUA, VIOLET]
    for (label, tr), col in zip(traces.items(), cols):
        de = np.abs(tr[:, 4] - tr[-1, 4]) * 27.211386 * 1e3
        ax.plot(tr[:, 0] - tr[0, 0] + 1, np.maximum(de, 1e-2), color=col, lw=1.6, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("time step of the relaxation run")
    ax.set_ylabel("|E(step) $-$ E(final)|  (meV)")
    ax.legend(fontsize=8)
    ax.set_title("Damped atomic relaxation with the constraint held", loc="left", fontsize=11)
    return fig


def plot_path_geometry(g, h_offsets_ang, max_pd_disp_ang, tri_disp_ang, ax=None):
    """Relaxed H offset from the straight O-T line and Pd displacements along the path."""
    _style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
    else:
        fig = ax.figure
    ax.plot(g, max_pd_disp_ang, "o-", color=BLUE, lw=1.8, ms=5, label="largest Pd displacement")
    ax.plot(
        g,
        tri_disp_ang,
        "s-",
        color=ORANGE,
        lw=1.8,
        ms=5,
        label="mean displacement of the 3 triangle Pd",
    )
    ax.plot(
        g,
        h_offsets_ang,
        "^-",
        color=AQUA,
        lw=1.8,
        ms=5,
        label="H offset from the straight O$\\to$T line",
    )
    ax.set_xlabel("reaction coordinate g")
    ax.set_ylabel("displacement (Å)")
    ax.legend(fontsize=8)
    ax.set_title("Lattice response along the path", loc="left", fontsize=11)
    return fig


def _wrap(points: np.ndarray, centre: np.ndarray, cell: float) -> np.ndarray:
    """Minimum-image positions relative to ``centre`` in a cubic cell of edge ``cell``."""
    d = points - centre
    return d - np.round(d / cell) * cell


def plot_path_structures(
    frames: list[tuple[np.ndarray, np.ndarray]], a: float, cell: float | None = None
):
    """Three-dimensional view of the relaxed H positions inside the O and T cages.

    The six Pd atoms of the octahedral cage and the four of the tetrahedral cage
    (three are shared) are taken from the relaxed structure at the tetrahedral
    end point. The shared triangle the H crosses is drawn as a filled face.

    Args:
        frames: List of ``(pd_positions, h_position)`` in angstrom for every path
            point, relaxed, in path order.
        a: Lattice constant in angstrom.
        cell: Supercell edge in angstrom (``2 a`` when ``None``).

    Returns:
        The Matplotlib figure.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: PLC0415

    _style()
    cell = cell or 2 * a
    fig = plt.figure(figsize=(6.4, 5.6))
    ax = fig.add_subplot(projection="3d")
    ax.grid(False)
    o_site = np.array([0.0, 0.0, 0.5]) * a
    t_site = np.array([0.25, 0.25, 0.75]) * a
    pd_final = frames[-1][0]
    rel = _wrap(pd_final, o_site, cell)
    d_o = np.linalg.norm(rel, axis=1)
    d_t = np.linalg.norm(rel + o_site - t_site, axis=1)
    cage_o = rel[d_o < 0.55 * a] + o_site
    cage_t = rel[d_t < 0.50 * a] + o_site
    cage = np.unique(np.vstack([cage_o, cage_t]).round(4), axis=0)
    ax.scatter(
        cage[:, 0],
        cage[:, 1],
        cage[:, 2],
        s=900,
        color="#d9d8d2",
        edgecolor=INK2,
        alpha=0.75,
        depthshade=False,
    )

    def edges(points, dmax):
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                if np.linalg.norm(points[i] - points[j]) < dmax:
                    ax.plot(*zip(points[i], points[j]), color=INK2, lw=0.7, alpha=0.6)

    edges(cage_o, 0.75 * a)
    edges(cage_t, 0.75 * a)
    shared = np.array([p for p in cage_o if any(np.linalg.norm(p - q) < 1e-3 for q in cage_t)])
    if len(shared) == 3:
        ax.add_collection3d(
            Poly3DCollection([shared], facecolor=ORANGE, alpha=0.25, edgecolor=ORANGE)
        )
    n = len(frames)
    hs = np.array([_wrap(h[None, :], o_site, cell)[0] + o_site for _, h in frames])
    ax.plot(hs[:, 0], hs[:, 1], hs[:, 2], color=BLUE, lw=1.6)
    for i, h in enumerate(hs):
        shade = plt.cm.Blues(0.35 + 0.6 * i / max(n - 1, 1))
        ax.scatter(*h, s=140, color=shade, edgecolor="white", depthshade=False)
    ax.text(*(hs[0] + [-0.55, -0.25, -0.45]), "O, g = 0", color=INK, fontsize=9)
    ax.text(*(hs[-1] + [0.15, 0.15, 0.25]), "T, g = 1", color=INK, fontsize=9)
    ax.text(
        *(shared.mean(axis=0) + [0.9, -0.2, -0.75]), "shared Pd triangle", color=ORANGE, fontsize=8
    )
    lo, hi = cage.min(axis=0) - 0.6, cage.max(axis=0) + 0.6
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(hi - lo)
    ax.view_init(elev=22, azim=-118)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title(
        "Relaxed H positions (light to dark, g = 0 to 1) in the Pd cages", loc="left", fontsize=11
    )
    return fig


def hero_figure(profile_kwargs: dict, arrhenius_kwargs: dict):
    """Side-by-side profile and Arrhenius plot used as the README hero."""
    _style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.4))
    plot_energy_profile(ax=ax1, **profile_kwargs)
    plot_arrhenius(ax=ax2, **arrhenius_kwargs)
    fig.tight_layout(w_pad=2.5)
    return fig
