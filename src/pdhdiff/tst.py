"""Transition-state theory in the one-dimensional harmonic approximation.

The jump rate of a thermally activated hop is::

    Gamma = (omega_0 / 2 pi) exp(-E_a / k_B T)

with the attempt frequency ``omega_0`` from the curvature of the energy profile at
the initial site and the activation energy ``E_a`` from the barrier height. The
diffusion constant follows from the master equation on the O and T sublattices,

    D = 1/2 sum_i P_i sum_j Gamma_(j<-i) (r_j - r_i) (x) (r_j - r_i),

where ``P_i`` are equilibrium site populations and the geometric sums run over
the hop targets of a site. The factor 1/2 is the Einstein relation
``<dr (x) dr> = 2 D t``; the project description prints the expression without
it and the code can reproduce that convention with ``einstein_half=False``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pdhdiff.constants import AMU_KG, EV_J, HBAR_EV_S, KB_EV, MASS_H_AMU
from pdhdiff.structure import HopGeometry, hop_geometry


def force_constant_si(k_g_ev: float, path_length_ang: float) -> float:
    """Convert a force constant in eV per ``g^2`` to N/m.

    Args:
        k_g_ev: Curvature of the profile in eV per unit reaction coordinate squared.
        path_length_ang: Length of the path (``g`` from 0 to 1) in angstrom.

    Returns:
        Force constant in N/m.
    """
    k_ev_ang2 = k_g_ev / path_length_ang**2
    return k_ev_ang2 * EV_J / 1e-20


def attempt_frequency(k_si: float, mass_amu: float = MASS_H_AMU) -> float:
    """Angular attempt frequency ``omega_0 = sqrt(k / m)`` in rad/s."""
    return float(np.sqrt(k_si / (mass_amu * AMU_KG)))


def vibrational_quantum_mev(omega: float) -> float:
    """``hbar omega`` in meV for an angular frequency in rad/s."""
    return HBAR_EV_S * omega * 1e3


def jump_rate(omega: float, ea_ev: float, temperature_k) -> np.ndarray:
    """Harmonic TST rate per pathway, ``(omega / 2 pi) exp(-E_a / k_B T)``, in 1/s."""
    t = np.asarray(temperature_k, dtype=float)
    return omega / (2 * np.pi) * np.exp(-ea_ev / (KB_EV * t))


@dataclass
class SiteModel:
    """The ingredients of the two-site diffusion model.

    Attributes:
        de_t_minus_o_ev: Energy of the tetrahedral minimum above the octahedral one.
        ea_o_to_t_ev: Barrier seen from the octahedral site.
        ea_t_to_o_ev: Barrier seen from the tetrahedral site.
        omega_o: Attempt frequency at O in rad/s.
        omega_t: Attempt frequency at T in rad/s.
        mult_o: Number of O sites per Pd atom (1).
        mult_t: Number of T sites per Pd atom (2).
    """

    de_t_minus_o_ev: float
    ea_o_to_t_ev: float
    ea_t_to_o_ev: float
    omega_o: float
    omega_t: float
    mult_o: int = 1
    mult_t: int = 2


def site_populations(
    model: SiteModel, temperature_k, mode: str = "harmonic1d"
) -> tuple[np.ndarray, np.ndarray]:
    """Fraction of time the hydrogen spends on the O and on the T sublattice.

    Args:
        model: Site energies and frequencies.
        temperature_k: Temperature or array of temperatures in kelvin.
        mode: ``"harmonic1d"`` (default in the analysis) uses the prefactor
            ``omega_O / omega_T`` of a one-dimensional harmonic well, which satisfies
            detailed balance together with the one-dimensional TST rates.
            ``"harmonic3d"`` applies the prefactor ``(omega_O / omega_T)^3`` of the
            project description, assuming isotropic wells with the path curvature.
            ``"boltzmann"`` uses the site energies only.

    Returns:
        ``(P_oct, P_tet)`` with ``P_oct + P_tet = 1``.
    """
    t = np.asarray(temperature_k, dtype=float)
    boltz = np.exp(-model.de_t_minus_o_ev / (KB_EV * t))
    if mode == "harmonic3d":
        pref = (model.omega_o / model.omega_t) ** 3
    elif mode == "harmonic1d":
        pref = model.omega_o / model.omega_t
    elif mode == "boltzmann":
        pref = 1.0
    else:
        raise ValueError(mode)
    w_o = model.mult_o * 1.0
    w_t = model.mult_t * pref * boltz
    return w_o / (w_o + w_t), w_t / (w_o + w_t)


@dataclass
class DiffusionResult:
    """Diffusion constant and its ingredients on a temperature grid.

    Attributes:
        temperature_k: Temperatures.
        gamma_t_from_o: Rate of a single O to T hop in 1/s.
        gamma_o_from_t: Rate of a single T to O hop in 1/s.
        p_oct: Population of the octahedral sublattice.
        p_tet: Population of the tetrahedral sublattice.
        d_m2_s: Diffusion constant (scalar, cubic symmetry) in m^2/s.
        term_o_m2_s: Contribution of hops out of O sites.
        term_t_m2_s: Contribution of hops out of T sites.
        detailed_balance_ratio: ``(P_oct n_TO Gamma_TO) / (P_tet n_OT Gamma_OT)``, the ratio
            of the forward and backward fluxes between the sublattices. It equals 1 when
            populations and rates are mutually consistent, which is the case for the
            ``"harmonic1d"`` populations combined with the one-dimensional rates.
    """

    temperature_k: np.ndarray
    gamma_t_from_o: np.ndarray
    gamma_o_from_t: np.ndarray
    p_oct: np.ndarray
    p_tet: np.ndarray
    d_m2_s: np.ndarray
    term_o_m2_s: np.ndarray
    term_t_m2_s: np.ndarray
    detailed_balance_ratio: np.ndarray


def diffusion_constant(
    model: SiteModel,
    temperature_k,
    geometry: HopGeometry | None = None,
    population_mode: str = "harmonic1d",
    einstein_half: bool = True,
) -> DiffusionResult:
    """Evaluate the diffusion constant of the two-sublattice hop model.

    Args:
        model: Site energies, barriers and attempt frequencies.
        temperature_k: Temperature grid in kelvin.
        geometry: Hop topology; computed for the fcc lattice constant when ``None``.
        population_mode: See :func:`site_populations`.
        einstein_half: Include the factor 1/2 of the Einstein relation. ``False``
            reproduces the expression as printed in the project description.

    Returns:
        A :class:`DiffusionResult` in SI units.
    """
    geo = geometry or hop_geometry()
    t = np.atleast_1d(np.asarray(temperature_k, dtype=float))
    g_to = jump_rate(model.omega_o, model.ea_o_to_t_ev, t)
    g_ot = jump_rate(model.omega_t, model.ea_t_to_o_ev, t)
    p_o, p_t = site_populations(model, t, population_mode)
    # the tensors are isotropic; take the trace / 3 and convert angstrom^2 to m^2
    f_o = np.trace(geo.factor_o_ang2) / 3 * 1e-20
    f_t = np.trace(geo.factor_t_ang2) / 3 * 1e-20
    pref = 0.5 if einstein_half else 1.0
    term_o = pref * p_o * g_to * f_o
    term_t = pref * p_t * g_ot * f_t
    # detailed balance between the sublattices: P_oct n_TO Gamma_TO = P_tet n_OT Gamma_OT
    ratio = (p_o * geo.n_t_around_o * g_to) / (p_t * geo.n_o_around_t * g_ot)
    return DiffusionResult(
        temperature_k=t,
        gamma_t_from_o=g_to,
        gamma_o_from_t=g_ot,
        p_oct=p_o,
        p_tet=p_t,
        d_m2_s=term_o + term_t,
        term_o_m2_s=term_o,
        term_t_m2_s=term_t,
        detailed_balance_ratio=ratio,
    )


def with_mass(
    model: SiteModel, mass_amu: float, reference_mass_amu: float = MASS_H_AMU
) -> SiteModel:
    """Classical isotope substitution: scale both attempt frequencies by ``sqrt(m_ref / m)``.

    The energy profile is mass independent in the Born-Oppenheimer approximation,
    so in classical harmonic TST only the prefactors change. Quantum effects
    (zero-point energy, tunnelling), which dominate the real isotope effect of
    hydrogen, are not included.

    Args:
        model: Site model of the reference isotope.
        mass_amu: Mass of the substituted isotope in atomic mass units.
        reference_mass_amu: Mass used for ``model``.

    Returns:
        A new :class:`SiteModel` with scaled frequencies.
    """
    scale = float(np.sqrt(reference_mass_amu / mass_amu))
    return SiteModel(
        de_t_minus_o_ev=model.de_t_minus_o_ev,
        ea_o_to_t_ev=model.ea_o_to_t_ev,
        ea_t_to_o_ev=model.ea_t_to_o_ev,
        omega_o=model.omega_o * scale,
        omega_t=model.omega_t * scale,
        mult_o=model.mult_o,
        mult_t=model.mult_t,
    )


def arrhenius_fit(temperature_k: np.ndarray, d_m2_s: np.ndarray) -> tuple[float, float]:
    """Fit ``ln D = ln D0 - E_a / (k_B T)`` and return ``(E_a in eV, D0 in m^2/s)``."""
    x = 1.0 / (KB_EV * np.asarray(temperature_k, dtype=float))
    y = np.log(np.asarray(d_m2_s, dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    return float(-slope), float(np.exp(intercept))
