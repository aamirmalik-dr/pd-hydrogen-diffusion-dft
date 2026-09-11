"""Physical constants and the fixed parameters of the Pd-H calculations.

All constants are CODATA 2018 values. Energies from CP-PAW protocols are in
Hartree, lengths in the structure files are in units of the lattice constant.
"""

HARTREE_EV = 27.211386245988
"""Hartree in electronvolt."""

BOHR_ANG = 0.529177210903
"""Bohr radius in angstrom."""

KB_EV = 8.617333262e-5
"""Boltzmann constant in eV per kelvin."""

EV_J = 1.602176634e-19
"""Electronvolt in joule."""

AMU_KG = 1.66053906660e-27
"""Atomic mass unit in kilogram."""

HBAR_EV_S = 6.582119569e-16
"""Reduced Planck constant in eV s."""

MASS_H_AMU = 1.00794
"""Mass of a hydrogen atom in atomic mass units (the real mass, not the fictitious
mass used for the Car-Parrinello dynamics in the CP-PAW input)."""

A_PD_ANG = 3.89
"""Lattice constant of fcc palladium used in every calculation, in angstrom."""

D_EXP_298K_M2_S = 3.248e-11
"""Experimental diffusion coefficient of H in Pd at 298 K in m^2/s, Powell and
Kirkpatrick, Phys. Rev. B 43, 6968 (1991), as quoted in the project description."""

EA_EXP_EV = 0.23
"""Commonly quoted experimental activation energy of H diffusion in Pd in eV
(Voelkl and Alefeld, Hydrogen in Metals I, 1978)."""
