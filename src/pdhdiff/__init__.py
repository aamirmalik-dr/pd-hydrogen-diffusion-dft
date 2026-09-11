"""pdhdiff: post-processing of CP-PAW calculations of hydrogen diffusion in fcc palladium.

The package reads CP-PAW protocol and structure files, assembles the energy
profile of a hydrogen atom moving from an octahedral to a tetrahedral
interstitial site, fits the profile, and evaluates jump rates and the
diffusion constant with one-dimensional harmonic transition-state theory.
"""

from pdhdiff import constants, cppaw_io, profile, structure, tst

__all__ = ["constants", "cppaw_io", "profile", "structure", "tst"]
__version__ = "0.1.0"
