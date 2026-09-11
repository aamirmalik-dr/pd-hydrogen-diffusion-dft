# Results

All values are produced by `python scripts/run_all.py` from the committed CP-PAW protocol files. Energies from CP-PAW are in Hartree (1 H = 27.2114 eV).

## Bulk palladium

| Quantity | Primitive cell | 2 x 2 x 2 supercell |
|---|---|---|
| Atoms | 1 | 32 |
| k-points (irreducible) | 63 | 8 |
| Bands | 10 | 166 |
| Plane-wave cutoff | 30 Ry | 30 Ry |
| Total energy | -30.0020008 H | -960.2619568 H |
| Energy per atom | -30.00200 H | -30.00819 H |
| Chemical potential | 21.435 eV | 21.654 eV |
| Steps to convergence | 251 | 457 |

The supercell energy per atom is 168 meV lower than the primitive-cell value. Both runs use the same k-point density, but they sample different points of the same band structure, and for a metal with a sharp d peak at E_F the coarse grids are not converged with respect to each other. Every quantity below is a difference between runs with the identical supercell and grid, so this offset cancels.

Electronic structure of the primitive cell (`electronic_structure.json`): d band from -5.14 eV to E_F, DOS at E_F 4.43 states/eV/atom (both spins), occupied total DOS integrates to 9.96 electrons.

## Hydrogen at the interstitial sites (unconstrained relaxations)

| Site | Total energy | Relative | Pd-H distances | Steps |
|---|---|---|---|---|
| Octahedral, H at (0, 0, 1/2) a | -960.8412761 H | 0 | 1.955 to 1.966 A (six neighbours, ideal 1.945 A) | 337 |
| Tetrahedral, H at (1/4, 1/4, 3/4) a | -960.8381902 H | +83.97 meV | 1.753 A (four neighbours, ideal 1.684 A) | 378 |

In both runs the H atom stays at the high-symmetry position and the energy is stationary at the end. The structures are converged in energy but not in force; see the section on convergence below before quoting the geometries to better than 0.02 A.

## Energy profile along the reaction coordinate

`profile.csv`, energies relative to the g = 0 point, all from the stage-2 relaxations.

| g | Constraint value (bohr) | E (H) | E - E_O (meV) | Max Pd displacement (A) | Steps |
|---|---|---|---|---|---|
| 0.0 | -4.90 | -960.8412762 | 0.0 | 0.016 | 337 |
| 0.1 | -4.17 | -960.8410758 | 5.5 | 0.021 | 355 |
| 0.2 | -3.43 | -960.8404399 | 22.8 | 0.031 | 348 |
| 0.3 | -2.70 | -960.8390136 | 61.6 | 0.051 | 392 |
| 0.4 | -1.96 | -960.8370711 | 114.4 | 0.068 | 368 |
| 0.5 | -1.23 | -960.8352050 | 165.2 | 0.084 | 364 |
| 0.6 | -0.49 | -960.8339712 | 198.8 | 0.095 | 361 |
| 0.7 | 0.245 | -960.8339921 | 198.2 | 0.098 | 360 |
| 0.8 | 0.98 | -960.8352027 | 165.3 | 0.094 | 368 |
| 0.9 | 1.72 | -960.8370149 | 116.0 | 0.083 | 375 |
| 1.0 | 2.45 | -960.8381737 | 84.4 | 0.069 | 379 |

Checks: the constraint value printed by CP-PAW equals a (g - 2/3) in bohr for every point; the end points agree with the unconstrained site energies to -0.003 meV (O) and +0.45 meV (T); the H atom deviates from the straight O to T line by less than 0.007 A.

## Fits

| Quantity | Value | Method |
|---|---|---|
| Transition state | g = 0.650, E = 203.8 meV | cubic polynomial through the four points around the maximum |
| Transition state, cross-check | g = 0.648, E = 203.0 meV | not-a-knot cubic spline through all eleven points |
| Highest calculated point | g = 0.6, E = 198.8 meV | raw |
| E_act (O to T) | 203.8 meV | |
| E_act (T to O) | 119.4 meV | |
| k_O | 0.968 eV/g^2 | E = k/2 g^2 + c g^4 on g <= 0.3, rms residual 0.2 meV; two-point estimate 1.09 |
| k_T | 7.31 eV/g^2 | E = E_T + k/2 d^2 + c d^3, d = g - 1, on g >= 0.7, rms residual 1.7 meV; two-point estimate 6.31 |

Path length |T - O| = sqrt(3) a / 4 = 1.6844 A converts k to force constants; the real hydrogen mass 1.00794 u gives the attempt frequencies.

| Site | k (N/m) | omega_0 (rad/s) | nu_0 (Hz) | hbar omega_0 (meV) |
|---|---|---|---|---|
| Octahedral | 5.47 | 5.72 x 10^13 | 9.10 x 10^12 | 37.6 |
| Tetrahedral | 41.3 | 1.57 x 10^14 | 2.50 x 10^13 | 103.3 |

For comparison, inelastic neutron scattering places the local H mode on the octahedral site of PdH_x near 60 to 70 meV; the value here is the soft direction towards the triangle face only.

## Hop geometry

Evaluated numerically from the fcc site lists (`pdhdiff.structure.hop_geometry`): every octahedral site has 8 tetrahedral neighbours and every tetrahedral site 4 octahedral neighbours, all at 1.684 A. The geometric factors sum_j (r_j - r_i)(x)(r_j - r_i) are isotropic, (a^2/2) 1 for hops out of O and (a^2/4) 1 for hops out of T, so D is a scalar.

## Diffusion constant

Populations from the one-dimensional harmonic wells with site multiplicities 1 (O) and 2 (T) per Pd atom; rates Gamma = (omega_0 / 2 pi) exp(-E_a / k_B T); D = 1/2 sum_i P_i sum_j Gamma_ji |r_j - r_i|^2 / 3.

| T (K) | Gamma(T <- O) per pathway (1/s) | Gamma(O <- T) per pathway (1/s) | P_tet | D (m^2/s) |
|---|---|---|---|---|
| 200 | 6.6 x 10^7 | 2.4 x 10^10 | 0.005 | 5.0 x 10^-12 |
| 300 | 3.4 x 10^9 | 2.5 x 10^11 | 0.027 | 2.5 x 10^-10 |
| 400 | 2.5 x 10^10 | 7.8 x 10^11 | 0.059 | 1.8 x 10^-9 |
| 600 | 1.8 x 10^11 | 2.5 x 10^12 | 0.125 | 1.2 x 10^-8 |
| 800 | 4.7 x 10^11 | 4.4 x 10^12 | 0.176 | 2.9 x 10^-8 |
| 1000 | 8.5 x 10^11 | 6.3 x 10^12 | 0.215 | 5.1 x 10^-8 |

| Quantity | Value |
|---|---|
| D(298 K) | 2.40 x 10^-10 m^2/s |
| D(298 K), experiment | 3.25 x 10^-11 m^2/s (Powell and Kirkpatrick 1991) |
| Ratio calculated / experiment | 7.4 |
| Arrhenius fit 200 to 1000 K | E_a = 0.198 eV, D_0 = 5.3 x 10^-7 m^2/s |
| Experimental activation energy (quoted) | about 0.23 eV |
| Detailed-balance ratio of the two sublattice fluxes at 300 K | 1.000 |

Alternative conventions at 298 K, for comparison with the seminar slides and the project description:

| Convention | D(298 K) (m^2/s) |
|---|---|
| Populations with the three-dimensional harmonic prefactor (omega_O/omega_T)^3 | 1.39 x 10^-10 |
| Populations from site energies only | 4.29 x 10^-10 |
| Three-dimensional prefactor and no factor 1/2 (as on the seminar slides) | 2.78 x 10^-10 |

## Sensitivity to the estimators

`metrics.json`, block `sensitivity`. Each fitted quantity has more than one reasonable estimator; the table gives D(298 K) with one estimator changed at a time and the range over all combinations.

| Quantity | Primary | Alternative | D(298 K) with the alternative |
|---|---|---|---|
| Barrier | local cubic, 203.8 meV | cubic spline, 203.0 meV | 2.47 x 10^-10 m^2/s |
| Barrier | | highest calculated point, 198.8 meV | 2.91 x 10^-10 m^2/s |
| k_O | least squares, 0.968 eV/g^2 | two-point, 1.091 eV/g^2 | 2.54 x 10^-10 m^2/s |
| k_T | least squares, 7.31 eV/g^2 | two-point, 6.31 eV/g^2 | 2.39 x 10^-10 m^2/s |
| All combinations | 2.39 x 10^-10 m^2/s | | 2.39 to 3.09 x 10^-10 m^2/s |

The estimator choice moves D(298 K) by at most 30 %, less than the spread between the population conventions (factor 1.8) and far less than the gap to experiment (factor 7).

## Classical isotope scaling

With the same profile and classical harmonic TST only the attempt frequencies change, by sqrt(m_H / m): D_D / D_H = 0.707 and D_T / D_H = 0.578 at every temperature (`metrics.json`, block `isotopes_classical`). The measured isotope effect of hydrogen in Pd is governed by zero-point energy and tunnelling, which this treatment does not contain, so these ratios are the classical reference point and not a prediction.

## Convergence of the relaxations

`relaxation_diagnostics.json` and `figures/relaxation_diagnostics.png`, produced by `scripts/relaxation_diagnostics.py` from the periodic atom lists in the protocols.

CP-PAW's autopilot ends a damped relaxation when the energy and the ionic kinetic energy stop changing, not when the forces vanish. The atom lists printed every 100 steps carry the forces that propagated the atoms, so they show what the runs looked like when they were stopped. Following the tetrahedral cage atom PD03 through the tetrahedral run:

| Step | Pd-H distance (A) | Radial force on PD03, outward (eV/A) | E - E_final (meV) |
|---|---|---|---|
| 645 (start of the relaxation, ideal lattice) | 1.684 | +1.36 | +159 |
| 700 | 1.687 | +1.17 | +144 |
| 800 | 1.703 | +0.65 | +70 |
| 900 | 1.729 | +0.04 | +8 |
| 1000 | 1.752 | -1.08 | 0 |
| 1022 (end) | 1.753 | -0.63 | 0 |

The cage opens, the force crosses zero near 1.73 A, the atoms overshoot while the second shell is still relaxing, and the autopilot's final friction phase brings them to rest with an inward force left on them. The linear fit of force against distance gives a radial stiffness of 31 eV/A^2 for the cage atom (with the surroundings relaxing at the same time, so this is the soft response).

Residual forces at the end of every stage-2 run (1 mH/bohr = 0.0514 eV/A):

| Run | Largest force, free atoms (mH/bohr) | Largest force, constrained atoms (mH/bohr) | Energy drift, last 30 steps (meV) | Ionic T at the end (K) | Harmonic bound on the unrelaxed energy (meV) |
|---|---|---|---|---|---|
| octahedral | 1.5 | | -0.35 | 3 | 0.5 |
| tetrahedral | 12.3 | | +0.22 | 1 | 26 |
| g = 0.0 | 1.5 | 1.4 | -0.35 | 3 | 0.5 |
| g = 0.1 | 8.3 | 11.8 | +0.14 | 2 | 26 |
| g = 0.2 | 1.7 | 3.0 | -0.41 | 5 | 1.4 |
| g = 0.3 | 1.9 | 4.4 | +0.41 | 1 | 0.4 |
| g = 0.4 | 4.0 | 16.3 | +1.23 | 1 | 34 |
| g = 0.5 | 4.3 | 17.1 | +0.33 | 1 | 37 |
| g = 0.6 | 4.8 | 18.6 | +0.38 | 1 | 44 |
| g = 0.7 | 5.6 | 19.6 | +0.38 | 2 | 49 |
| g = 0.8 | 6.0 | 15.6 | +0.19 | 1 | 32 |
| g = 0.9 | 7.5 | 13.6 | +0.16 | 1 | 26 |
| g = 1.0 | 7.2 | 6.3 | +0.03 | 1 | 7 |

The constrained atoms are the H atom and the three triangle atoms; their residual force includes the constraint reaction only for the H atom (the triangle coefficients are small), so the 12 to 20 mH/bohr on the triangle atoms are genuine unbalanced forces, up to 1 eV/A. The last column is the harmonic energy sum F^2 / 2k over the first-shell atoms with the 31 eV/A^2 stiffness, an upper estimate of how far each energy could sit above the fully relaxed value.

That bound is far from tight, and the profile itself shows it. The residual forces vary erratically from point to point (12 mH/bohr at g = 0.1, 3 at g = 0.2, 4 at g = 0.3, 16 at g = 0.4) while the energies of these points lie on the octahedral parabola with an rms residual of 0.2 meV. An energy error proportional to F^2 would put a 25 meV kink between g = 0.1 and g = 0.2; there is none. The same holds at the barrier, where the four points around the maximum agree with a cubic spline through all points to 1 meV. The residual forces therefore lie along stiff coordinates of the cage and cost little energy: the energies of the profile are converged to the meV level, the geometries are not converged to better than about 0.02 A, and the printed forces must not be used as such.

Two consequences. First, the Lagrange multiplier of the constraint, which equals dE/dVAL only at a constrained minimum, undershoots the slope of the fitted profile by 20 to 25 % on the octahedral side (right panel of the figure); its thermodynamic integration gives 47 meV at g = 1 against 84 meV directly. It is a check of the relaxation, not an independent route to the profile. Second, the agreement of the path end points with the unconstrained site runs (0.003 and 0.45 meV) is a reproducibility check of the constrained setup rather than a convergence proof: those pairs of runs start from the same structure and follow the same protocol.

A force-converged rerun of the thirteen stage-2 calculations, continued from their restart files with a stricter stop criterion, is the first thing to do before quoting these energies to better than a few meV.

## What limits the accuracy

In order of expected size: the PBE barrier (tens of meV), the missing zero-point energy (hbar omega differs by 66 meV between the O and T wells along the path alone), the one-dimensional harmonic prefactor, the coarse k-point grid and the 32-atom cell, the incomplete force convergence of the relaxations (meV level in the energies, see above), and the neglect of correlated O to T to O return jumps. None of these was converged or corrected in the project; they are listed so the numbers above are read with the right error bars.
