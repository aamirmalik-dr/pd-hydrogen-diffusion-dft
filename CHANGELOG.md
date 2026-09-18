# Changelog

## 0.3.0 (2026-09-18)

- Quantum estimates from the one-dimensional profile: the curvature at the top of the barrier (local cubic and spline) gives the imaginary frequency of the unstable mode, the Wigner tunnelling factor and the crossover temperature; the octahedral well frequency gives the quantum correction of the path mode (`wigner_correction`, `crossover_temperature`, `quantum_well_factor` in `pdhdiff.tst`, block `quantum_estimates` in `results/metrics.json`). Both factors exceed one, so the write-up no longer lists tunnelling as a candidate for the overestimated diffusion constant.
- Comparison with experiment split into a barrier part and a prefactor part (block `comparison_with_experiment`), and the Arrhenius figure shows the measured 298 K value continued with the quoted experimental activation energy.
- Arrhenius figure: legend and summary text no longer overlap the curves. Diagnostics figure: step labels of the octahedral cage atom no longer overlap.
- Walkthrough notebook extended with the comparison above, imports cleaned, re-executed; notebooks are now linted in CI.
- Documentation corrected where it had fallen behind the code: test count in the README, the description of the stage-2 stop criterion in `docs/cppaw_workflow.md`.

## 0.2.0 (2026-09-11)

- Relaxation diagnostics: the periodic atom lists of the protocols are parsed (`parse_reports`), and a new script follows the force on the cage atoms through each run, tabulates the residual forces, the energy drift and a harmonic bound on the unrelaxed energy, and compares the printed Lagrange multiplier with the slope of the profile. Finding: the runs are converged in energy (meV level) but not in force (up to 1 eV/A on the atoms next to H); the write-up now says so and no longer calls the end points "true minima".
- Sensitivity analysis: spread of the barrier and curvature estimators propagated to the diffusion constant, reported in `results/metrics.json` and `results/RESULTS.md`.
- Classical isotope scaling of the attempt frequencies (H, D, T) added to the analysis output with its caveat.
- Trajectory file written in extended-xyz form readable by ASE.
- Test suite extended: xyz files against protocol positions, input generator round trip, figure rendering, sensitivity block.
- Continuous integration on Linux and Windows with a reproduction check of the committed result files.
- Names and places written with their diacritics; package metadata and citation file completed.
- Linter versions pinned to those the code was formatted with (ruff 0.16.7, black 26.5.1); numpy 2 required.

## 0.1.0 (2026-09-11)

- First public release: CP-PAW inputs and protocols, `pdhdiff` package, results, figures, seminar slides.
