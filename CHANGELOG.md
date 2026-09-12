# Changelog

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
