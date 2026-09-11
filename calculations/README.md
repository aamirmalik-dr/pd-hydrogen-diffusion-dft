# CP-PAW calculations

Complete inputs and protocol output of every run used in the project. Protocol files are gzipped (`.prot.gz`); the parsers in `pdhdiff.cppaw_io` read them directly. Restart files, trajectories, projector data and band data are not committed (they are binary and large); everything needed to redo a run from its input is here, and everything needed to reproduce the analysis is in the protocols.

| Directory | Content |
|---|---|
| `01_pd_primitive/` | Bulk Pd in the one-atom fcc cell: `pd.strc`, `pd.cntl`, `pd.prot.gz`; DOS control and protocol (`pd.dcntl`, `pd.dprot`) with the `dos/` output of `paw_dos.x` (total, s, p, d, per spin); band-structure control and protocol (`pd.bcntl`, `pd.bprot`) with `pd_bands.dat` along Gamma-X-W-L-Gamma-K-X |
| `02_pd_supercell/` | 2 x 2 x 2 cubic supercell, 32 Pd, no hydrogen, plus the final structure as `pd_super.xyz` |
| `03_h_octahedral/` | One H at (0, 0, 1/2) a. `stage1_wavefunctions/` optimises the wave functions with frozen atoms, `stage2_relaxation/` restarts from it and relaxes all atoms. The stage-2 protocol contains two consecutive runs appended by the code; the parsers use the last one |
| `03_h_tetrahedral/` | Same for H at (1/4, 1/4, 3/4) a |
| `04_path_octa_to_tetra/g_0.0` to `g_1.0` | Eleven constrained calculations with H placed at O + g (T - O) and the linear constraint held (`MOVE=F`), two stages each as above, with the relaxed structure as `pd_super_octa.xyz` |

## Settings shared by all runs

```
!CONTROL
  !FOURIER EPWPSI=30. CDUAL=2.0 !END          plane-wave cutoff 30 Ry
  !DFT TYPE=10 !END                             PBE
  !PSIDYN ... !AUTO ... !END !END               damped wave-function dynamics with automatic friction
  !MERMIN T[K]=0. TETRA+=T ADIABATIC=T RETARD=20. STOP=T STARTTYPE='N' !END
!END
```

Stage 2 adds `!RDYN FRIC=0.01 STOP=T !AUTO ... !END !END` (atoms propagated) and `START=F` (restart from the stage-1 wave functions). The structure files set `LUNIT[AA]=3.89`, `!KPOINTS R=20.`, `EMPTY=5`, the `NDLSS` augmentation setups for Pd (10 valence electrons) and H (fictitious mass `M=2.` for the dynamics), and the atom list in units of the lattice constant.

## Scrubbing

Two strings were replaced before committing: the cluster account name under which CP-PAW was compiled (now `cluster-user`) and the compute-node name in the MPI report (now `compute-node`). Nothing else in the files was altered; energies, positions and timings are as written by the code.

## Provenance

CP-PAW development version, git hash 460e590 of https://github.com/cp-paw/cp-paw, compiled 31 August 2026; runs executed 7 to 8 September 2026 with 16 MPI tasks, 2 to 3 minutes of wall-clock time per 32-atom relaxation.
