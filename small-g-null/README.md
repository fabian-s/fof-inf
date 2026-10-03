# Twenty curves and a zero surface (S1.2, S1.8)

Two blocks of `design_table()` (`R/design.R`), run with the study's code:
- `small_G` (cells 124–144): the main factorial (Gaussian low/mid/high, Poisson and binary mid/high;
  errors iid, ou, fpc) with G = 20 curves, the first 20 of each replicate's draws (common random numbers
  with the G = 40 and G = 100 settings).
- `null_beta` (cells 145–152): β ≡ 0 (`truth = "null"`, `R/truth.R`), otherwise the smooth-truth setting
  at the default signal (γ, the functional intercept and the Gaussian error SD unchanged; the binary
  intercept is set again for the 0.2 response level); Gaussian and binary, errors iid and fpc, G = 40 and 100.

Runs (R = 200): `sbatch --array=... slurm/array.slurm <chunk> <run> 124:152 ../r4-results` for the study's
tasks, and `satterthwaite/satt-ext.slurm` for their Satterthwaite degrees of freedom.

Analysis: `Rscript small-g-null/analyse.R <results_dir> <satt_results_dir> small-g-null/summaries`
(`cells.csv`, `g-ladder.csv`, `g-ladder-satt.csv`, `null.csv`, `failures.csv`); relative error with
`RELERR_DIRS=<results_dir> Rscript analysis/relative-error.R synthetic small-g-null/summaries/relative-error.csv`.
