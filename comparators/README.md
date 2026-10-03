# Comparators: refund's remedies for dependent residuals vs REML + CL2 (S3.2)

Six Gaussian main-factorial settings (iid/ou/fpc × G 40/100, M1, D = 61) and twelve plasmode cells (six
datasets with independent curves, MID truth, REML residuals, curve flips, G = all/40), with the study's
replicates; REML fits are checked against the study records.

- `comparators-lib.R`: arms. REML model-based and CL2; `pcre()` with FPCA eigenfunctions of the REML
  residual curves (pve 0.95); GLS with the algorithm of `pffr_gls()` (Σ̂ = raw residual covariance or FPCA
  covariance); curve bootstrap with `pffr_coefboot()`'s resampling and refit (percentile and bootstrap-SE
  Wald), B = 199 on replicates 1–100.
- `run-comparators.R`, `comparators.slurm`: worker with lock-directory claiming
  (`sbatch --array=1-N comparators/comparators.slurm <out> <keys|all> <max_rep> [boot_syn] [boot_plas]`).
- `check-coefboot.R`: this bootstrap vs `refund::pffr_coefboot()` on one replicate (same seed and B).
- `analyse-comparators.R`, `analyse.slurm` → `summaries/comparators-cells.csv` (cell × estimand × arm),
  `comparators-contrasts.csv` (paired vs REML + CL2), `comparators-fits.csv`, `comparators-failures.csv`,
  `SUMMARY.md`.
- CL2 on the comparator fits, GLS (raw or FPCA Σ̂) + CL2 and pcre + CL2, on the same replicates
  (`cl2-arms-lib.R`, `run-cl2-arms.R`, `cl2-arms.slurm`, `check-cl2-arms.R`, `analyse-cl2-arms.R`,
  `cl2-arms-analyse.slurm`) → `summaries/cl2-arms-*.csv`, `cl2-arms-SUMMARY.md`.
