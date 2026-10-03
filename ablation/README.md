# Covariance components (Appendix A.6 of the paper)

What each part of CL2 contributes, on the 18 main-factorial settings at the default signal (G = 40 and
100, R = 200, the study's replicates): CR1 vs CL2 with the I − H approximation vs exact CL2; Bayesian vs
frequentist form; z vs t_{G−1} critical values; and the cost of per-replicate smoothing-parameter
selection (a refit at the setting's geometric-mean REML smoothing parameters). CR1 is refund's
`sandwich = "cluster"` without a small-sample factor.

- `make-sp.R`: geometric-mean REML smoothing parameters per penalty from the study summary
  (`summaries/final/study.rds`, written by `analysis/run-analysis.R`) → `sp-fixed.csv`.
- `ablation-lib.R`, `run-ablation.R`: worker; claims tasks through lock directories and checks every REML
  fit against the study's stored record of the same replicate.
- `ablation.slurm`: `sbatch --array=1-60 ablation/ablation.slurm` (study task records in `../results/<run>`).
- `analyse-ablation.R`, `analyse.slurm` → `summaries/`: `ablation-cells.csv` (cell × estimand × arm),
  `ablation-contrasts.csv` (paired contrasts), `ablation-failures.csv`, `SUMMARY.md` (generated tables).
