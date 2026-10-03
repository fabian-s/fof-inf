# Satterthwaite degrees of freedom for CL2 intervals (Appendix A.7 of the paper; S1.2, S1.7, S2.1)

Critical values for REML + CL2 (exact, Bayesian) pointwise intervals: z, t_{G−1}, and per-point
Satterthwaite (Bell–McCaffrey) degrees of freedom from refund 79a346fb (`df_gram` full and diagonal; the
full df rescaled to the Bayesian variance). Main factorial at the default signal, R = 200, the study's
replicates; design, reproduction check and task claiming as in `../ablation/`.

- `satt-lib.R`: degrees of freedom through `refund:::pffr_df_context()` / `pffr_influence_df()` on the
  study's evaluation matrices; checks that the df context reproduces the CL2 covariance.
- `run-satt.R`, `satt.slurm`: workers (`sbatch --array=1-50 satterthwaite/satt.slurm`).
- `check-coef.R`: the df on the study's contrasts equal the df column of
  `coef.pffr(crit = "satterthwaite")` (`summaries/check-coef.txt`).
- `analyse-satt.R`, `analyse.slurm` → `summaries/satt-cells.csv`, `satt-contrasts.csv`, `satt-df.csv`,
  `satt-failures.csv`, `SUMMARY.md`.

## Further settings and the plasmode cells

The same records and critical values for every synthetic setting outside the main factorial (cells
43–152, including the 20-curve and zero-surface settings of `../small-g-null/`) and for the 120 plasmode
cells with curve sign flips and attached residuals (nine datasets, G = 40 and all subjects), R = 200 each.

- `run-satt-ext.R`, `satt-ext.slurm`: workers, modes `syn` (cells as an R expression) and `plasmode`
  (datasets). Each record is checked against the study's stored record where one exists.
- `analyse-satt-ext.R`, `r4-analyse.slurm` → `summaries/satt-{syn,plasmode}-*.csv`: cells, contrasts, df,
  failures, check (z arm vs the study's REML + CL2 summary), info.
