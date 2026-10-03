# Truths outside the spline space (S5.1)

Do REML + CL2 and model-based intervals cover the true surface when β is not in the spline space, and
does a richer basis help? Gaussian, R² 0.5, D = 61, M1, R = 200.

- Truths (`misspec-lib.R`): T0 = the study's projected smooth β; T1 = T0 + a two-dimensional Gaussian bump
  (sd 0.05, centre (0.3, 0.6)); T2 = T0 + Σ_{j=12..16} √2 cos(jπs) · sin(πt) (weakly identified
  directions). The added features carry 20% of ∫Var(∫Xβ), then β is rescaled to the study's signal; α and
  γ as in the study. Scored against the exact functions.
- Cells: truth × {iid G = 100, fpc G = 100, fpc G = 40} × basis {default, xlarge}; arms REML model-based,
  REML + CL2 (Bayesian and frequentist), NCV + CL2. `make_truth()` is masked so that the study's
  `simulate()` builds all truths; T0 is checked against the study's records by common random numbers
  (`summaries/misspec-crn.csv`). T0 with the xlarge basis at G = 100 equals study cells 113/115: replicates
  1–4 are recomputed as a check, the rest are read from the study's records. NCV in the xlarge cells only
  for replicates 1–100.
- `run-misspec.R`, `misspec.slurm`: worker, replicate-major tasks claimed through lock directories.
- `analyse-misspec.R`, `analyse.slurm` → `summaries/`: `misspec-cells.csv`, `misspec-contrasts.csv`,
  `misspec-truths.csv`, `misspec-crn.csv`, `misspec-failures.csv`, `misspec-fits.csv`, `SUMMARY.md`.
