# mgcv's covariance for curve-blocked NCV fits (S4.2)

mgcv offers a neighbourhood-aware covariance for NCV fits (`nei$jackknife = TRUE`). For refund's
curve-blocked neighbourhoods (one neighbourhood per curve) mgcv 1.9-5 does not compute it correctly: with
`jackknife` unset it returns the Bayesian `Vp`, and with `jackknife = TRUE` it reads the matrix of
point-deletion changes with the wrong row stride, so the returned covariance is about `Vp − Ve` (SEs 0.5–0.8
times those of `Vp`). This folder evaluates both mgcv's result and the intended estimator, computed from
the block fit.

- `jack.R`: mgcv's covariance with `nei$jackknife = TRUE` on main-factorial settings 12/24/36 (Gaussian,
  Poisson, binary; smooth errors, G = 100, default signal), replicates 1–20: same NCV fit, intervals from
  `Vp`, mgcv's jackknife `Vp` and `Vc`, and CL2 → `jack-reps.csv`.
- `jackfix-lib.R`: `jackfix_vcov()` computes from the block-NCV fit what mgcv returns for the equivalent
  per-point neighbourhood form (each point predicted with its whole curve left out; same NCV criterion,
  same λ̂ and coefficients, at D times the cost per NCV evaluation): single-point Newton perturbations,
  deviance-residual ratios at the leave-curve-out fit, `Vcv = Σ_g (Σ_{i∈g} rat_i dd_i)(…)' n/(n − edf)`,
  `alpha = max(tr Vraw₊ / tr Ve, 1)`, `Vj = Vcv + alpha (Vp − Ve)` (psd projection), and mgcv's
  `jackknife = TRUE` result for the per-point form. About 1 s per study fit.
- `check-pointwise-pffr.R`: study-scale check, setting 12, replicate 1 refitted by pffr with per-point
  neighbourhoods: same log λ and coefficients, `Vp` and `Vc` equal to `jackfix_vcov()` up to 1.3e-13
  (`summaries/jackfix-pointwise-check-c12-r1.csv`).
- `run-jackfix.R`, `jackfix.slurm`, `jackfix-pm.slurm`: workers on the study's replicates (settings
  12/24/36, replicates 1–100) and on 42 plasmode cells (six datasets with independent curves, all
  subjects, curve flips, attached residuals, every truth × residual pair, replicates 1–100). Arms on the
  same NCV fit: `ncv_vp` (Bayesian `Vp`), `ncv_cl2` (CL2), `jk_nei` (neighbourhood-aware `Vp`,
  corrected), `jk_nei_vc` (its `Vc`, with the λ̂ correction), `jk_cv` (its sampling part `Vcv` alone),
  `jk_raw` (mgcv's `jackknife = TRUE` for the per-point form); `reml_cl2` (REML + CL2 on the same data).
- `analyse-jackfix.R`, `jackfix-pm-an.slurm` → `summaries/jackfix-{syn,pm}-*.csv` (cells, per-dataset
  aggregates, reproduction check of the study records, failures, run information).

Software: the study library (mgcv 1.9-5 patched build, refund 79a346fb).
