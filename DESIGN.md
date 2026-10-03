# Study design in the code

The design is described in Appendix B of the paper and Section S2 of the supplement. This file maps it
to the code.

## Synthetic study

**Settings.** `design_table()` (`R/design.R`) has 152 rows (cells), a core factorial plus blocks that
change one factor away from it. Defaults: Gaussian, error `fpc`, G = 100, default signal, smooth truth,
D = 61, model M1, rich covariate, default bases.

| block | factors | cells |
|---|---|---|
| `core` | family {gaussian, poisson, binomial} × error {iid, ou, fpc} × G {40, 100} × signal (3 Gaussian, 2 GLM levels) | 1–42 |
| `warp` | misregistration × family × G {40, 100} × truth {smooth, wiggly}, high signal | 43–54 |
| `dense_grid` | D = 241 × family × {iid, ou, fpc} | 55–63 |
| `rough_truth` | `wiggly` truth × family × {iid, fpc} × signal | 64–77 |
| `term_type` | model M2 (f(x, t) instead of ff) × family × {iid, ou, fpc} | 78–86 |
| `families` | {scat, betar, nb} × {iid, fpc}; nb also fitted as Poisson | 87–94 |
| `ar1_home` | Gaussian, error `ou0` (exactly AR(1) on the grid) | 95 |
| `oscillating` | error `osc` × family | 96–98 |
| `warp_ar1` | Gaussian `warp`, wiggly truth, high signal | 99 |
| `heteroskedastic` | Gaussian × {fpc_raw, fpc_het, fpc_raw_het} × G {40, 100} | 100–105 |
| `lowrank_covariate` | low-rank covariate × family × {iid, fpc} | 106–111 |
| `basis_size` | bases {large, xlarge} × family × {iid, fpc} | 112–123 |
| `small_G` | the core factorial at G = 20 | 124–144 |
| `null_beta` | β ≡ 0 × {gaussian, binomial} × {iid, fpc} × G {40, 100} | 145–152 |

Extra fits inside existing cells (column `extra`): pointwise NCV and AR(1) working models in Gaussian
cells at G = 100 (`mark_extra_fits()`).

**Constants** (`R/constants.R`): grids (t: 241 points, subsampled to D; s: 51 points), bases (`BASES`:
`ps`, k = 12 for the t-direction and the intercept, ff 8 × 10; `large` and `xlarge` halve and quarter
the knot spacing), signal levels per family (`SIGNAL_LEVELS`), error kernels, the test cohort (N_TEST =
50 covariates), evaluation grids and G_MAX.

**Error processes** (`R/dgp.R`). `iid`; `ou` (Ornstein–Uhlenbeck plus 20% white noise); `fpc` (the
residual covariance of the DTI application, frozen in `data/fpc-kernel.rds` by `scripts/make-fpc-kernel.R`);
`ou0` (pure AR(1)); `osc` (sign-changing correlations); `fpc_raw`, `fpc_het`, `fpc_raw_het`
(heteroskedastic along t, across subjects, both); `warp` (misregistration). Non-Gaussian responses use the
kernel as the latent layer of a Gaussian copula, so every marginal is exactly the fitted family.

**Truths** (`R/truth.R`): fixed formulas projected onto the fitted bases (exactly representable),
`smooth` or `wiggly`, scaled to the cell's signal; `null` sets β to zero.

**Fits and intervals** (`R/fit.R`, `R/summarise.R`). Every dataset is fitted by REML and by curve-blocked
NCV (`pffr(method = "NCV")`, one curve per neighbourhood), always with `sandwich = "none"`. Covariances:
model-based (`Vc` for REML, `Vp` for NCV) and CL2 (`sandwich = "cl2"`, exact leverage adjustment, Bayesian
form; the frequentist form is stored too). Estimands: the conditional mean of a frozen test cohort
(response scale), β(s, t) or f(x, t), the full intercept α(t) and γ(t), each as L θ̂ with its truth.
All intervals are pointwise at nominal 95%.

**Metrics** (`R/summarise.R`, `analysis/contrasts.R`): grid-averaged coverage and its 5% quantile over
points, width, interval score, MSE split into squared bias and variance, per-term EDF and log smoothing
parameters, fit times; paired per-replicate contrasts between arms and between cells with Monte Carlo SEs.
The recommendation rule (paired calibration error and interval-score ratio over a pool of cells) is in
`R/summarise.R` and `analysis/contrasts.R`; `analysis/rule-variants.R` evaluates its variants.

**Replicates.** R = 200 per cell, run as two runs of 100. All random draws depend on the replicate only
and are made at full size, so all cells share common random numbers. Tasks are enumerated
replicate-major and looked up by id; each result file records the refund SHA and the mgcv version.

## Plasmode study

`plasmode/plasmode.R` (cells, truths, `plasmode_simulate()`), `plasmode/sources.R` (one loader per
dataset), `plasmode/rule.R` (recommendation rule per dataset), `plasmode/analyse-plasmode.R`
(aggregates). The model is fitted once to each real dataset by REML and by NCV; the fitted means (and
their log-midpoint, `MID`) are frozen as truths, and replicates add the dataset's own residual curves with
random signs, attached to their own subject or permuted. Nine datasets: six with independent curves
(`ecg_gls4b`, `af`, `running`, `dti`, `gait_vc`, `ecg_roahd`) and three with dependent curves (`ocean`,
`aemet`, `adelaide`), with block sign flips and a block-NCV arm for the latter. A further 16-cell design
on ECG strain (`ecg_gls4`) and running includes the AR(1) working model.

## Code rules

- A block is only rows of the design table; there is no block-specific code path.
- Variation is handled by small lookup tables (families, error kernels, models, bases).
- The generating family and the fitted family are separate columns.
- One estimand interface: an evaluation matrix L plus the aligned truth.
- Every constant lives in `R/constants.R`.
- Validation: `tests/` and `scripts/check-algebra.R` (Gaussian, fixed penalties, known Σ: the CL2 and
  model-based SEs against the exact sampling SE).
