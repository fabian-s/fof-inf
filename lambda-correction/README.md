# Smoothing-parameter uncertainty terms for CL2 (S4.1 (e)–(f), S4.3)

CL2 treats the selected smoothing parameters as fixed. This folder adds a cluster-robust term for the
variability of ρ̂ = log λ̂ to the CL2 covariance, built from analytic curve scores of the smoothing
criterion (infinitesimal jackknife, IJ, with and without the cross term between θ̂ and ρ̂), and compares
it with a one-step and an exact curve jackknife, mgcv's own correction, and an oracle that refits at
fixed smoothing parameters. Derivations: `NOTE.qmd` (REML), `NOTE-ncv.qmd` (Gaussian NCV),
`NOTE-ncvg.qmd` (Poisson and binomial NCV).

**REML** (S4.3):
- `lc-lib.R`: `lc_pieces()` (fit pieces, checked against mgcv's `Vp` and `db.drho`), `ij_scores()`
  (analytic curve scores of the REML gradient, Gaussian/Poisson/binomial, canonical links),
  `onestep_loo()`, `jackknife_refits()` (exact curve jackknife), `theta_fixed_sp()`, `lc_covariances()`
  (all arms for one `sandwich = "none"` REML fit). `lc-run-lib.R`: task and worker loop.
- `validate.R`: single-replicate checks (finite differences in curve weights; IJ and one-step vs exact
  jackknife) → `validation/`.
- `run-lc.R`, `lc.slurm` (`lc-v2.slurm`: without the exact jackknife): 18 main-factorial settings at the default signal × 200
  replicates, common random numbers with the study records; fixed-λ arm with
  `../ablation/sp-fixed.csv`. `analyse-lc.R` → `summaries/` (cells, contrasts, calibration, targeting,
  regime, times, `SUMMARY.md`); `summaries-v2/`: a second run of the six dependent G = 40 settings with
  regularised arms (replicates 1–100).
- Plasmode: `run-lc-plasmode.R`, `lc-plasmode*.slurm`, `analyse-lc-plasmode.R` → `summaries-plasmode/`
  (IJ and one-step arms on the plasmode cells, replicates 1–200).
- `flag-check.R`: does a large per-fit SD(ρ̂) or a small λ̂ flag replicates with a high relative error
  of β or CL2 misses (cited in the write-up).
- `make-report.R` → `REPORT.md`.

**NCV** (S4.1 (e)–(f)):
- `ncv-lib.R` (Gaussian: criterion, exact gradient and Hessian, analytic curve scores, IJ, one-step, exact
  re-minimisation), `ncvg-lib.R` (Poisson and binomial: one Newton step per left-out curve, directional
  derivatives in ρ and the curve weights; components with negligible curvature treated as fixed);
  checks `ncv-validate.R`, `ncvg-validate.R`, `ncvg-sensitivity.R`.
- `make-ncv-sp.R` → `ncv-sp-fixed.csv` (geometric-mean NCV smoothing parameters for the fixed-λ arm).
- `run-lc-ncv.R`, `lc-ncv*.slurm` → `summaries-ncv/` (6 Gaussian settings × 200);
  `run-lc-ncvg.R`, `lc-ncvg*.slurm` → `summaries-ncvg/` (12 Poisson and binomial settings × 200);
  `analyse-lc-ncv.R` for both. Plasmode: `run-lc-ncv-plasmode.R`, `lc-ncvp*.slurm`,
  `analyse-lc-ncv-plasmode.R` → `summaries-ncv-plasmode/` (130 cells, replicates 1–100).
- `make-report-ncv.R`, `make-report-ncvg.R` → `REPORT-ncv.md`, `REPORT-ncvg.md`.

The report scripts also read the per-replicate reproduction checks against the study records
(`*-crn.csv`, written by the analysis scripts). These are not included; their maxima are stated in the
reports and in the `SUMMARY*.md` files.
