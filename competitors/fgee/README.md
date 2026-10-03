# Functional GEE (fastFGEE) vs REML + CL2 (S3.3)

Comparison with the penalized functional GEE of Loewinger, Levis, Cui & Pereira (arXiv 2506.20437v3;
R package fastFGEE 0.2.2 from CRAN). Three experiments, R = 200 per cell:

- **A**, their simulation design (their §4, their data generator): function-on-scalar model, Gaussian,
  Poisson and binomial responses through a Gaussian copula with AR(1) × exchangeable latent correlation,
  N ∈ {25, 100, 500} subjects × n_i ∈ {5, 1} curves, their `pffr` call. Every interval (their pffr
  sandwich baselines, refund's cluster sandwiches and CL2, fastFGEE's one-step intervals) is computed on
  the same fit and scored through one evaluation matrix per coefficient function, clustered by subject.
  Their pffr-sandwich baseline adds mgcv's smoothing parameters (on the summed scale) to a bread averaged
  over clusters; arm `fgee_base_fixed` rescales the penalty.
- **B**, the study's 18 main-factorial settings at the default signal (the study's replicates): fastFGEE
  one-step started from the study's REML fit, cluster = curve, working correlations independence, AR(1)
  or FPCA; sandwich Wald intervals from its covariance for β and the mean, and its calibrated critical
  values for α and γ; paired estimation comparison with the study's REML and NCV fits on the same data.
- **P**, plasmode: the four public datasets with independent curves (running, DTI, gait, 8-lead ECG), MID
  truth, REML residuals attached, curve flips, G = all and 40; the two cardiology datasets (de-identified
  arrays, run on the cluster only) in `summaries/cardio/`.

Files:
- `fetch-upstream.sh`: downloads the authors' data generator (github.com/gloewing/fgee_onestep,
  `simulations/functions/sim_functional_data.R` at commit de65400; not redistributed here) to
  `upstream-sim_functional_data.R`; `fgee-dgp.R`: experiment A's design and fits.
- `run.R`: worker for A, B and P (`Rscript competitors/fgee/run.R <A|B|P> <worker> <n_workers> <out_dir>
  [max_rep]`; tasks claimed through lock directories); `fgee.slurm`, `fgee-cardio.slurm`. fastFGEE is
  loaded from `~/R/fastfgee-lib`, refund and mgcv from the study library.
- `smoke-a.R`, `smoke-b.R`: single-data-set checks (the evaluation matrices reproduce fastFGEE's fitted
  curves and SEs).
- `analyse.R` → `summaries/{A,B,P}-summary.csv` (coverage, width, MSE, SE/SD per method and estimand,
  with MC SEs), `A-pooled.csv`, `*-timing.csv`, `*-errors.csv`, `*-method-errors.csv`; `report.R` →
  `summaries/REPORT.md`.
- `pair-ncv.R` (per-replicate pairing with the study's REML and NCV fits on identical data sets; checks
  that the REML estimates agree) and `paired-summary.R` (paired MSE ratios with bootstrap intervals) →
  `summaries/*-paired-ncv-summary.csv`.
