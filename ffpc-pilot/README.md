# A principal component term for the functional covariate (S1.9)

Inference for β(s,t) from pffr's `ffpc()` term (regression on the leading FPC scores of X,
β(s,t) = Σ_k ψ_k(s) b_k(t)), fitted by REML with exact Bayesian CL2, next to the study's ff fits on the
same simulated data (common random numbers). Rules for the number of components K: ffpc's
variance-explained rule (pve 0.95, 0.90; at most 15 components) and fixed K = 3, 5.

- `ffpc-pilot.R`: first run, settings 11, 12 (Gaussian, smooth errors, G = 40 and 100) and 107 (low-rank
  covariate), 20 replicates, pve 0.99 / 0.95 / 0.90 → `ffpc-pilot-summary.csv`.
- `run-ffpc.R`, `ffpc.slurm`: workers for settings 5, 11, 12, 35, 67, 107, 126 and 129 × rules
  {pve 0.95, pve 0.90, K = 3, K = 5}, R = 200; the study's stored ff REML and NCV records of the same
  replicates are read from `FFPC_REF_DIRS`.
- `summarise.R` → `ffpc-lrz-summary.csv` (from the cluster run): retained components, relative RMSE of β, coverage of the true β and of its projection onto the
  estimated FPC span, half-width relative to ff REML + CL2, detection.
