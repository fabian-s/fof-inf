# Intervals centred at the NCV fit (S4.1)

Under within-curve dependence curve-blocked NCV estimates β with a fraction of REML's error, but NCV + CL2
intervals undercover. This folder diagnoses why and evaluates constructions that try to repair it
(Supplement S4.1 (a)–(g)); none of them gives intervals around the NCV estimate that hold their level on
the plasmode datasets.

Notation: A = Ve Vp⁻¹ (for a penalized Gaussian fit at fixed λ, E θ̂ = Aθ); V_f / V_b = frequentist /
Bayesian CL2 covariance (V_b = V_f + Vp − Ve).

| arm | centre | covariance |
|---|---|---|
| `rbc_f` | twiced, θ̃ = (2I − A)θ̂ | M V_f M', M = 2I − A |
| `rbc_b` | twiced | M V_b M' |
| `rbc_p` | twiced | M V_f M' + (I − A)(Vp − Ve)(I − A)' |
| `rbc2_f` | two steps, M₂ = 3I − 3A + A² | M₂ V_f M₂' |
| `ncv_ijx` | NCV | CL2 plus the cluster-robust λ̂-uncertainty term with cross term (`../lambda-correction/`) |
| `rbc_p_ijx`, `rbc2_p_ijx` | twiced, two steps | the λ̂ term propagated through M or M₂ (Gaussian only) |
| hybrid | NCV | REML + CL2's pointwise half-width |

Diagnostics from the stored study records (no refits): `analysis/ncv-interval-diagnosis.R` →
`summaries/final/ncv-diag-*.csv`, `summaries/plasmode/ncv-diag-*.csv`.

Code:
- `rbc-lib.R`: twicing and its covariances. `combo-lib.R`: combination with the λ̂ term
  (uses `lambda-correction/lc-lib.R` and `ncv-lib.R`).
- Workers (task claiming through lock directories; plasmode task files stay on the cluster):
  `worker-rbc.R` (synthetic; `rbc.slurm`), `worker-rbc-plasmode.R` (`rbc-pm.slurm`),
  `worker-fixlambda-plasmode.R` (refit at the median NCV λ̂; `fixlam-pm.slurm`), `worker-combo.R`
  (modes `synthetic`, `plasmode`, `misspec`; `combo.slurm`).
- Analysis: `analyse-rbc.R` (coverage, 5% quantile of pointwise coverage, half-width, relative error,
  detection, interval score, the hybrid arm, matched-coverage view), run by `rbc-pm-an.slurm`,
  `fixlam-an.slurm`, `combo-an.slurm`, `combo2-an.slurm` and `score-an.slurm`.

Summaries (`summaries/`):
- `rbc-cells.csv` (18 main-factorial settings at the default signal and 14 rough-truth settings,
  replicates 1–200), `rbc-fresh-cells.csv` (same settings, fresh replicates 201–400, REML + CL2 refitted);
- `rbc-plasmode-cells.csv` (60 plasmode cells: nine datasets × truth × residual source, all subjects,
  replicates 1–100), `fixlam-plasmode-cells.csv` (the same cells refitted at the median NCV λ̂);
- `combo-syn-cells.csv`, `combo-plasmode-cells.csv` (twicing combined with the λ̂ term, Gaussian);
  `combo2-{syn,pm,ms}-cells.csv` (the second combination run with the two-step arms: synthetic, plasmode
  and misspecified-truth cells of `../misspec/`);
- `score-*-cells.csv`: interval scores and matched coverage on paired replicates.
