# Estimation and inference in penalized function-on-function regression with dependent residual curves

Code and aggregate results for the paper by Fabian Scheipl (LMU Munich, MCML).
Paper, supplement and the detailed write-up of all results: <https://fabian-s.github.io/fof-inf/>.

The study evaluates point estimates and pointwise confidence intervals of `pffr()` (R package refund)
when residual curves are dependent within curves: a simulation study (152 settings, 200 replicates each)
and a plasmode study on nine real datasets, plus follow-up experiments on covariance components,
degrees of freedom, competing methods and intervals centred at the NCV fit.

## Contents

| folder | content |
|---|---|
| `R/` | simulation code: design table, data generation, fits, interval inputs, scoring, summaries |
| `scripts/` | setup, the array runner, the algebra check, construction of the frozen `fpc` kernel |
| `data/` | the frozen `fpc` error kernel (`fpc-kernel.rds`) |
| `analysis/` | study summaries, recommendation rule, figures, DTI application, diagnostics |
| `plasmode/` | plasmode study: data sources, cell design, analysis, recommendation rule |
| `slurm/` | SLURM jobs for the main runs and the study library on the cluster |
| `tests/` | unit tests of the simulation code (`tests/run-tests.R`) |
| `summaries/` | aggregate results of the main and plasmode studies |
| `paper/` | the manuscript and the supplement (Quarto) with their helper code |
| `writeup/` | the detailed write-up of all results (`RESULTS.qmd`) with its figures and tables |
| `software/` | the two mgcv 1.9-5 source builds used by the study and the patch between them |
| `docs/` | the project page (GitHub Pages) |
| follow-up experiments | `ablation/` (covariance components), `satterthwaite/` (degrees of freedom), `small-g-null/` (20 curves, zero surface), `misspec/` (truths outside the spline space), `comparators/` (pcre, GLS, curve bootstrap), `competitors/fgee/` (functional GEE), `ncv-rbc/` (intervals centred at the NCV fit), `lambda-correction/` (smoothing-parameter uncertainty), `ncv-jack-recheck/` (mgcv's NCV covariance), `ffpc-pilot/` (principal component term); each with its own README |

`DESIGN.md` maps the study design to the code. Every entry script is run from the repository root.
Raw per-task records are not included (they are large and, for the cardiology datasets, derived from
patient data); the summary files are. The four largest summary tables are stored gzipped (`*.csv.gz`;
the analysis scripts write them as `.csv`). Per-replicate reproduction checks against the study's
records (`*-crn.csv`) are included only where a document reads them; the other experiments report
their maxima in their `SUMMARY*.md` files.

## Rendering the documents

The paper, the supplement and the write-up read only the summary files in this repository (no model
fitting). They need R with readr, dplyr, tidyr, tibble, knitr, kableExtra, ggplot2, patchwork and withr,
Quarto (≥ 1.7) and XeLaTeX with the Libertinus fonts:

```
cd paper && quarto render paper.qmd && quarto render supplement.qmd
cd writeup && quarto render RESULTS.qmd   # HTML and PDF; also writes writeup/figures/ and writeup/tables/
```

## Software of the study

R 4.6.1 with single-threaded BLAS (`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`); the full package list is
`paper/data/restart-lib-versions.txt`.
- refund: development branch `pffr-refactor` at commit `79a346fb8cc1220a2a6fd2dc8dd3e33166a21061`
  (`remotes::install_github("refunders/refund@79a346fb")`). Scripts check the build: set
  `REFUND_SHA` to this commit, because a GitHub installation does not record it.
- mgcv: pre-release 1.9-5 in two builds, `software/mgcv_1.9-5.tar.gz` (main synthetic factorial) and
  `software/mgcv_1.9-5-patched.tar.gz` (all later runs); `software/mgcv-1.9-5-ncv-vc-and-cxwyd.patch` is
  the difference; SHA-256 checksums in `paper/data/mgcv-checksums.txt`. Supplement S6 says which runs used
  which build; the patch does not touch the REML, NCV and CL2 computations of the study.
- fastFGEE 0.2.2 (CRAN) for `competitors/fgee/`.

Install, for example, with `R CMD INSTALL software/mgcv_1.9-5-patched.tar.gz` and then refund from
the commit above. `slurm/install-restart-lib.R` builds the same library into `~/R/restart-lib`
on the cluster (set `RESTART_SRC` to a directory with the mgcv tarball and a refund source tarball
of the commit); entry scripts put the library in `RESTART_LIB` first on the library path.

## Reproducing the results

The runs used the LRZ Linux cluster (SLURM, serial partition, one core and 2.5 GB per task). Jobs run
from a copy of the repository inside a project directory; raw per-task records go to directories beside
that copy (`../results`, `../<experiment>-results`) and are not part of this repository. The job files
in `slurm/` and the experiment folders load R 4.6.1 from a conda environment (`$HOME/envs/r461`) and
the study library from `$HOME/R/restart-lib`; adapt these lines to another cluster.

1. **Synthetic study.** `sbatch --array=1-<n> slurm/array.slurm <chunk> <run> [cells]`, two runs of
   100 replicates per setting (`scripts/run-array.R`). One setting locally, for example setting 12:
   ```
   export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 REFUND_SHA=79a346fb8cc1220a2a6fd2dc8dd3e33166a21061
   for run in 1 2; do RESTART_CELLS=12 SLURM_ARRAY_TASK_ID=1 Rscript scripts/run-array.R 100 $run results-12; done
   Rscript -e 'source("scripts/setup.R"); readr::write_csv(summarise_study("results-12")$cells, "cells-12.csv")'
   ```
2. **Summaries.** `Rscript analysis/run-analysis.R <results dirs, comma-separated> summaries/final`, then
   `analysis/rule-variants.R`, `analysis/concurvity.R`, `analysis/relative-error.R`,
   `analysis/log-sp-spread.R`, `analysis/ncv-interval-diagnosis.R`, `analysis/timing-benchmark.R`,
   `analysis/scat-convergence.R` and `analysis/write-report.R` (each documents its call in its header).
3. **Plasmode study.** Public datasets: `plasmode/export-public.R` writes them as plain `.rds` files
   (`PLASMODE_PUBLIC_DATA`, default `../plasmode-public-data`); then per dataset
   `slurm/plasmode-public.slurm <app> <n_workers>` (a preparation job, then the array) and
   `plasmode/analyse-plasmode.R` over the nine datasets (→ `summaries/plasmode/`). The two-dataset design:
   `slurm/plasmode-running.slurm`. Fitted truths and task records go to `PLASMODE_CACHE` (default
   `../plasmode-cache`), outside the repository.
4. **NCV-centred hybrid intervals and coverage maps.** `analysis/hybrid-intervals.R` (→ `summaries/hybrid/`).
5. **DTI application.** `analysis/dti-application.R`, `analysis/dti-application-case.R` (→ `summaries/dti*/`).
6. **Studentized errors and design effects** for the write-up and the supplement:
   `writeup/lrz-zstats.slurm`, `writeup/lrz-describe-cells.slurm` (→ `paper/data/`, `writeup/tables/`).
7. **Follow-up experiments:** see each folder's README.

## Data

The DTI, running, gait, 8-lead ECG, ocean, weather and electricity data are public (sources in Appendix
B.2 of the paper and in `plasmode/export-public.R`). The two cardiology datasets (ECG strain study and
atrial-fibrillation trial, N. Bouchahda's group) are patient and volunteer data and are not shared; this
repository contains only aggregate results for them, so their plasmode results cannot be reproduced from
public data. Their loaders in `plasmode/sources.R` read the locations given by `PLASMODE_ECG_DIR` and
`PLASMODE_AF_FILE`.

## How a task flows through the code

A *cell* (a setting) is one row of `design_table()` (`R/design.R`): generating family, fitted family,
error process, G, signal level, truth, D, model (M1 = ff + z, M2 = f(x,t) + z), extra fits.
A *task* is (cell, replicate); `task_table()` enumerates them replicate-major. `run_task(task_id)`:

1. `simulate(cell, rep)` (`R/dgp.R`) draws the data. All random draws depend on `rep` only
   and are made at full size (G_MAX curves, dense t-grid), then subset to the cell's G and D:
   common random numbers across cells, so contrasts between cells are paired.
2. `fit_model(sim, cell, method)` (`R/fit.R`) fits pffr by REML and NCV (plus the cell's extra
   fits), always with `sandwich = "none"`.
3. `interval_inputs(fits, sim, cell)` (`R/fit.R`) evaluates every estimand as L θ̂ with
   model-based and CL2 standard errors per fit (link scale).
4. `score_task(...)` (`R/score.R`) packs point-level records with the truth; `run_task`
   saves them atomically as `task-c<cell>-r<rep>.rds`, with the refund SHA and mgcv version.
5. `summarise_study(dir)` (`R/summarise.R`) builds intervals, coverage, widths, interval scores,
   MSE decompositions, paired contrasts and MC SEs from the saved records.

The plasmode study reuses this path with `run_task(..., simulate_fun = plasmode_simulate)`.

## Contracts between files

`simulate(cell, rep)` returns a list:
- `data`: data.frame for pffr. `Y` is a G × D matrix (`I()`); M1 has `X` (G × 51 matrix, `I()`)
  and `z`; M2 has `x` and `z`.
- `yind` (the cell's t-grid, length D), `xind` (`S_GRID`).
- `test`: data.frame of the frozen test cohort (N_TEST rows; `X`/`x` and `z`, no `Y`).
- `truth`: named list of numeric vectors, one per estimand of the cell's model, each in the
  order of `eval_points(estimand)`:
  - `mean`: test-cohort conditional mean on the **response** scale, N_TEST × EVAL_T;
  - `beta` (M1): β(s,t) on `EVAL_S × EVAL_T`; `f` (M2): f(x,t) on `EVAL_X × EVAL_T`, centred on
    this dataset's x_i (Σ_i f(x_i, t) = 0, as pffr fits it);
  - `alpha` (M1 only): full intercept on EVAL_T; `gamma`: γ(t) on EVAL_T.
  For `warp` cells these are the observed-clock targets.
- `diag`: small list of per-dataset diagnostics.

`eval_points(estimand)` (`R/truth.R`) returns the data.frame of evaluation points in vector
order (first column varies fastest): `mean` → (curve, t); `beta` → (s, t); `f` → (x, t);
`alpha`, `gamma` → t.

## Tests

`Rscript tests/run-tests.R` from the repository root (testthat; synthetic fixtures, no data needed).

## Citation

Scheipl, F. (2026). Estimation and inference in penalized function-on-function regression with
dependent residual curves. Manuscript. <!-- TODO: journal / arXiv reference -->
