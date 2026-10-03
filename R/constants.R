# Every number the study uses. Nothing else in R/ contains magic numbers.

# Grids -------------------------------------------------------------------------

# Data are generated on the dense t-grid; D = 61 keeps every 4th point (nested).
T_DENSE <- seq(0, 1, length.out = 241)
T_SUBSAMPLE <- list(`61` = seq(1, 241, by = 4), `241` = seq_len(241))
S_GRID <- seq(0, 1, length.out = 51) # functional covariate argument
t_grid <- function(D) T_DENSE[T_SUBSAMPLE[[as.character(D)]]]

# Estimands are evaluated on coarser nested grids, identical for every D.
EVAL_T <- seq(0, 1, length.out = 31)
EVAL_S <- seq(0, 1, length.out = 26)
EVAL_X <- seq(0.05, 0.95, length.out = 19)

# Fitted bases -------------------------------------------------------------------

# The default bases (every cell except the basis-size block). The truths are always
# projected onto these (R/truth.R).
BASIS_T <- list(bs = "ps", k = 12, m = c(2, 1)) # bs.yindex and bs.int
BASIS_FF <- list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(8, 10)) # s, t
BASIS_X <- list(bs = "ps", k = 8, m = c(2, 1)) # M2's s(x); t margin is BASIS_T

# Basis-size extension (2026-09-28), cell column `basis` (M1 only). mgcv's cubic "ps"
# basis of dimension k has k - 3 equal knot intervals on the (slightly widened) data
# range [0, 1], so k' = 2k - 3 halves every interval: the spline spaces are nested,
# default in large in xlarge, and the default-basis truths stay exactly representable.
# Only the over-parameterization changes. s-margins 13 and 23 vs the rich covariate's
# centred rank 24: 1.5 k_s = 19.5 (large) and 34.5 (xlarge, weakly identified in s).
BASES <- list(
  default = list(t = BASIS_T, ff = BASIS_FF),
  large = list(
    t = modifyList(BASIS_T, list(k = 21)),
    ff = modifyList(BASIS_FF, list(k = c(13, 17)))
  ),
  xlarge = list(
    t = modifyList(BASIS_T, list(k = 39)),
    ff = modifyList(BASIS_FF, list(k = c(23, 31)))
  )
)
# Plasmode ECG 8-lead: with k = 12 neither the functional
# intercept nor the LBBB effect gamma(t) can follow the QRS complex (share of the REML
# residual SS that is mean-curve misfit: 22% overall, 28% by group). Intercept k = 40 and
# bs.yindex (gamma) k = 30 bring it to 0.2% overall and 2.1% by group, so mean structure
# does not leak into the frozen truth's residuals. `int` is the intercept basis (default `t`).
BASES$ecg8 <- list(
  t = modifyList(BASIS_T, list(k = 30)),
  ff = BASIS_FF,
  int = modifyList(BASIS_T, list(k = 40))
)

# Software (asserted in every task) -------------------------------------------------

REFUND_SHA_EXPECTED <- "79a346fb8cc1220a2a6fd2dc8dd3e33166a21061" # pffr-refactor
MGCV_VERSION <- "1.9.5" # Simon Wood's pre-release tarball
MGCV_TARBALL_MD5 <- "d00c4b40db25559417129ba744f9e16a"

# Replicates and seeds -----------------------------------------------------------

N_REP <- 200 # two runs of 100
G_MAX <- 100 # data are drawn for G_MAX curves; G = 40 uses the first 40
N_TEST <- 50 # frozen test cohort
SEED_BASE <- 20260925 # rep r uses SEED_BASE + r in every cell (common random numbers)
SEED_TEST <- 4711 # test cohort
SEED_POPULATION <- 4712 # population sample for GLM intercepts
N_POPULATION <- 2000
LEVEL <- 0.95 # nominal level of every interval

# Concurvity check (R/fit.R): run on the REML fit of these replicates only, one per LRZ run.
CONCURVITY_REPS <- c(1, N_REP / 2 + 1)
# Warn above this worst-case concurvity: 0.8 is the usual "serious" level (the term is then
# nearly a linear combination of the others); the designs' values are far below it.
CONCURVITY_MAX <- 0.8

# Numerical settings ---------------------------------------------------------------

PROJECTION_N <- 101 # grid points for least-squares projection of truths
RANK_PVE <- 0.995 # effective rank = components explaining this share (centred)
INTERCEPT_RANGE <- c(-10, 10) # root search for GLM intercepts

# Recommendation rule and summaries ---------------------------------------

PROPOSAL <- "ncv_cl2_bias"
FALLBACK <- "reml_cl2"
ADEQUATE <- 0.02 # calibration error of an adequate recipe
PRACTICAL_RATIO <- 0.95 # interval-score ratio that counts as a practical difference
N_BOOT <- 1000
SEED_BOOT <- 1

# Covariates ----------------------------------------------------------------------

# X_i(s) = xi_0 + sum_{j=1}^{29} xi_j sqrt(2) cos(j pi s), xi_j ~ N(0, (j+1)^-2).
# The constant term keeps a constant-in-s part of beta identifiable.
X_N_TERMS <- 30
X_SD <- 1 / seq_len(X_N_TERMS)
# Low-rank covariate (extension 2026-09-26): variance (j+1)^-3, matching the public running
# data's knee angle (first PC 83%, 99.5% rank 9 < 1.5 k_s), so beta is weakly identified.
X_SD_LOWRANK <- seq_len(X_N_TERMS)^-1.5

# Truths --------------------------------------------------------------------------

SIGNAL_SHARE_MAIN <- 0.8 # share of S from the ff (M1) or f(x,t) (M2) term; gamma gets the rest
BUMP_CENTER <- 0.5
BUMP_SD <- 0.1
BUMP_HEIGHT <- 1 # in units of sqrt(S)

# Signal levels: R^2 for Gaussian and scat, sqrt(S) on the link scale otherwise.
SIGNAL_LEVELS <- list(
  gaussian = c(low = 0.2, mid = 0.5, high = 0.8),
  scat = c(mid = 0.5),
  poisson = c(mid = 0.5, high = 1),
  binomial = c(mid = 0.5, high = 1),
  betar = c(mid = 0.5),
  nb = c(mid = 0.5)
)

# Average response level each GLM's intercept is shifted to.
MEAN_LEVEL <- c(poisson = 3, binomial = 0.2, betar = 0.3, nb = 3)
SCAT_DF <- 4
BETAR_PHI <- 20
NB_SIZE <- 5

# Error processes (correlation kernels on [0, 1]) ---------------------------------

OU_RANGE <- 0.089 # ou: DE 8 at D = 61 with 20% nugget
OU_NUGGET <- 0.2
OU0_RANGE <- 0.071 # ou0: DE 8 at D = 61, no nugget (exactly AR(1) on the grid)
# fpc: residual FPCs of the DTI application (scripts/make-fpc-kernel.R): 24 components with
# a long, irregular eigenvalue tail, plus the estimated white noise; DE 10.7 on DTI's grid.
# osc (extension, 2026-09-25): damped cosine, correlations change sign with lag (half of all
# pairs negative, minimum -0.42); true SD / model SE for beta 1.8 at fixed penalties.
OSC_RANGE <- 0.3
OSC_PERIOD <- 0.4
OSC_NUGGET <- 0.2
# Heteroskedastic extension (2026-09-26, Gaussian only): fpc_raw keeps the DTI residual
# variance profile along t (rescaled to average variance 1); *_het multiplies subject i's
# errors by exp(HET_SD z_i - HET_SD^2) (average variance 1 for z ~ N(0, 1)).
HET_SD <- 0.5
FPC_KERNEL_FILE <- normalizePath("data/fpc-kernel.rds", mustWork = FALSE) # sourced from the repository root

# Misregistration: h(t) = t + a sin(pi t) / pi, a ~ N(0, WARP_SD^2) truncated.
# 10% maximum displacement SD (WARP_SD / pi), at high signal, because 5% at the
# default signal is a near-null cell (true SD / model SE 0.93-1.03 for beta;
# 1.9 at this setting, fixed-penalty Monte Carlo over the errors).
WARP_SD <- 0.314
WARP_MAX <- 0.9
WARP_QUAD_N <- 201 # quadrature nodes over a for observed-clock targets

# AR(1) comparator: rho profiled over this grid, then refined.
AR_RHO_GRID <- c(seq(0, 0.9, by = 0.1), 0.95, 0.99)
