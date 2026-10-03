# Experiment A: fastFGEE's own simulation design (Loewinger et al., arXiv 2506.20437v3, §4),
# one replicate = one data set, scored for every interval method on the SAME pffr fit.
#
# Data: their GenerateData() (upstream-sim_functional_data.R, downloaded by fetch-upstream.sh from
# github.com/gloewing/fgee_onestep simulations/functions/sim_functional_data.R) with
# their settings: L = 50, n_i curves per subject, sigma_sq = 10, divFactor = 3,
# AR1(0.75) functional x exchangeable/AR1(0.75) longitudinal Gaussian copula.
#
# Estimands: beta_0(s) (full intercept: scalar + smooth), beta_1(s), beta_2(s) on their
# grid s. All methods are scored through one evaluation matrix L_r per coefficient
# function applied to the coefficient vector in pffr's (constrained) parametrisation,
# so centring constraints cannot differ between methods: fastFGEE works on pffr's own
# design matrix (fit$X) and returns coefficients in the same order (checked below).

source_upstream <- function(dir) {
  env <- new.env()
  sys.source(file.path(dir, "upstream-sim_functional_data.R"), envir = env)
  env$GenerateData
}

# Their settings (fGEE_sim.R): z11 = 5, z12 = 3 (unused in the Kronecker branch),
# sigma_sq = 10, divFactor = 3, rho = 0.75, fn.corr = "ar1".
simulate_fgee <- function(
  GenerateData,
  family,
  N,
  n_i,
  long_corr,
  seed,
  L = 50
) {
  set.seed(seed)
  sim <- suppressMessages(GenerateData(
    Nsubj = N,
    numFunctPoints = L,
    min_visit = n_i,
    max_visit = n_i,
    numLongiPoints = n_i,
    sigma_sq = 10,
    sigma_z11 = 5,
    sigma_z12 = 3,
    sigma_z21 = 2,
    sigma_z22 = 1,
    divFactor = 3,
    corstr.long = long_corr,
    corstr.fn = "ar1",
    rho.fn = 0.75,
    rho.long = 0.75,
    family = family,
    true.cov = FALSE
  ))
  Y <- as.matrix(sim$data$Y)
  dimnames(Y) <- NULL
  d <- data.frame(ID = sim$data$subjID, sim$data$Cov, time = sim$time)
  d$Y <- I(Y)
  list(data = d, s = sim$data$funcArg, beta = t(sim$beta)) # beta: 3 x L
}

# Their initial fit (fastFGEE .fgee_fit_internal / Appendix A.9): bam, fREML, discrete,
# ps basis with knots + 1 = 16 functions and m = c(2, 1) for yindex and intercept.
# reml = TRUE: the same model and bases fitted as in the study (gam, REML).
fit_initial <- function(sim, family, reml = FALSE) {
  # do.call with evaluated arguments: pffr re-evaluates parts of its call elsewhere,
  # where local variables such as `reml` are not visible.
  do.call(
    refund::pffr,
    list(
      formula = Y ~ X1 + X2,
      yind = sim$s,
      data = sim$data,
      family = family,
      algorithm = if (reml) "gam" else "bam",
      method = if (reml) "REML" else "fREML",
      discrete = !reml,
      bs.yindex = list(bs = "ps", k = 16, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 16, m = c(2, 1)),
      sandwich = "none"
    )
  )
}

# Evaluation matrices for beta_0 (full intercept), beta_1, beta_2 at the grid s.
coef_matrices <- function(fit, s) {
  p <- length(fit$coefficients)
  one <- \(sm, by = NULL) {
    grid <- data.frame(s)
    names(grid) <- sm$term
    if (sm$by != "NA") grid[[sm$by]] <- 1
    L <- matrix(0, length(s), p)
    L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(sm, grid)
    L
  }
  bys <- vapply(fit$smooth, \(sm) sm$by, "")
  L0 <- one(fit$smooth[[which(bys == "NA")]])
  L0[, names(fit$coefficients) == "(Intercept)"] <- 1
  list(
    beta0 = L0,
    beta1 = one(fit$smooth[[which(bys == "X1")]]),
    beta2 = one(fit$smooth[[which(bys == "X2")]])
  )
}

se_of <- \(L, V) sqrt(pmax(rowSums((L %*% V) * L), 0))

# fastFGEE's pffr-sandwich baseline (paper: "pffr (z)"/"pffr (Wild)"; fgee(gee.fit =
# FALSE)): reachable in 0.2.2 only through the internal entry point.
fgee_baseline <- function(fit, sim, family, joint_ci = "wild", boot = 3000) {
  fastFGEE:::.fgee_fit_internal(
    formula = Y ~ X1 + X2,
    data = sim$data,
    cluster = "ID",
    family = family,
    corr_fn = "independent",
    corr_long = "independent",
    time = "time",
    pffr.mod = fit,
    gee.fit = FALSE,
    joint.CI = joint_ci,
    boot.samps = boot,
    var.type = "sandwich",
    bs = "ps",
    knots = 15
  )
}

# fastFGEE's proposed one-step estimator (public interface), working correlation
# corr_fn x corr_long.
fgee_onestep <- function(fit, sim, family, corr_fn, corr_long, boot = 3000) {
  fastFGEE::fgee(
    formula = Y ~ X1 + X2,
    data = sim$data,
    cluster = "ID",
    family = family,
    corr_fn = corr_fn,
    corr_long = corr_long,
    time = "time",
    pffr.mod = fit,
    joint.CI = "wild",
    boot.samps = boot,
    bs = "ps",
    knots = 15,
    verbose.tuning = FALSE
  )
}

# The scale mgcv's smoothing parameters refer to: mgcv penalises the deviance for
# known-scale families and the (unscaled) RSS for the Gaussian, so the penalty on
# the score scale is S_lambda / scale with scale = sig2 (1 for Poisson/binomial).
penalty_scale <- function(fit) if (fit$scale.estimated) fit$sig2 else 1
