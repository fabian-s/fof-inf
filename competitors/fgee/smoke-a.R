# Smoke test for experiment A: one data set, all covariances on one pffr fit.
# Run from competitors/fgee with the BLAS pinned.
.libPaths(c(path.expand("~/R/fastfgee-lib"), .libPaths()))
suppressPackageStartupMessages({
  library(refund)
  library(fastFGEE)
})
source("fgee-dgp.R")
GenerateData <- source_upstream(".")

args <- commandArgs(TRUE)
fam_name <- args[1] %||% "gaussian"
N <- as.integer(args[2] %||% 25)
family <- switch(
  fam_name,
  gaussian = gaussian(),
  poisson = poisson(),
  binomial = binomial()
)

sim <- simulate_fgee(
  GenerateData,
  fam_name,
  N = N,
  n_i = 5,
  long_corr = "exchangeable",
  seed = 1
)
t0 <- proc.time()[[3]]
fit <- fit_initial(sim, family)
cat(
  "pffr fit",
  round(proc.time()[[3]] - t0, 1),
  "s; p =",
  length(coef(fit)),
  "sig2 =",
  fit$sig2,
  "\n"
)

base <- fgee_baseline(fit, sim, family, boot = 500)
cat(
  "baseline coefficients == pffr coefficients:",
  isTRUE(all.equal(unname(as.numeric(base$beta)), unname(fit$coefficients))),
  "\n"
)

Ls <- coef_matrices(fit, sim$s)
ID <- sim$data$ID
V <- list(
  model = fit$Vp,
  fgee_base = base$vb,
  refund_cluster = refund:::pffr_vcov(fit, sandwich = "cluster", cluster = ID),
  refund_cl2 = refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    cluster = ID
  )
)
# Baseline with the penalty put on the averaged score scale: P / (N * scale).
w <- base$working0
meat <- w$dd_sum - w$N * tcrossprod(w$d_bar)
A_fix <- w$W_bar + base$pen.mat / (w$N * penalty_scale(fit))
Ai <- solve(A_fix)
V$fgee_base_fixed <- Ai %*% meat %*% Ai / (w$N * (w$N - 1))
# Check: the baseline's own bread reproduces vb
A_raw <- w$W_bar + base$pen.mat
Ar <- solve(A_raw)
cat(
  "reconstructed baseline == vb:",
  isTRUE(all.equal(
    Ar %*% meat %*% Ar / (w$N * (w$N - 1)),
    base$vb,
    tolerance = 1e-6
  )),
  "\n"
)

for (r in names(Ls)) {
  L <- Ls[[r]]
  ses <- sapply(V, \(v) se_of(L, v))
  cat("\n", r, ": median SE ratio to refund_cluster\n")
  print(round(apply(ses / ses[, "refund_cluster"], 2, median), 3))
  est <- drop(L %*% fit$coefficients)
  truth <- sim$beta[match(r, names(Ls)), ]
  cat("  pointwise coverage (1 data set):\n")
  print(round(colMeans(abs(est - truth) <= 1.96 * ses), 2))
}
# Their CI centres/widths vs ours (pointwise, from their CI object)
ci <- base$crit$ci %||% base$model$crit$ci
str(ci, max.level = 1)

# Mapping checks: their SE / fit curves vs L %*% (vb, coefficients)
for (k in seq_along(Ls)) {
  their_se <- ci$se[[k]]
  their_fit <- ci$fit[[k]]
  cat(
    names(Ls)[k],
    ": max |their_se - ours| / ours =",
    signif(
      max(abs(their_se - se_of(Ls[[k]], base$vb)) / se_of(Ls[[k]], base$vb)),
      3
    ),
    "; max |their_fit - ours| =",
    signif(max(abs(their_fit - drop(Ls[[k]] %*% fit$coefficients))), 3),
    "\n"
  )
}
cat(
  "max rel diff reconstructed vs vb:",
  signif(
    max(abs(Ar %*% meat %*% Ar / (w$N * (w$N - 1)) - base$vb)) /
      max(abs(base$vb)),
    3
  ),
  "\n"
)
str(ci$df)
str(ci$info)
