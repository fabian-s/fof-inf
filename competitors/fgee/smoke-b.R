# Smoke test for experiment B: fastFGEE one-step on the study's own data sets
# (M1: ff(X) + z, one curve per cluster). Run from the repository root with the BLAS pinned.
.libPaths(c(path.expand("~/R/fastfgee-lib"), .libPaths()))
source("scripts/setup.R")
suppressPackageStartupMessages(library(fastFGEE))

args <- commandArgs(TRUE)
fam <- args[1] %||% "gaussian"
model <- args[2] %||% "M1"
cells <- design_table()
cell <- cells[
  cells$block == ifelse(model == "M1", "core", "term_type") &
    cells$family == fam &
    cells$error == "fpc" &
    cells$G == 100L &
    cells$signal == "mid",
][1, ]
print(as.data.frame(cell))
sim <- simulate(cell, 1L)
fit <- fit_model(sim, cell, "REML")
cat("REML fit ok, p =", length(fit$coefficients), "\n")
d <- sim$data
d$ID <- seq_len(nrow(d))
xind <- sim$xind
ff <- refund::ff
form <- if (model == "M1") {
  Y ~ ff(X, xind = xind, splinepars = BASIS_FF) + z
} else {
  Y ~ s(x) + z
}
for (cf in c("independent", "ar1", "fpca")) {
  t0 <- proc.time()[[3]]
  res <- tryCatch(
    fastFGEE::fgee(
      formula = form,
      data = d,
      cluster = "ID",
      family = FAMILIES[[cell$fit_family]]$fit(),
      corr_fn = cf,
      corr_long = "independent",
      pffr.mod = fit,
      joint.CI = "wild",
      boot.samps = 500,
      verbose.tuning = FALSE
    ),
    error = \(e) e
  )
  if (inherits(res, "error")) {
    cat(cf, "ERROR:", conditionMessage(res), "\n")
    next
  }
  cat(
    cf,
    "ok in",
    round(proc.time()[[3]] - t0, 1),
    "s; length(beta) =",
    length(res$beta),
    "; lambda =",
    signif(res$lambda, 3),
    "\n"
  )
  for (est in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, est)
    e_reml <- drop(L %*% fit$coefficients)
    e_fgee <- drop(L %*% res$beta)
    se_f <- sqrt(rowSums((L %*% res$vb) * L))
    cat(sprintf(
      "  %-6s cor(est) = %.3f  median SE fgee = %.3g\n",
      est,
      cor(e_reml, e_fgee),
      median(se_f)
    ))
  }
}
