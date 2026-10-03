# Algebra check: Gaussian, fixed design, fixed penalties, known error
# covariance. At fixed penalties the pffr estimate is linear, beta_hat = B X'y with
# B = (X'X + S)^-1, so its mean is B X' mu and its covariance B X'(I_G x Sigma)X B.
# The Monte Carlo mean and covariance must match; the MC mean of CL2 and of the
# model-based covariance are compared with the truth (reported, not gated).
# Run from the repository root: Rscript scripts/check-algebra.R [n_mc]
suppressPackageStartupMessages(library(refund))
source("scripts/setup.R")
args <- commandArgs(trailingOnly = TRUE)
n_mc <- if (length(args)) as.integer(args[1]) else 300L
stopifnot("n_mc must be a positive integer" = !is.na(n_mc) && n_mc > 0)

cells <- design_table()
cell <- as.list(cells[
  cells$block == "core" &
    cells$family == "gaussian" &
    cells$error == "fpc" &
    cells$G == 40L &
    cells$signal == "mid",
])
sim <- simulate(cell, 1)
fit0 <- fit_model(sim, cell, "REML")

# Fixed pieces: mean, error covariance, design, penalty, and the beta evaluation matrix.
truth <- make_truth(cell$model, cell$truth)
cov <- list(X = unclass(sim$data$X), z = sim$data$z)
mu <- as.vector(t(linear_predictor(truth, cov, sim$yind, 1, 0))) # curve-major stacking
sigma <- signal_scale(cell)$sigma
t_index <- T_SUBSAMPLE[[as.character(cell$D)]]
Sigma <- sigma^2 * correlation_kernel("fpc", T_DENSE)[t_index, t_index]
X <- predict(fit0, type = "lpmatrix", reformat = FALSE)
S <- matrix(0, ncol(X), ncol(X))
j <- 0L
for (sm in fit0$smooth) {
  ii <- sm$first.para:sm$last.para
  for (P in sm$S) {
    j <- j + 1L
    S[ii, ii] <- S[ii, ii] + fit0$sp[j] * P
  }
}
B <- solve(crossprod(X) + S)
XSX <- Reduce(
  `+`,
  lapply(seq_len(cell$G), \(g) {
    rows <- (g - 1) * length(t_index) + seq_along(t_index)
    crossprod(X[rows, ], Sigma %*% X[rows, ])
  })
)
V_true <- B %*% XSX %*% B
L <- evaluation_matrix(fit0, sim, "beta")

# Monte Carlo over error draws only (design and penalties fixed).
chol_S <- chol(Sigma)
draws <- withr::with_seed(
  99,
  lapply(seq_len(n_mc), \(r) {
    sim$data$Y <- I(
      matrix(mu, cell$G, byrow = TRUE) +
        matrix(rnorm(cell$G * length(t_index)), cell$G) %*% chol_S
    )
    fit <- pffr(
      Y ~ ff(X, xind = S_GRID, splinepars = BASIS_FF) + z,
      data = sim$data,
      yind = sim$yind,
      bs.yindex = BASIS_T,
      bs.int = BASIS_T,
      sandwich = "none",
      sp = fit0$sp
    )
    list(
      coef = fit$coefficients,
      v_cl2 = refund:::pffr_vcov(
        fit,
        sandwich = "cl2",
        cl2_adjustment = "exact"
      ),
      v_model = refund:::pffr_vcov(fit, sandwich = "none")
    )
  })
)

coefs <- do.call(rbind, lapply(draws, `[[`, "coef"))
mean_mc <- colMeans(coefs)
mean_exact <- drop(B %*% crossprod(X, mu))
se_mean <- sqrt(diag(V_true) / n_mc)
var_of <- \(V) rowSums((L %*% V) * L) # variance of beta(s, t) at the eval points
avg_v <- \(name) Reduce(`+`, lapply(draws, `[[`, name)) / n_mc
ratio <- \(v) summary(sqrt(v / var_of(V_true)))

cat(sprintf(
  "n_mc = %d, G = %d, D = %d, cell %d\n",
  n_mc,
  cell$G,
  cell$D,
  cell$cell
))
cat(sprintf(
  "mean: max |MC - exact| / MC-SE = %.2f (should be < ~4)\n",
  max(abs(mean_mc - mean_exact) / se_mean)
))
cat("SE ratio for beta(s,t), empirical / exact:\n")
print(ratio(var_of(cov(coefs))))
cat("SE ratio, mean CL2 / exact:\n")
print(ratio(var_of(avg_v("v_cl2"))))
cat("SE ratio, mean model-based / exact:\n")
print(ratio(var_of(avg_v("v_model"))))
