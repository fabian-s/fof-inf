# Realistic smooth error process: residual functional principal components of the DTI
# application (baseline visit, rcst ~ ff(cca)), estimated once and frozen in
# data/fpc-kernel.rds. Run from the repository root.
suppressPackageStartupMessages(library(refund))
source("R/constants.R")

dti <- refund::DTI
base <- dti$visit == 1 & complete.cases(dti$cca) & complete.cases(dti$rcst)
Y <- dti$rcst[base, ]
X <- dti$cca[base, ]
t_grid <- seq(0, 1, length.out = ncol(Y))
s_grid <- seq(0, 1, length.out = ncol(X))

fit <- pffr(
  Y ~ ff(X, xind = s_grid, splinepars = BASIS_FF),
  yind = t_grid,
  bs.yindex = BASIS_T,
  bs.int = BASIS_T,
  sandwich = "none"
)
res <- matrix(residuals(fit), nrow = nrow(Y))

# Smoothed residual covariance (FACE) -> eigenfunctions, eigenvalues, white-noise variance.
fpca <- fpca.face(res, argvals = t_grid, pve = 0.99, var = TRUE)
kernel <- list(
  grid = t_grid,
  efunctions = fpca$efunctions / sqrt(diff(t_grid[1:2])), # L2[0,1]-normalized
  evalues = fpca$evalues * diff(t_grid[1:2]),
  sigma2 = fpca$sigma2,
  source = "refund::DTI, visit 1, rcst ~ ff(cca), fpca.face(pve = 0.99)",
  n_curves = nrow(Y)
)
saveRDS(kernel, "data/fpc-kernel.rds")
share <- kernel$evalues / sum(kernel$evalues)
cat(
  "curves",
  nrow(Y),
  "components",
  length(share),
  "\nshares",
  round(share, 3),
  "\nnoise share of total variance",
  round(kernel$sigma2 / (kernel$sigma2 + sum(kernel$evalues)), 3),
  "\n"
)
