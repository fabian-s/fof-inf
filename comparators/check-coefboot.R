# Equivalence check: the BOOT arm (comparators-lib.R: curve_bootstrap) against refund's
# pffr_coefboot(method = "resample", type = "percent") on one synthetic replicate, same
# seed and B. pffr_coefboot needs a fit whose stored call it can re-evaluate, so the
# reference fit is made with every pffr argument evaluated (do.call). Compared: the
# 95% percentile limits of gamma(t) on EVAL_T (coef.pffr's n1 = 31 grid on [0, 1]).
# Run from the repository root: Rscript comparators/check-coefboot.R [cell] [B]

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")
suppressPackageStartupMessages(library(refund))

args <- commandArgs(trailingOnly = TRUE)
cell_id <- if (length(args) >= 1) as.integer(args[1]) else 11L
B <- if (length(args) >= 2) as.integer(args[2]) else 19L
seed <- 1L

cell <- comparator_syn_cells()
cell <- cell[cell$cell == cell_id, ]
sim <- simulate(cell, 1)
fit <- fit_model(sim, cell, "REML")
Ls <- reml_matrices(fit, sim, cell)

# reference fit with an evaluable call (same model as fit_pffr)
bases <- cell_bases(cell)
xind <- sim$xind
ff_basis <- bases$ff
ff <- refund::ff
ref <- do.call(
  refund::pffr,
  list(
    formula = Y ~ ff(X, xind = xind, splinepars = ff_basis) + z,
    data = sim$data,
    yind = sim$yind,
    bs.yindex = bases$t,
    bs.int = bases$t,
    family = gaussian(),
    sandwich = "none",
    algorithm = "gam",
    method = "REML"
  )
)
stopifnot(max(abs(ref$coefficients - fit$coefficients)) < 1e-8)

set.seed(seed)
cb <- refund::pffr_coefboot(
  ref,
  n1 = length(EVAL_T),
  B = B,
  conf = LEVEL,
  type = "percent",
  method = "resample",
  showProgress = FALSE
)
gam_terms <- cb$smterms
idx <- grep("z", names(gam_terms))
cat("pffr_coefboot terms:", names(gam_terms), "\n")
g <- gam_terms[[idx]]$coef
print(head(g))

ours <- curve_bootstrap(sim, cell, Ls, B, seed)
lo <- ours$values$gamma$lo
hi <- ours$values$gamma$hi
cat("columns:", names(g), "\n")
lower_ref <- g[[ncol(g) - 1]] # coef.pffr layout: ..., value, lower, upper
upper_ref <- g[[ncol(g)]]
grid_ok <- isTRUE(all.equal(g[[1]], EVAL_T))
cat(sprintf(
  "grid matches EVAL_T: %s; max |lower diff| %.2e; max |upper diff| %.2e; failed %d\n",
  grid_ok,
  max(abs(lower_ref - lo)),
  max(abs(upper_ref - hi)),
  ours$n_failed
))
cat("info:", ours$row$info, "\n")
