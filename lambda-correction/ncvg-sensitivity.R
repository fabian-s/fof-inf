# Threshold sensitivity of the boundary rule + sign check of the cross term (GLM NCV), local:
#   Rscript lambda-correction/ncvg-sensitivity.R <cell> <reps>   (from the repository root; writes validation/)
source("scripts/setup.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")
source("lambda-correction/ncvg-lib.R")
a <- commandArgs(TRUE)
cid <- as.integer(a[1])
reps <- eval(parse(text = a[2]))
cell <- design_table()
cell <- cell[cell$cell == cid, ]
out <- list()
for (r in reps) {
  sim <- simulate(cell, r)
  fit <- fit_model(sim, cell, "NCV")
  P <- ncvg_pieces(fit)
  hs <- ncvg_hessian_scores(P)
  core <- refund:::pffr_influence(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact"
  )
  L <- evaluation_matrix(fit, sim, "beta")
  se <- \(V) sqrt(pmax(rowSums((L %*% V) * L), 0))
  se0 <- se(core$correction * tcrossprod(core$K) + core$B2)
  # sign check: fixed-rho curve influence dbeta/dw_h = H^-1 X_h' u_h vs CL2 influence columns
  st <- ncvg_state(P)
  dbw <- sapply(
    P$idx,
    \(ix) drop(st$Hinv %*% crossprod(P$X[ix, , drop = FALSE], st$u[ix]))
  )
  sgn <- median(sapply(seq_along(P$idx), \(h) cor(dbw[, h], core$K[, h])))
  xj <- sqrt(colSums((P$X %*% P$J)^2))
  xj <- xj / max(xj)
  dh <- diag(hs$H) / max(diag(hs$H))
  rules <- list(
    c1e2 = dh > 1e-2,
    c1e3 = dh > 1e-3,
    c1e4 = dh > 1e-4,
    xj1e2 = xj > 1e-2 & diag(hs$H) > 0,
    xj1e3 = xj > 1e-3 & diag(hs$H) > 0
  )
  for (nm in names(rules)) {
    free <- rules[[nm]]
    D <- matrix(0, length(P$groups), length(P$sp))
    if (any(free))
      D[, free] <- ncv_ij_influence(
        P,
        hs$H[free, free, drop = FALSE],
        hs$s[, free, drop = FALSE]
      )
    R <- rho_part_linear(P, D)
    out[[length(out) + 1]] <- data.frame(
      cell = cid,
      rep = r,
      rule = nm,
      n_fixed = sum(!free),
      infl_ij = mean(se(add_rho_part(core, R, "ij")) / se0),
      infl_ijx = mean(se(add_rho_part(core, R, "ij", cross = TRUE)) / se0),
      sign_cor = sgn,
      min_rel_curv = min(dh),
      xj = paste(signif(xj, 2), collapse = "/"),
      curv = paste(signif(dh, 2), collapse = "/")
    )
  }
  cat(r, "\n")
}
saveRDS(
  do.call(rbind, out),
  sprintf("lambda-correction/validation/ncvg-sens-%d.rds", cid)
)
