# Satterthwaite / Bell-McCaffrey df for the CL2 intervals (core cells):
# shared definitions. Sourced after scripts/setup.R and ablation/ablation-lib.R, whose
# cells (ablation_cells), rep-major tasks (ablation_tasks), CRN columns and interval
# scoring (interval_scores) are reused unchanged.

# Per-point df of the CL2-exact sampling quadratic form for the study's contrasts --------

# refund 79a346fb's coef.pffr(crit = "satterthwaite", df_gram = ...) computes
# pffr_df_context(object, "cl2", cl2_adjustment = <the covariance's adjustment>) and
# pffr_df_from_context(ctx, Xp) with Xp the full-coefficient-space contrasts of its own
# grid. Here the same two internals are applied to the study's evaluation matrices L.
# The context's influence core is the one pffr_vcov(sandwich = "cl2") builds its
# covariance from (same cache key), which is checked: the frequentist CL2 covariance
# rebuilt from the context's core must equal pffr_vcov(..., freq = TRUE).
df_context <- function(fit, V_freq) {
  ctx_full <- refund:::pffr_df_context(
    fit,
    "cl2",
    cl2_adjustment = "exact",
    df_gram = "full"
  )
  if (!isTRUE(ctx_full$ok)) stop("refund returned no df context")
  if (!identical(ctx_full$core$adjustment, "exact")) {
    stop("df context does not use the exact CL2 adjustment")
  }
  V_core <- refund:::pffr_influence_vcov(ctx_full$core, freq = TRUE)
  list(
    ctx = ctx_full,
    G = ctx_full$G,
    check = max(abs(V_core - V_freq)) / max(abs(V_freq))
  )
}

# df columns of one estimand: full (default) and diagonal Gram, plus the context's
# expected sampling variance (working model) of the full Gram.
df_columns <- function(ctx, L) {
  full <- refund:::pffr_influence_df(ctx$core, L, df_gram = "full")
  diag <- refund:::pffr_influence_df(ctx$core, L, df_gram = "diagonal")
  list(
    df_full = full$df,
    df_diag = diag$df,
    esv = full$expected_sampling_variance
  )
}

# Critical values ----------------------------------------------------------------------

# z, t_{G-1}, per-point Satterthwaite (full Gram = refund's default; diagonal Gram), and
# satt_b2: the full-Gram df rescaled to the Bayesian covariance. refund's df describes
# the sampling quadratic form of the frequentist CL2 variance; the paper's Bayesian CL2
# variance adds the fixed term Vp - Ve. Satterthwaite for V_B = V_S + const gives
# nu_B = nu_S * (E V_B / E V_S)^2, here with the plug-in ratio (se_b / se_f)^4.
CRITS <- c("z", "t", "satt", "satt_diag", "satt_b2")

critical_values <- function(crit, G, rec) {
  q <- 1 - (1 - LEVEL) / 2
  df <- switch(
    crit,
    z = Inf,
    t = G - 1,
    satt = rec$df_full,
    satt_diag = rec$df_diag,
    satt_b2 = rec$df_full * (rec$se_cl2e_b / rec$se_cl2e_f)^4,
    stop("unknown critical value: ", crit)
  )
  # As coef.pffr: qt at max(df, 1); undefined df (NA) gives missing limits there.
  # Here missing df are counted and scored with z (see n_df_na).
  out <- ifelse(is.finite(df), qt(q, pmax(df, 1)), qnorm(q))
  if (length(out) > 1) dim(out) <- dim(rec$df_full)
  out
}

# REML + CL2 records of one replicate (moved from run-satt.R, 2026-10-03) ---------------

cl2_cov <- function(fit, freq) {
  V <- refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    freq = freq,
    cl2_adjustment = "exact"
  )
  if (!identical(attr(V, "cl2_adjustment"), "exact")) {
    stop("refund did not use the exact CL2 adjustment")
  }
  V
}

satt_records <- function(fit, sim, cell) {
  V_model <- if (!is.null(fit$Vc)) fit$Vc else fit$Vp
  V_b <- cl2_cov(fit, FALSE)
  V_f <- cl2_cov(fit, TRUE)
  start <- proc.time()[["elapsed"]]
  dc <- df_context(fit, V_f)
  se <- \(L, V) sqrt(rowSums((L %*% V) * L))
  records <- lapply(setNames(nm = model_estimands(cell$model)), \(estimand) {
    L <- evaluation_matrix(fit, sim, estimand)
    df <- df_columns(dc$ctx, L)
    cbind(
      truth = sim$truth[[estimand]],
      est_REML = drop(L %*% fit$coefficients),
      se_model = se(L, V_model),
      se_cl2e_b = se(L, V_b),
      se_cl2e_f = se(L, V_f),
      df_full = df$df_full,
      df_diag = df$df_diag,
      esv = df$esv
    )
  })
  list(
    records = records,
    G = dc$G,
    core_check = dc$check,
    df_time = proc.time()[["elapsed"]] - start
  )
}
