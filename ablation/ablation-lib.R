# Covariance ablation (core cells): shared definitions. Sourced after
# scripts/setup.R (which provides the study's design, data generation and fit code).

# Cells, tasks, arms ------------------------------------------------------------------

# The 18 core cells at the default signal level (M1, D = 61, default bases).
ablation_cells <- function(cells = design_table()) {
  core <- cells[cells$block == "core" & cells$signal == "mid", ]
  stopifnot(
    nrow(core) == 18,
    all(core$model == "M1"),
    all(core$D == 61),
    all(core$basis == "default")
  )
  core
}

# Rep-major task table over the core cells; tasks are looked up by id.
ablation_tasks <- function(cells = ablation_cells(), reps = seq_len(N_REP)) {
  tasks <- tidyr::expand_grid(rep = reps, cell = cells$cell)
  tasks$task <- seq_len(nrow(tasks))
  dplyr::relocate(tasks, task)
}

# Covariance arms. fit: "REML" (smoothing parameters selected per replicate) or "FIX"
# (fixed at the cell's geometric-mean REML values). sandwich "none" = model-based
# (Vc where available, else Vp, as in R/fit.R). "cluster" = refund's cluster sandwich
# with its default dof_correction = "none" (no small-sample factor); "cl2" with the
# diagonal-leverage ("shortcut") or full-block ("exact") adjustment. freq = FALSE is the
# Bayesian form (+ Vp - Ve), freq = TRUE the frequentist form.
COV_ARMS <- tibble::tribble(
  ~arm,
  ~fit,
  ~sandwich,
  ~adjustment,
  ~freq,
  "model",
  "REML",
  "none",
  NA,
  FALSE,
  "cr1_b",
  "REML",
  "cluster",
  NA,
  FALSE,
  "cr1_f",
  "REML",
  "cluster",
  NA,
  TRUE,
  "cl2s_b",
  "REML",
  "cl2",
  "shortcut",
  FALSE,
  "cl2s_f",
  "REML",
  "cl2",
  "shortcut",
  TRUE,
  "cl2e_b",
  "REML",
  "cl2",
  "exact",
  FALSE,
  "cl2e_f",
  "REML",
  "cl2",
  "exact",
  TRUE,
  "fix_model",
  "FIX",
  "none",
  NA,
  FALSE,
  "fix_cl2e_b",
  "FIX",
  "cl2",
  "exact",
  FALSE
)

# Stored study record columns that the REML fit must reproduce (common random numbers):
# ablation column = study column.
CRN_COLUMNS <- c(
  truth = "truth",
  est_REML = "est_REML",
  se_model = "se_model_REML",
  se_cl2e_b = "se_cl2_REML",
  se_cl2e_f = "se_cl2f_REML"
)

# Covariances ---------------------------------------------------------------------------

# Covariance of one arm from a sandwich = "none" fit, via the calls R/fit.R uses.
arm_covariance <- function(fit, spec) {
  if (spec$sandwich == "none") {
    return(if (!is.null(fit$Vc)) fit$Vc else fit$Vp)
  }
  V <- refund:::pffr_vcov(
    fit,
    sandwich = spec$sandwich,
    freq = spec$freq,
    cl2_adjustment = if (is.na(spec$adjustment)) NULL else spec$adjustment
  )
  if (
    spec$sandwich == "cl2" &&
      !identical(attr(V, "cl2_adjustment"), spec$adjustment)
  ) {
    stop("refund did not use the ", spec$adjustment, " CL2 adjustment")
  }
  V
}

# Fixed smoothing parameters of a cell, named and in the order of `sp_names` (fit$sp).
fixed_sp <- function(sp_table, cell_id, sp_names) {
  rows <- sp_table[sp_table$cell == cell_id, ]
  sp <- setNames(exp(rows$mean_log_sp), rows$penalty)
  if (!setequal(names(sp), sp_names) || length(sp) != length(sp_names)) {
    stop("fixed sp penalties do not match fit$sp for cell ", cell_id)
  }
  sp[sp_names]
}

# Scoring --------------------------------------------------------------------------------

# Point-level scores of one interval (link-scale interval; for `mean` the estimate and
# both endpoints are transformed to the response scale, as in R/summarise.R).
interval_scores <- function(est, se, truth, crit, estimand, linkinv) {
  lower <- est - crit * se
  upper <- est + crit * se
  if (estimand == "mean") {
    lower <- linkinv(lower)
    upper <- linkinv(upper)
  }
  alpha <- 1 - LEVEL
  list(
    covered = lower <= truth & truth <= upper,
    width = upper - lower,
    score = (upper - lower) +
      2 /
        alpha *
        ((lower - truth) * (truth < lower) + (truth - upper) * (truth > upper))
  )
}
