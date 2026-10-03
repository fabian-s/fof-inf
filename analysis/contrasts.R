# Uncertainty for pooled and cross-cell quantities. All cells share their random draws
# per replicate (common random numbers), so a replicate bootstrap resamples replicate
# ids jointly across cells and arms: every statistic below is computed from task-level
# scores weighted by w[rep], and its interval from the resampled weights.

# Resampling ------------------------------------------------------------------------------

METRICS <- c("coverage", "width", "score", "sq_error", "deviance")

# Task-level metrics as (cell x estimand x arm) by replicate matrices, so that cell means
# under replicate weights are matrix-vector products.
prepare_resampling <- function(tasks) {
  keys <- dplyr::distinct(tasks, cell, estimand, arm)
  reps <- sort(unique(tasks$rep))
  row <- match(
    paste(tasks$cell, tasks$estimand, tasks$arm),
    paste(keys$cell, keys$estimand, keys$arm)
  )
  col <- match(tasks$rep, reps)
  matrices <- lapply(setNames(METRICS, METRICS), \(m) {
    x <- matrix(NA_real_, nrow(keys), length(reps))
    x[cbind(row, col)] <- tasks[[m]]
    x
  })
  list(keys = keys, reps = reps, matrices = matrices)
}

# Cell means (and the coverage MC SE) under replicate weights w (frequency weights);
# missing replicates of a cell (failed tasks) get weight 0 in that cell.
cell_means <- function(prep, w = rep(1, length(prep$reps))) {
  means <- lapply(prep$matrices, \(x) {
    present <- !is.na(x)
    x[!present] <- 0
    list(
      mean = drop(x %*% w) / drop(present %*% w),
      sq = drop(x^2 %*% w),
      n = drop(present %*% w)
    )
  })
  cov <- means$coverage
  var <- (cov$sq / cov$n - cov$mean^2) * cov$n / (cov$n - 1)
  out <- prep$keys
  for (m in METRICS) out[[m]] <- means[[m]]$mean
  out$coverage_se <- sqrt(pmax(var, 0) / cov$n)
  out
}

# Estimate and percentile interval of stat(cell_means) under the replicate bootstrap.
# `stat` returns a named numeric vector.
boot_stat <- function(tasks, stat, n_boot = N_BOOT, seed = SEED_BOOT) {
  prep <- prepare_resampling(tasks)
  n_rep <- length(prep$reps)
  estimate <- stat(cell_means(prep))
  draws <- withr::with_seed(
    seed,
    replicate(n_boot, {
      stat(cell_means(prep, tabulate(sample.int(n_rep, replace = TRUE), n_rep)))
    })
  )
  draws <- matrix(draws, nrow = length(estimate))
  tibble::tibble(
    name = names(estimate),
    estimate = unname(estimate),
    se = apply(draws, 1, sd, na.rm = TRUE),
    lower = apply(draws, 1, quantile, 0.025, na.rm = TRUE),
    upper = apply(draws, 1, quantile, 0.975, na.rm = TRUE)
  )
}

# A pooled statistic needs every cell of its prescribed set.
require_cells <- function(tasks, required, what) {
  missing <- setdiff(required, unique(tasks$cell))
  if (length(missing)) stop(what, ": no results for cells ", toString(missing))
}

straddles <- function(lower, upper, threshold)
  lower < threshold & upper > threshold

# The recommendation with uncertainty --------------------------------------------

# Calibration errors and interval-score ratio of proposal vs fallback on the pool, with
# replicate-bootstrap intervals; a criterion whose interval straddles its threshold is
# flagged unresolved. The decision itself is recommend() on the point estimates.
recommend_with_uncertainty <- function(
  tasks,
  summaries,
  cells,
  family,
  estimand,
  proposal = PROPOSAL,
  signal = NULL
) {
  pool <- recommendation_pool(cells, family, signal)
  sub <- dplyr::filter(
    tasks,
    cell %in% pool,
    estimand == !!estimand,
    arm %in% c(proposal, FALLBACK)
  )
  if (!setequal(unique(sub$cell), pool)) {
    return(list(
      decision = recommend(
        summaries,
        cells,
        family,
        estimand,
        proposal,
        signal
      ),
      criteria = tibble::tibble(name = NA_character_, estimate = NA_real_)
    ))
  }
  stat <- \(m) {
    get_arm <- \(a) dplyr::arrange(dplyr::filter(m, arm == a), cell)
    p <- get_arm(proposal)
    f <- get_arm(FALLBACK)
    c(
      ce_proposal = calibration_error(p$coverage, p$coverage_se),
      ce_fallback = calibration_error(f$coverage, f$coverage_se),
      score_ratio = exp(mean(log(p$score / f$score)))
    )
  }
  intervals <- boot_stat(sub, stat)
  flags <- with(
    intervals,
    dplyr::case_when(
      name == "score_ratio" ~
        straddles(lower, upper, PRACTICAL_RATIO) |
          straddles(lower, upper, 1 / PRACTICAL_RATIO),
      TRUE ~ straddles(lower, upper, ADEQUATE)
    )
  )
  list(
    decision = recommend(summaries, cells, family, estimand, proposal, signal),
    criteria = dplyr::mutate(intervals, unresolved = flags)
  )
}

# Question 1: does CL2 replace model-based intervals? Calibration error of both on the
# pool (for REML and NCV fits), and the coverage cost of CL2 under iid (core iid cells,
# G = 100, all signal levels).
question1 <- function(tasks, cells, family, estimand) {
  pool <- recommendation_pool(cells, family)
  iid <- cells$cell[
    cells$block == "core" &
      cells$error == "iid" &
      cells$G == 100L &
      cells$family == family
  ]
  arms <- c("reml_model", "reml_cl2", "ncv_model", "ncv_cl2")
  sub <- dplyr::filter(
    tasks,
    cell %in% c(pool, iid),
    estimand == !!estimand,
    arm %in% arms
  )
  stat <- \(m) {
    ce <- sapply(arms, \(a) {
      x <- dplyr::filter(m, cell %in% pool, arm == a)
      calibration_error(x$coverage, x$coverage_se)
    })
    iid_cov <- \(a) mean(dplyr::filter(m, cell %in% iid, arm == a)$coverage)
    c(
      setNames(ce, paste0("ce_", arms)),
      # paired improvement of CL2 over model-based (recommended if >= 2 pp)
      ce_gain_reml = ce[["reml_model"]] - ce[["reml_cl2"]],
      ce_gain_ncv = ce[["ncv_model"]] - ce[["ncv_cl2"]],
      iid_cost_reml = iid_cov("reml_model") - iid_cov("reml_cl2"),
      iid_cost_ncv = iid_cov("ncv_model") - iid_cov("ncv_cl2")
    )
  }
  require_cells(sub, c(pool, iid), "question 1")
  out <- boot_stat(sub, stat)
  gain <- grepl("^ce_gain", out$name)
  cost <- grepl("^iid_cost", out$name)
  dplyr::mutate(
    out,
    unresolved = (gain & straddles(lower, upper, ADEQUATE)) |
      (cost & straddles(lower, upper, ADEQUATE))
  )
}

# Profiles over the core factorial ---------------------------------------------------------

# For one family and estimand: the mean over core cells of a per-cell value, at each
# level of `factor` (averaged over the other core factors), with bootstrap intervals.
# `value` maps the cell means of one cell (all arms) to a number, e.g. an arm's
# coverage or the log MSE ratio of two arms.
factor_profile <- function(tasks, cells, family, estimand, factor, value) {
  core <- cells[cells$block == "core" & cells$family == family, ]
  sub <- dplyr::filter(tasks, cell %in% core$cell, estimand == !!estimand)
  require_cells(sub, core$cell, paste("profile", family, estimand))
  levels <- unique(core[[factor]])
  stat <- \(m) {
    per_cell <- sapply(split(m, m$cell), value)
    level_of <- core[[factor]][match(as.integer(names(per_cell)), core$cell)]
    setNames(sapply(levels, \(l) mean(per_cell[level_of == l])), levels)
  }
  dplyr::mutate(boot_stat(sub, stat), factor = factor, .before = 1)
}

# Per-cell value functions; NA when the cell lacks the arm (e.g. no replicates yet).
one_or_na <- function(x) if (length(x) == 1) x else NA_real_
arm_value <- function(metric, arm) \(m) one_or_na(m[[metric]][m$arm == arm])
log_value <- function(metric, arm) \(m) log(arm_value(metric, arm)(m))
log_ratio <- function(metric, arm_a, arm_b) {
  \(m)
    log(
      one_or_na(m[[metric]][m$arm == arm_a]) /
        one_or_na(m[[metric]][m$arm == arm_b])
    )
}

# Matched cells (extensions vs core) ---------------------------------------------------------

# Paired difference of one per-cell value between two cells (e.g. dense grid vs its
# core cell at D = 61), with a replicate-bootstrap interval.
matched_contrast <- function(tasks, cell_a, cell_b, estimand, value) {
  sub <- dplyr::filter(
    tasks,
    cell %in% c(cell_a, cell_b),
    estimand == !!estimand
  )
  stat <- \(m)
    c(diff = value(m[m$cell == cell_a, ]) - value(m[m$cell == cell_b, ]))
  boot_stat(sub, stat)
}

# The core cell an extension cell is compared with: same family and error, and the core
# defaults for everything the extension changed (basis_size cells: the default-basis
# core cell with the same family, error, G and signal).
matching_core_cell <- function(cells, cell_id) {
  x <- cells[cells$cell == cell_id, ]
  hit <- cells$block == "core" &
    cells$family == x$family &
    cells$G == x$G &
    cells$signal == x$signal &
    cells$error == x$error
  # warp cells: the iid cell with the same family, G, signal and truth (core or rough
  # truth block); none exists for wiggly truths at G = 40
  if (x$block == "warp") {
    hit <- cells$block %in%
      c("core", "rough_truth") &
      cells$family == x$family &
      cells$G == x$G &
      cells$signal == x$signal &
      cells$truth == x$truth &
      cells$error == "iid"
  }
  cells$cell[hit]
}
