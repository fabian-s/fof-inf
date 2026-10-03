# From saved task records to cell summaries, paired contrasts and the §5 rule.

# Interval types ----------------------------------------------------------------------

# Each interval type: its estimate and SE columns in a task record (link scale). The
# bias-aware interval adds the NCV-REML difference in quadrature to the NCV CL2 variance.
# ncv_cl2f_bias uses the frequentist CL2 (without refund's Bayesian Vp - Ve term), so
# smoothing bias is not counted twice. ncv_point has no interval; it enters
# only the MSE comparisons.
ARMS <- data.frame(
  arm = c(
    "reml_model",
    "reml_cl2",
    "ncv_model",
    "ncv_cl2",
    "ncv_cl2_bias",
    "ncv_cl2f_bias",
    "ar1_model",
    "ncv_point"
  ),
  est = c(
    "est_REML",
    "est_REML",
    "est_NCV",
    "est_NCV",
    "est_NCV",
    "est_NCV",
    "est_AR1",
    "est_NCV_point"
  ),
  se = c(
    "se_model_REML",
    "se_cl2_REML",
    "se_model_NCV",
    "se_cl2_NCV",
    "se_cl2_NCV",
    "se_cl2f_NCV",
    "se_model_AR1",
    NA
  ),
  bias_aware = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
)

# Plasmode block arm: CL2 clustered on the blocks ("b"),
# or on the alternative partition ("a"), for the REML, NCV and block-NCV fits (NCV_block
# "ncvb", NCV_alt "ncva"); bias-aware intervals use the REML fit as reference, as above.
# Arms whose columns a record lacks are skipped.
ARMS <- rbind(
  ARMS,
  data.frame(
    arm = c(
      "reml_cl2b",
      "ncv_cl2b",
      "ncv_cl2b_bias",
      "ncvb_model",
      "ncvb_cl2",
      "ncvb_cl2b",
      "ncvb_cl2_bias",
      "ncvb_cl2b_bias",
      "reml_cl2a",
      "ncva_model",
      "ncva_cl2",
      "ncva_cl2a",
      "ncva_cl2_bias",
      "ncva_cl2a_bias"
    ),
    est = c(
      "est_REML",
      "est_NCV",
      "est_NCV",
      rep("est_NCV_block", 5),
      "est_REML",
      rep("est_NCV_alt", 5)
    ),
    se = c(
      "se_cl2b_REML",
      "se_cl2b_NCV",
      "se_cl2b_NCV",
      "se_model_NCV_block",
      "se_cl2_NCV_block",
      "se_cl2b_NCV_block",
      "se_cl2_NCV_block",
      "se_cl2b_NCV_block",
      "se_cl2a_REML",
      "se_model_NCV_alt",
      "se_cl2_NCV_alt",
      "se_cl2a_NCV_alt",
      "se_cl2_NCV_alt",
      "se_cl2a_NCV_alt"
    ),
    bias_aware = c(
      FALSE,
      FALSE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      TRUE,
      TRUE
    )
  )
)

# Point-level scores of one interval type: intervals are built on the link scale; for
# `mean` the estimate and both endpoints are transformed to the response scale.
point_scores <- function(record, arm, estimand, linkinv) {
  spec <- ARMS[ARMS$arm == arm, ]
  est <- record[, spec$est]
  se <- if (is.na(spec$se)) NA_real_ else record[, spec$se]
  if (spec$bias_aware) se <- sqrt(se^2 + (est - record[, "est_REML"])^2)
  z <- qnorm(1 - (1 - LEVEL) / 2)
  lower <- est - z * se
  upper <- est + z * se
  if (estimand == "mean") {
    lower <- linkinv(lower)
    upper <- linkinv(upper)
    est <- linkinv(est)
  }
  truth <- record[, "truth"]
  alpha <- 1 - LEVEL
  list(
    covered = lower <= truth & truth <= upper,
    width = upper - lower,
    score = (upper - lower) +
      2 /
        alpha *
        ((lower - truth) * (truth < lower) + (truth - upper) * (truth > upper)),
    error = est - truth
  )
}

# Task-level summaries -------------------------------------------------------------------

# One row per (estimand, arm) of one task: grid averages, plus the point vectors
# (kept as list columns) needed for pointwise coverage and bias/variance.
# Expected excess deviance of a predicted mean mu_hat when the truth is mu, E_Y d(Y,
# mu_hat) - E_Y d(Y, mu) = d(mu, mu_hat) for exponential families; NA otherwise.
excess_deviance <- function(mu, mu_hat, family) {
  if (!family %in% c("gaussian", "poisson", "binomial")) return(NA_real_)
  mean(FAMILIES[[family]]$fit()$dev.resids(mu, mu_hat, 1))
}

# Estimands that are scored in a cell. Under misregistration with a non-identity link,
# the observed-clock mean is not a GLM mean, so only the response-scale mean (whose
# target is exact) is scored there.
scored_estimands <- function(cell, estimands) {
  if (cell$error == "warp" && cell$fit_family != "gaussian") return("mean")
  estimands
}

score_task_records <- function(result, cell) {
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  rows <- list()
  for (estimand in scored_estimands(cell, names(result$records))) {
    record <- result$records[[estimand]]
    available <- ARMS$est %in%
      colnames(record) &
      (is.na(ARMS$se) | ARMS$se %in% colnames(record))
    for (arm in ARMS$arm[available]) {
      s <- point_scores(record, arm, estimand, linkinv)
      deviance <- if (estimand == "mean") {
        excess_deviance(
          record[, "truth"],
          record[, "truth"] + s$error,
          cell$fit_family
        )
      } else {
        NA_real_
      }
      rows[[length(rows) + 1]] <- tibble::tibble(
        task = result$task,
        cell = cell$cell,
        rep = result$rep,
        estimand = estimand,
        arm = arm,
        coverage = mean(s$covered),
        width = mean(s$width),
        score = mean(s$score),
        sq_error = mean(s$error^2),
        deviance = deviance,
        covered = list(s$covered),
        error = list(s$error)
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Loading ---------------------------------------------------------------------------------

# Reads task files; separates successes from failures.
load_results <- function(files) {
  results <- lapply(files, readRDS)
  failed <- vapply(results, \(r) !identical(r$status, "success"), TRUE)
  list(
    ok = results[!failed],
    failures = dplyr::bind_rows(lapply(results[failed], \(r) r$failure))
  )
}


# Scores of all tasks; each task's stored cell specification must match the design
# table (results made with a different design would otherwise be mislabelled).
task_scores <- function(results, cells) {
  keys <- c(
    "family",
    "fit_family",
    "error",
    "G",
    "signal",
    "truth",
    "D",
    "model",
    "basis",
    "app",
    "residual",
    "flip"
  )
  # results and designs from before the basis column have the default bases; before the
  # plasmode extension, residuals came from the truth's own fit and signs were per curve
  # (the synthetic study has no app: NA on both sides)
  spec <- \(x) {
    x <- as.list(x)
    defaults <- list(
      basis = "default",
      app = NA_character_,
      residual = x$truth,
      flip = "curve"
    )
    unlist(modifyList(defaults, x)[keys])
  }
  dplyr::bind_rows(lapply(results, \(r) {
    cell <- as.list(cells[cells$cell == r$cell, ])
    if (!identical(spec(cell), spec(r$cell_spec))) {
      stop("task ", r$task, ": stored cell spec differs from the design table")
    }
    score_task_records(r, cell)
  }))
}

# Whole study, one cell at a time (point-level vectors never all in memory): cell
# summaries, task-level scalar scores (for paired contrasts), per-fit diagnostics, per-term
# EDFs and log sp, and failures. `dirs` may hold several result directories; files are
# assigned to cells by the cell stored in each file. Differing software versions are
# reported, not refused.
summarise_study <- function(dirs, cells = design_table()) {
  files <- unlist(lapply(
    dirs,
    list.files,
    pattern = "^task-.*\\.rds$",
    full.names = TRUE
  ))
  cores <- getOption("restart.cores", 1L)
  index <- unlist(parallel::mclapply(files, mc.cores = cores, \(f) {
    r <- readRDS(f)
    if (identical(r$status, "success")) r$cell else NA_integer_
  }))
  failed <- files[is.na(index)]
  per_cell <- parallel::mclapply(
    split(files[!is.na(index)], index[!is.na(index)]),
    mc.cores = cores,
    \(f) {
      res <- load_results(f)
      scores <- task_scores(res$ok, cells)
      keyed <- \(part) {
        dplyr::bind_rows(lapply(res$ok, \(r) {
          dplyr::mutate(r[[part]], task = r$task, cell = r$cell, rep = r$rep)
        }))
      }
      list(
        software = dplyr::bind_rows(lapply(res$ok, \(r) {
          tibble::as_tibble(
            r$software %||% r$stamp[c("refund_sha", "mgcv_version")]
          )
        })),
        summary = if (nrow(scores)) summarise_cells(scores),
        tasks = if (nrow(scores)) dplyr::select(scores, -covered, -error),
        fits = keyed("fits"),
        terms = keyed("terms")
      )
    }
  )
  parts <- c(
    cells = "summary",
    tasks = "tasks",
    fits = "fits",
    terms = "terms",
    software = "software"
  )
  out <- lapply(parts, \(part) dplyr::bind_rows(lapply(per_cell, `[[`, part)))
  out$failures <- dplyr::bind_rows(lapply(failed, \(f) readRDS(f)$failure))
  # local runs record an abbreviated refund SHA, LRZ runs the full one
  out$software$refund_sha <- substr(out$software$refund_sha, 1, 8)
  out$software <- dplyr::distinct(out$software)
  if (nrow(out$software) > 1)
    warning("results come from several software versions")
  out
}

# Cell summaries ---------------------------------------------------------------------------

# Per (cell, estimand, arm): coverage, width, interval score, MSE with its split into
# integrated squared bias and variance, and the 5% quantile of pointwise coverage.
# MC SEs over replicates. SDs are computed before any column is overwritten.
summarise_cells <- function(scores) {
  scores |>
    dplyr::group_by(cell, estimand, arm) |>
    dplyr::summarise(
      n_rep = dplyr::n(),
      coverage_se = sd(coverage) / sqrt(n_rep),
      width_se = sd(width) / sqrt(n_rep),
      score_se = sd(score) / sqrt(n_rep),
      mse_se = sd(sq_error) / sqrt(n_rep),
      coverage_mean = mean(coverage),
      width_mean = mean(width),
      score_mean = mean(score),
      mse = mean(sq_error),
      deviance = mean(deviance),
      point_stats(do.call(cbind, error), do.call(cbind, covered)),
      .groups = "drop"
    ) |>
    dplyr::mutate(variance = mse - bias2)
}

# Integrated squared bias and the 5% quantile of pointwise coverage over the grid, with
# bootstrap SEs over replicates (columns of the points x replicates matrices). All
# bootstrap draws at once: W holds the resampled replicate counts, one column per draw.
point_stats <- function(
  errors,
  covered,
  n_boot = N_BOOT / 5,
  seed = SEED_BOOT
) {
  n <- ncol(errors)
  draws <- withr::with_seed(seed, {
    replicate(n_boot, tabulate(sample.int(n, replace = TRUE), n))
  })
  W <- cbind(1, matrix(draws, nrow = n)) # stays n x (1 + n_boot) even for n = 1
  m <- errors %*% W / n
  v <- (errors^2 %*% W / n - m^2) * n / (n - 1)
  # first column = the estimate, the rest = bootstrap draws (distinct names: tibble()
  # would otherwise see its own new columns)
  b2 <- colMeans(m^2 - v / n)
  q05 <- apply((covered %*% W) / n, 2, quantile, 0.05, na.rm = TRUE)
  tibble::tibble(
    bias2 = b2[1],
    bias2_se = sd(b2[-1]),
    pointwise_q05 = q05[1],
    pointwise_q05_se = sd(q05[-1])
  )
}

# Integrated squared bias from a (points x replicates) error matrix: the squared mean
# error overstates bias^2 by var/R, which is subtracted.
unbiased_bias2 <- function(errors) {
  m <- rowMeans(errors)
  v <- apply(errors, 1, var)
  mean(m^2 - v / ncol(errors))
}

# Paired contrast of two arms within cells: per-replicate differences of a task-level
# metric, and the ratio of aggregates with a replicate-bootstrap interval.
paired_contrast <- function(scores, arm_a, arm_b, metric, n_boot = N_BOOT) {
  wide <- scores |>
    dplyr::filter(arm %in% c(arm_a, arm_b)) |>
    dplyr::select(cell, estimand, rep, arm, value = dplyr::all_of(metric)) |>
    tidyr::pivot_wider(names_from = arm, values_from = value) |>
    dplyr::filter(!is.na(.data[[arm_a]]), !is.na(.data[[arm_b]])) # complete pairs only
  wide |>
    dplyr::group_by(cell, estimand) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      diff_se = sd(.data[[arm_a]] - .data[[arm_b]]) / sqrt(dplyr::n()),
      diff = mean(.data[[arm_a]] - .data[[arm_b]]),
      ratio = sum(.data[[arm_a]]) / sum(.data[[arm_b]]),
      ratio_ci = list(boot_ratio_ci(.data[[arm_a]], .data[[arm_b]], n_boot)),
      .groups = "drop"
    )
}

boot_ratio_ci <- function(a, b, n_boot, seed = SEED_BOOT) {
  ratios <- withr::with_seed(
    seed,
    replicate(n_boot, {
      i <- sample.int(length(a), replace = TRUE)
      sum(a[i]) / sum(b[i])
    })
  )
  unname(quantile(ratios, c(0.025, 0.975)))
}

# The recommendation rule --------------------------------------------------------

# Pool: core cells at G = 100 with ou/fpc errors (all signal levels), plus the
# Gaussian warp cell at G = 100 with the smooth truth.
recommendation_pool <- function(cells, family, signal = NULL) {
  core <- cells$block == "core" & cells$error %in% c("ou", "fpc")
  warp <- cells$block == "warp" &
    cells$family == "gaussian" &
    cells$truth == "smooth"
  in_pool <- (core | warp) & cells$G == 100L & cells$family == family
  if (!is.null(signal)) in_pool <- in_pool & cells$signal == signal
  cells$cell[in_pool]
}

# Root-mean-square undercoverage over cells: cells at or above nominal count as zero (width,
# and so overcoverage, is priced by the interval score); each undercovering cell's squared
# error is reduced by its MC variance so that MC noise alone does not count.
calibration_error <- function(coverage, coverage_se) {
  under <- coverage < LEVEL
  sqrt(max(0, mean(under * ((coverage - LEVEL)^2 - coverage_se^2))))
}

recommend <- function(
  summaries,
  cells,
  family,
  estimand,
  proposal = PROPOSAL,
  signal = NULL
) {
  pool <- recommendation_pool(cells, family, signal)
  arm_rows <- \(a) {
    dplyr::filter(summaries, cell %in% pool, estimand == !!estimand, arm == a)
  }
  both <- dplyr::inner_join(
    arm_rows(proposal),
    arm_rows(FALLBACK),
    by = "cell",
    suffix = c("_p", "_f")
  )
  if (!setequal(both$cell, pool)) {
    return(tibble::tibble(
      family = family,
      estimand = estimand,
      decision = "incomplete pool"
    ))
  }
  prop <- dplyr::rename_with(
    both,
    \(n) sub("_p$", "", n),
    dplyr::ends_with("_p")
  )
  fall <- dplyr::rename_with(
    both,
    \(n) sub("_f$", "", n),
    dplyr::ends_with("_f")
  )
  ce <- c(
    proposal = calibration_error(prop$coverage_mean, prop$coverage_se),
    fallback = calibration_error(fall$coverage_mean, fall$coverage_se)
  )
  ratio <- exp(mean(log(prop$score_mean / fall$score_mean)))
  adequate <- ce <= ADEQUATE
  decision <- if (sum(adequate) == 1) {
    names(ce)[adequate]
  } else if (!any(adequate)) {
    paste("neither adequate; less bad:", names(which.min(ce)))
  } else if (ratio <= PRACTICAL_RATIO) {
    "proposal"
  } else if (ratio >= 1 / PRACTICAL_RATIO) {
    "fallback"
  } else {
    "fallback (practically equivalent; simpler)"
  }
  tibble::tibble(
    family = family,
    estimand = estimand,
    n_cells = length(pool),
    ce_proposal = ce[["proposal"]],
    ce_fallback = ce[["fallback"]],
    score_ratio = ratio,
    decision = decision
  )
}
