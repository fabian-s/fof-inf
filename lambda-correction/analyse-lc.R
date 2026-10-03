# Lambda-correction study: cell summaries, paired contrasts, smoothing-parameter
# calibration and SUMMARY.md from the task records of run-lc.R. Run from the repository
# directory (code copy):
#   Rscript lambda-correction/analyse-lc.R <results_dir> <ablation_results_dir> <out_dir>
# The fixed-lambda CL2 arm is read from the covariance ablation's records (same cells,
# replicates and data). Scoring definitions as in ablation/analyse-ablation.R.
# No number in SUMMARY.md is typed by hand.

source("scripts/setup.R")
source("ablation/ablation-lib.R")

CRITS <- c("z", "t")
LC_ARMS <- c(
  "cl2",
  "naive",
  "vc",
  "ij_lin",
  "ij_linx",
  "os_lin",
  "os_linx",
  "os_nl",
  "os_nlx",
  "ij_reg",
  "ij_regx",
  "os_reg",
  "os_regx",
  "jk_lin",
  "jk_linx",
  "jk_nl",
  "jk_nlx",
  "jk_full"
)
ARM_LABELS <- c(
  cl2 = "CL2",
  naive = "+mgcv J",
  vc = "+mgcv Vc",
  ij_lin = "IJ",
  ij_linx = "IJ×",
  os_lin = "OS",
  os_linx = "OS×",
  os_nl = "OSnl",
  os_nlx = "OSnl×",
  jk_lin = "JK",
  jk_linx = "JK×",
  jk_nl = "JKnl",
  jk_nlx = "JKnl×",
  jk_full = "JKfull",
  fix = "fixλ",
  ij_reg = "IJr",
  ij_regx = "IJr×",
  os_reg = "OSr",
  os_regx = "OSr×"
)

critical_value <- function(crit, G) {
  q <- 1 - (1 - LEVEL) / 2
  if (crit == "z") qnorm(q) else qt(q, G - 1)
}

# Identifying columns of a cell (locals first: tibble() evaluates sequentially, so a new
# `cell` column would shadow the cell row).
cell_keys <- function(cell) {
  id <- cell$cell
  fam <- cell$family
  err <- cell$error
  g <- cell$G
  tibble::tibble(cell = id, family = fam, error = err, G = g)
}

# Per-replicate scores and the summary row of one cell x estimand x arm x crit.
score_arm <- function(est, se, truth, reps, cell, estimand, arm, crit) {
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  s <- interval_scores(
    est,
    se,
    truth,
    critical_value(crit, cell$G),
    estimand,
    linkinv
  )
  per_rep <- tibble::tibble(
    rep = reps,
    coverage = colMeans(s$covered),
    width = colMeans(s$width),
    score = colMeans(s$score)
  )
  n <- length(reps)
  row <- dplyr::bind_cols(
    cell_keys(cell),
    tibble::tibble(
      estimand = estimand,
      arm = arm,
      crit = crit,
      n_rep = n,
      coverage = mean(per_rep$coverage),
      coverage_se = sd(per_rep$coverage) / sqrt(n),
      pointwise_q05 = unname(quantile(rowMeans(s$covered), 0.05)),
      width = mean(per_rep$width),
      score = mean(per_rep$score),
      se_sd_ratio = mean(sqrt(rowMeans(se^2)) / apply(est, 1, sd)),
      # robust: median SE over replicates / MC SD (RMS is driven by a few very large SEs)
      se_sd_median = mean(apply(se, 1, median) / apply(est, 1, sd))
    )
  )
  list(row = row, per_rep = dplyr::mutate(per_rep, arm = arm, crit = crit))
}

# Fixed-lambda columns of the matching ablation record (or NULL).
ablation_fix <- function(abl_dir, cell, rep) {
  path <- file.path(abl_dir, sprintf("task-c%03d-r%03d.rds", cell, rep))
  if (!file.exists(path)) return(NULL)
  r <- readRDS(path)
  if (!identical(r$status, "success")) return(NULL)
  r$records
}

summarise_cell <- function(files, cell, abl_dir) {
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  res <- all[ok]
  fix <- lapply(res, \(r) ablation_fix(abl_dir, r$cell, r$rep))
  rows <- list()
  contrasts <- list()
  fix_truth_diff <- 0
  for (estimand in model_estimands(cell$model)) {
    reps <- vapply(res, `[[`, 0, "rep")
    pick <- \(col)
      do.call(cbind, lapply(res, \(r) r$records[[estimand]][, col]))
    est <- pick("est_REML")
    truth <- pick("truth")
    per_rep <- list()
    for (arm in intersect(
      LC_ARMS,
      sub("^se_", "", colnames(res[[1]]$records[[estimand]]))
    )) {
      se <- pick(paste0("se_", arm))
      for (crit in CRITS) {
        out <- score_arm(est, se, truth, reps, cell, estimand, arm, crit)
        rows[[length(rows) + 1]] <- out$row
        per_rep[[length(per_rep) + 1]] <- out$per_rep
      }
    }
    has_fix <- !vapply(fix, is.null, TRUE)
    if (any(has_fix)) {
      fpick <- \(col)
        do.call(cbind, lapply(fix[has_fix], \(f) f[[estimand]][, col]))
      fix_truth_diff <- max(
        fix_truth_diff,
        max(abs(fpick("truth") - truth[, has_fix]))
      )
      for (crit in CRITS) {
        out <- score_arm(
          fpick("est_FIX"),
          fpick("se_fix_cl2e_b"),
          fpick("truth"),
          reps[has_fix],
          cell,
          estimand,
          "fix",
          crit
        )
        rows[[length(rows) + 1]] <- out$row
        per_rep[[length(per_rep) + 1]] <- out$per_rep
      }
    }
    # paired coverage contrasts of every arm against CL2 (same crit)
    pr <- dplyr::bind_rows(per_rep)
    base <- dplyr::filter(pr, arm == "cl2")
    contrasts[[estimand]] <- pr |>
      dplyr::filter(arm != "cl2") |>
      dplyr::inner_join(base, by = c("rep", "crit"), suffix = c("", "_cl2")) |>
      dplyr::group_by(arm, crit) |>
      dplyr::summarise(
        n_rep = dplyr::n(),
        diff = mean(coverage - coverage_cl2),
        diff_se = sd(coverage - coverage_cl2) / sqrt(dplyr::n()),
        .groups = "drop"
      ) |>
      dplyr::mutate(estimand = estimand, .before = 1) |>
      (\(d) dplyr::bind_cols(cell_keys(cell)[rep(1, nrow(d)), ], d))()
  }
  # beta coverage by tertile (within cell) of the ff penalty's log lambda-hat (first ff penalty)
  regime <- targeting <- NULL
  ff <- grep("L\\.X1$", names(res[[1]]$omega$hat))
  if (length(ff) == 1 && "beta" %in% names(res[[1]]$records)) {
    q <- critical_value("z", cell$G)
    arms_b <- intersect(
      LC_ARMS,
      sub("^se_", "", colnames(res[[1]]$records$beta))
    )
    per <- do.call(
      rbind,
      lapply(res, \(r) {
        m <- r$records$beta
        cov <- vapply(
          arms_b,
          \(a)
            mean(
              abs(m[, "est_REML"] - m[, "truth"]) <= q * m[, paste0("se_", a)]
            ),
          0
        )
        wid <- vapply(
          arms_b,
          \(a) mean(m[, paste0("se_", a)] / m[, "se_cl2"]),
          0
        )
        miss_cl2 <- mean(
          abs(m[, "est_REML"] - m[, "truth"]) > q * m[, "se_cl2"]
        )
        data.frame(
          rep = r$rep,
          rho_ff = r$omega$hat[[ff]],
          arm = arms_b,
          coverage = cov,
          se_ratio = wid,
          miss_cl2 = miss_cl2
        )
      })
    )
    # targeting: does an arm inflate the SE more where CL2 misses more? (within-cell Spearman)
    targeting <- per |>
      dplyr::filter(arm != "cl2") |>
      dplyr::group_by(arm) |>
      dplyr::summarise(
        spearman = suppressWarnings(cor(
          se_ratio,
          miss_cl2,
          method = "spearman"
        )),
        .groups = "drop"
      ) |>
      (\(d) dplyr::bind_cols(cell_keys(cell)[rep(1, nrow(d)), ], d))()
    reps <- unique(per[, c("rep", "rho_ff")])
    reps$tertile <- cut(
      rank(reps$rho_ff, ties.method = "first"),
      3,
      labels = c("low", "mid", "high")
    )
    per <- merge(per, reps[, c("rep", "tertile")])
    regime <- per |>
      dplyr::group_by(tertile, arm) |>
      dplyr::summarise(
        n = dplyr::n(),
        rho_ff = mean(rho_ff),
        coverage_se = sd(coverage) / sqrt(dplyr::n()),
        coverage = mean(coverage),
        se_ratio = mean(se_ratio),
        .groups = "drop"
      ) |>
      (\(d) dplyr::bind_cols(cell_keys(cell)[rep(1, nrow(d)), ], d))()
  }
  # smoothing-parameter calibration: MC SD of omega-hat vs estimated SDs
  om_hat <- do.call(rbind, lapply(res, \(r) r$omega$hat))
  colnames(om_hat)[colnames(om_hat) %in% c("", NA)] <- "log_scale" # unnamed in early records
  calib <- dplyr::bind_rows(lapply(
    c("naive", "ij", "os", "jk", "ij_reg", "os_reg"),
    \(src) {
      s <- do.call(rbind, lapply(res, \(r) r$omega[[src]]))
      if (is.null(s)) return(NULL)
      dplyr::bind_cols(
        cell_keys(cell)[rep(1, ncol(om_hat)), ],
        tibble::tibble(
          par = colnames(om_hat) %||% paste0("omega", seq_len(ncol(om_hat))),
          source = src,
          mc_sd = apply(om_hat, 2, sd),
          mc_mad = apply(om_hat, 2, mad),
          mc_iqr = apply(om_hat, 2, IQR) / 1.349,
          rms_est_sd = sqrt(colMeans(s^2)),
          median_est_sd = apply(s, 2, median)
        )
      )
    }
  ))
  list(
    cells = dplyr::bind_rows(rows),
    contrasts = dplyr::bind_rows(contrasts),
    calib = calib,
    regime = regime,
    targeting = targeting,
    crn = dplyr::bind_rows(lapply(
      res,
      \(r) if (!is.null(r$crn)) dplyr::mutate(r$crn, cell = r$cell, rep = r$rep)
    )),
    failures = dplyr::bind_rows(lapply(
      all[!ok],
      \(r)
        tibble::tibble(
          cell = r$cell,
          rep = r$rep,
          stage = r$stage,
          message = r$message
        )
    )),
    times = dplyr::bind_rows(lapply(
      res,
      \(r)
        tibble::tibble(
          cell = r$cell,
          rep = r$rep,
          total = r$time,
          fit = r$fit_time,
          !!!as.list(unlist(r$time_parts))
        )
    )),
    jk_nonconv = sum(
      vapply(res, \(r) r$n_jk_nonconv %||% NA_real_, 0),
      na.rm = TRUE
    ),
    n_files = length(files),
    fix_truth_diff = fix_truth_diff
  )
}

# Markdown helpers ----------------------------------------------------------------------

md_table <- function(df, digits = 1) {
  fmt <- \(x)
    if (is.numeric(x)) formatC(x, format = "f", digits = digits) else
      as.character(x)
  body <- vapply(
    seq_len(nrow(df)),
    \(i) paste0("| ", paste(vapply(df[i, ], fmt, ""), collapse = " | "), " |"),
    ""
  )
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body,
    ""
  )
}

wide <- function(cells, value, crit_ = "z", estimand_, arms, scale = 1) {
  cells |>
    dplyr::filter(crit == crit_, estimand == estimand_, arm %in% arms) |>
    dplyr::mutate(
      value = .data[[value]] * scale,
      arm = factor(ARM_LABELS[arm], ARM_LABELS[arms])
    ) |>
    dplyr::select(family, error, G, arm, value) |>
    tidyr::pivot_wider(names_from = arm, values_from = value) |>
    dplyr::arrange(
      factor(family, c("gaussian", "poisson", "binomial")),
      factor(error, c("iid", "ou", "fpc")),
      G
    )
}

with_mean <- function(df) {
  m <- dplyr::summarise(df, dplyr::across(where(is.numeric), mean))
  m$G <- NA
  dplyr::bind_rows(
    dplyr::mutate(df, G = as.character(G)),
    dplyr::mutate(m, family = "mean", error = "", G = "")
  )
}

# Main ----------------------------------------------------------------------------------

# (analyse-lc-plasmode.R sources the functions above only)
if (!exists("LC_ANALYSE_LIB_ONLY")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3)
    stop("usage: analyse-lc.R results_dir ablation_results_dir out_dir")
  results_dir <- args[1]
  abl_dir <- args[2]
  out_dir <- args[3]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cells <- ablation_cells()
  files_all <- list.files(
    results_dir,
    "^task-c\\d+-r\\d+\\.rds$",
    full.names = TRUE
  )
  out <- lapply(seq_len(nrow(cells)), \(i) {
    files <- grep(
      sprintf("task-c%03d-", cells$cell[i]),
      files_all,
      value = TRUE
    )
    if (!length(files)) return(NULL)
    summarise_cell(files, cells[i, ], abl_dir)
  })
  out <- out[!vapply(out, is.null, TRUE)]
  bind <- \(name) dplyr::bind_rows(lapply(out, `[[`, name))
  res <- list(
    cells = bind("cells"),
    contrasts = bind("contrasts"),
    calib = bind("calib"),
    regime = bind("regime"),
    targeting = bind("targeting"),
    crn = bind("crn"),
    failures = bind("failures"),
    times = bind("times")
  )
  for (nm in names(res))
    write.csv(
      res[[nm]],
      file.path(out_dir, paste0("lc-", nm, ".csv")),
      row.names = FALSE
    )

  main_arms <- c(
    "cl2",
    "naive",
    "vc",
    "ij_lin",
    "os_lin",
    "jk_lin",
    "os_nl",
    "jk_nl",
    "jk_full",
    "fix"
  )
  cross_arms <- c(
    "cl2",
    "ij_linx",
    "os_linx",
    "jk_linx",
    "os_nlx",
    "jk_nlx",
    "fix"
  )
  md <- c(
    "# Cluster-robust smoothing-parameter correction for CL2 (core cells, mid signal)",
    "",
    "Generated by `lambda-correction/analyse-lc.R`. Coverage of nominal 95% pointwise intervals (%),",
    "grid-averaged. All arms Bayesian CL2 (exact adjustment) plus a smoothing-parameter term:",
    "+mgcv J = J [H^-1]_rho J' (mgcv's first-order term); +mgcv Vc = Vc - Vp; IJ / OS / JK =",
    "Var(rho-hat) from the infinitesimal jackknife / one-step leave-curve-out / exact curve",
    "jackknife, mapped linearly through J; nl = curve-wise rho-induced change of theta-hat",
    "computed by refitting theta at the leave-out rho; × = with the cross term (total curve",
    "influence); JKfull = refit jackknife of theta-hat + (Vp - Ve); fixλ = ablation's CL2 at",
    "fixed (cell geometric-mean) lambda.",
    "",
    "## Run and checks",
    "",
    sprintf(
      "- Task files: %d; failures: %d; fits with jackknife refits not fully converged: %d refits in total.",
      sum(vapply(out, `[[`, 0, "n_files")),
      nrow(res$failures),
      sum(vapply(out, `[[`, 0, "jk_nonconv"))
    ),
    sprintf(
      "- Replicates per cell: %d to %d.",
      min(res$cells$n_rep),
      max(res$cells$n_rep)
    ),
    sprintf(
      "- CRN vs study records: max |diff| truth %.1e, est %.1e, CL2 SE %.1e (relative %.1e).",
      max(res$crn$max_abs[res$crn$column == "truth"]),
      max(res$crn$max_abs[res$crn$column == "est_REML"]),
      max(res$crn$max_abs[res$crn$column == "se_cl2"]),
      max(res$crn$max_rel[res$crn$column == "se_cl2"])
    ),
    sprintf(
      "- Fixed-lambda records joined from the ablation: max |truth difference| %.1e.",
      max(vapply(out, `[[`, 0, "fix_truth_diff"))
    ),
    sprintf(
      "- Time per task: median %.0f s (fit %.0f s, jackknife %.0f s), total %.1f core-hours.",
      median(res$times$total),
      median(res$times$fit),
      median(res$times$jackknife, na.rm = TRUE),
      sum(res$times$total) / 3600
    ),
    ""
  )
  for (est in c("beta", "alpha", "gamma", "mean")) {
    md <- c(
      md,
      paste0("## ", est),
      "",
      "Coverage (%), z:",
      "",
      md_table(with_mean(wide(
        res$cells,
        "coverage",
        "z",
        est,
        main_arms,
        100
      ))),
      "Coverage (%), t_{G-1}:",
      "",
      md_table(with_mean(wide(
        res$cells,
        "coverage",
        "t",
        est,
        main_arms,
        100
      ))),
      "With cross term, coverage (%), z:",
      "",
      md_table(with_mean(wide(
        res$cells,
        "coverage",
        "z",
        est,
        cross_arms,
        100
      ))),
      "5% quantile of pointwise coverage (%), z:",
      "",
      md_table(with_mean(wide(
        res$cells,
        "pointwise_q05",
        "z",
        est,
        main_arms,
        100
      ))),
      "SE/SD:",
      "",
      md_table(
        with_mean(wide(
          res$cells,
          "se_sd_ratio",
          "z",
          est,
          c(main_arms, cross_arms[-c(1, 7)])
        )),
        2
      ),
      "Width relative to CL2, z:",
      "",
      md_table(
        with_mean(
          wide(
            res$cells,
            "width",
            "z",
            est,
            c(main_arms, cross_arms[-c(1, 7)])
          ) |>
            dplyr::mutate(dplyr::across(-c(family, error, G), \(x) x / CL2))
        ),
        3
      ),
      "Interval score relative to CL2, z:",
      "",
      md_table(
        with_mean(
          wide(
            res$cells,
            "score",
            "z",
            est,
            c(main_arms, cross_arms[-c(1, 7)])
          ) |>
            dplyr::mutate(dplyr::across(-c(family, error, G), \(x) x / CL2))
        ),
        3
      )
    )
  }
  cal <- res$calib |>
    dplyr::mutate(ratio = rms_est_sd / mc_sd) |>
    dplyr::select(family, error, G, par, source, ratio) |>
    tidyr::pivot_wider(names_from = source, values_from = ratio)
  md <- c(
    md,
    "## Calibration of the smoothing-parameter variance",
    "",
    "RMS of the estimated SD of each omega-hat component (log sp; Gaussian: also log scale)",
    "divided by its Monte Carlo SD over replicates:",
    "",
    md_table(cal, 2)
  )
  writeLines(md, file.path(out_dir, "SUMMARY.md"))
  cat("wrote", file.path(out_dir, "SUMMARY.md"), "\n")
} # main
