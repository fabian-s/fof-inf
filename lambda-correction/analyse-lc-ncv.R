# NCV lambda correction: cell summaries, paired contrasts, calibration, targeting
# and SUMMARY-ncv.md from the task records of run-lc-ncv.R / run-lc-ncvg.R. Run from the repository root:
#   Rscript lambda-correction/analyse-lc-ncv.R <results_dir> <out_dir> [families]
# families: comma-separated core families (default gaussian).

LC_ANALYSE_LIB_ONLY <- TRUE
source("lambda-correction/analyse-lc.R")

NCV_ARMS <- c(
  cl2 = "CL2",
  naive_mgcv = "+mgcv J",
  naive = "+J H^-1 J'",
  ij_lin = "IJ",
  os_lin = "OS",
  ij_reg = "IJr",
  os_reg = "OSr",
  ij_linx = "IJ×",
  os_linx = "OS×",
  ij_regx = "IJr×",
  fix = "fixλ (oracle)"
)

ncv_cell <- function(files, cell) {
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  res <- all[ok]
  rows <- per_rep <- list()
  for (estimand in model_estimands(cell$model)) {
    pick <- \(rs, col)
      do.call(cbind, lapply(rs, \(r) r$records[[estimand]][, col]))
    reps <- vapply(res, `[[`, 0, "rep")
    have <- sub("^se_", "", colnames(res[[1]]$records[[estimand]]))
    for (arm in intersect(names(NCV_ARMS), have)) {
      rs <- if (arm == "fix") Filter(\(r) isTRUE(r$fixed_ok), res) else res
      est <- pick(rs, if (arm == "fix") "est_FIX" else "est_NCV")
      se <- pick(rs, paste0("se_", arm))
      for (crit in CRITS) {
        out <- score_arm(
          est,
          se,
          pick(rs, "truth"),
          vapply(rs, `[[`, 0, "rep"),
          cell,
          estimand,
          arm,
          crit
        )
        rows[[length(rows) + 1]] <- out$row
        per_rep[[length(per_rep) + 1]] <- dplyr::mutate(
          out$per_rep,
          estimand = estimand
        )
      }
    }
  }
  pr <- dplyr::bind_rows(per_rep)
  base <- dplyr::filter(pr, arm == "cl2")
  contrasts <- pr |>
    dplyr::filter(arm != "cl2") |>
    dplyr::inner_join(
      base,
      by = c("rep", "crit", "estimand"),
      suffix = c("", "_cl2")
    ) |>
    dplyr::group_by(estimand, arm, crit) |>
    dplyr::summarise(
      n_rep = dplyr::n(),
      diff = mean(coverage - coverage_cl2),
      diff_se = sd(coverage - coverage_cl2) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    (\(d) dplyr::bind_cols(cell_keys(cell)[rep(1, nrow(d)), ], d))()
  # calibration of log lambda-hat SDs and targeting (beta)
  om_hat <- do.call(rbind, lapply(res, \(r) r$omega$hat))
  calib <- dplyr::bind_rows(lapply(
    setdiff(names(res[[1]]$omega), c("hat", "boundary")),
    \(src) {
      s <- do.call(rbind, lapply(res, \(r) r$omega[[src]]))
      if (is.null(s)) return(NULL)
      dplyr::bind_cols(
        cell_keys(cell)[rep(1, ncol(om_hat)), ],
        tibble::tibble(
          par = colnames(om_hat),
          source = src,
          mc_sd = apply(om_hat, 2, sd),
          mc_iqr = apply(om_hat, 2, IQR) / 1.349,
          median_est_sd = apply(s, 2, median)
        )
      )
    }
  ))
  q <- critical_value("z", cell$G)
  have_b <- sub("^se_", "", colnames(res[[1]]$records$beta))
  tgt <- dplyr::bind_rows(lapply(
    setdiff(intersect(names(NCV_ARMS), have_b), c("cl2", "fix")),
    \(a) {
      x <- t(vapply(
        res,
        \(r) {
          m <- r$records$beta
          c(
            ratio = mean(m[, paste0("se_", a)] / m[, "se_cl2"]),
            miss = mean(abs(m[, "est_NCV"] - m[, "truth"]) > q * m[, "se_cl2"])
          )
        },
        numeric(2)
      ))
      tibble::tibble(
        arm = a,
        spearman = suppressWarnings(cor(
          x[, "ratio"],
          x[, "miss"],
          method = "spearman"
        )),
        infl_median = median(x[, "ratio"]),
        infl_q90 = unname(quantile(x[, "ratio"], 0.9)),
        infl_max = max(x[, "ratio"])
      )
    }
  )) |>
    (\(d) dplyr::bind_cols(cell_keys(cell)[rep(1, nrow(d)), ], d))()
  list(
    cells = dplyr::bind_rows(rows),
    contrasts = contrasts,
    calib = calib,
    targeting = tgt,
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
          converged = isTRUE(r$converged),
          total = r$time,
          fit = r$fit_time,
          !!!as.list(unlist(r$time_parts))
        )
    )),
    boundary = if (!is.null(res[[1]]$omega$boundary)) {
      b <- do.call(rbind, lapply(res, \(r) r$omega$boundary))
      dplyr::bind_cols(
        cell_keys(cell)[rep(1, ncol(om_hat)), ],
        tibble::tibble(
          par = colnames(om_hat),
          n_rep = nrow(b),
          rate = colMeans(b)
        )
      )
    },
    n_fixed_failed = sum(!vapply(res, \(r) isTRUE(r$fixed_ok), TRUE)),
    max_grad = max(vapply(res, \(r) max(abs(r$grad)), 0)),
    max_check = max(vapply(res, \(r) max(c(0, r$check)), 0))
  )
}

if (!exists("LC_NCV_LIB_ONLY")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2)
    stop("usage: analyse-lc-ncv.R results_dir out_dir [families]")
  fams <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else "gaussian"
  dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
  cells <- ablation_cells()
  cells <- cells[cells$family %in% fams, ]
  files_all <- list.files(
    args[1],
    "^task-c\\d+-r\\d+\\.rds$",
    full.names = TRUE
  )
  out <- lapply(seq_len(nrow(cells)), \(i) {
    f <- grep(sprintf("task-c%03d-", cells$cell[i]), files_all, value = TRUE)
    if (length(f)) ncv_cell(f, cells[i, ])
  })
  out <- out[!vapply(out, is.null, TRUE)]
  bind <- \(nm) dplyr::bind_rows(lapply(out, `[[`, nm))
  res <- list(
    cells = bind("cells"),
    contrasts = bind("contrasts"),
    calib = bind("calib"),
    targeting = bind("targeting"),
    boundary = bind("boundary"),
    crn = bind("crn"),
    failures = bind("failures"),
    times = bind("times")
  )
  for (nm in names(res))
    write.csv(
      res[[nm]],
      file.path(args[2], paste0("ncv-", nm, ".csv")),
      row.names = FALSE
    )

  tab <- \(value, crit_ = "z", ratio = FALSE, scale = 100, digits = 1) {
    d <- dplyr::filter(res$cells, crit == crit_)
    if (ratio)
      d <- d |>
        dplyr::group_by(cell, estimand) |>
        dplyr::mutate(
          !!value := .data[[value]] / .data[[value]][arm == "cl2"]
        ) |>
        dplyr::ungroup()
    d |>
      dplyr::mutate(
        v = .data[[value]] * scale,
        arm = factor(NCV_ARMS[arm], NCV_ARMS)
      ) |>
      dplyr::mutate(G = as.character(G)) |>
      dplyr::select(estimand, family, error, G, arm, v) |>
      tidyr::pivot_wider(names_from = arm, values_from = v) |>
      dplyr::arrange(
        factor(estimand, c("beta", "alpha", "gamma", "mean")),
        family,
        factor(error, c("iid", "ou", "fpc")),
        as.integer(G)
      ) |>
      dplyr::mutate(dplyr::across(
        where(is.numeric),
        \(x) formatC(x, format = "f", digits = digits)
      )) |>
      as.data.frame() |>
      md_table(digits)
  }
  md <- c(
    sprintf(
      "# Lambda correction for curve-blocked NCV (core %s cells, mid signal)",
      paste(fams, collapse = "/")
    ),
    "",
    "Generated by `lambda-correction/analyse-lc-ncv.R`. Arms: CL2 as in the study (on Vp);",
    "+mgcv J = J [H^-1] J' with mgcv's NCV Hessian; +J H^-1 J' with the exact NCV Hessian;",
    "IJ / OS = Var(rho-hat) from the NCV curve scores / one-step leave-curve-out NCV; r = H",
    "regularised (prior 1/10); × = with the cross term; fixλ = CL2 of the fit at the cell's",
    "geometric-mean NCV lambda (oracle).",
    "",
    sprintf(
      "- Task files: %d; failures: %d; fixed-lambda refits failed: %d; replicates per cell %d-%d.",
      length(files_all),
      nrow(res$failures),
      sum(vapply(out, `[[`, 0, "n_fixed_failed")),
      min(res$cells$n_rep[res$cells$arm == "cl2"]),
      max(res$cells$n_rep[res$cells$arm == "cl2"])
    ),
    sprintf(
      "- CRN vs study NCV records: max relative difference est %.1e, CL2 SE %.1e.",
      max(res$crn$max_rel[res$crn$column == "est_NCV"]),
      max(res$crn$max_rel[res$crn$column == "se_cl2"])
    ),
    sprintf(
      "- Max |NCV gradient| at the fits: %.1e. Time per task: median %.0f s, total %.1f core-hours.",
      max(vapply(out, `[[`, 0, "max_grad")),
      median(res$times$total),
      sum(res$times$total) / 3600
    ),
    sprintf(
      "- Max consistency-check discrepancy (Vp, J, NCV value and gradient vs mgcv; GLM runs): %.1e.",
      max(vapply(out, `[[`, 0, "max_check"))
    ),
    "",
    "## Coverage (%), z",
    "",
    tab("coverage"),
    "## Coverage (%), t_{G-1}",
    "",
    tab("coverage", "t"),
    "## 5% quantile of pointwise coverage (%), z",
    "",
    tab("pointwise_q05"),
    "## SE/SD",
    "",
    tab("se_sd_ratio", scale = 1, digits = 2),
    "## Width relative to CL2",
    "",
    tab("width", ratio = TRUE, scale = 1, digits = 3),
    "## Interval score relative to CL2",
    "",
    tab("score", ratio = TRUE, scale = 1, digits = 3)
  )
  writeLines(md, file.path(args[2], "SUMMARY-ncv.md"))
  cat("wrote", file.path(args[2], "SUMMARY-ncv.md"), "\n")
} # main
