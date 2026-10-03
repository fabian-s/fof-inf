# CL2-arms summaries: GLS (raw / fpc) + CL2 and pcre + CL2 against REML + CL2, NCV + CL2
# and the plain comparator arms, on common replicates. Run on LRZ from the restart
# directory (code copy; the cardiology task files never leave LRZ):
#   Rscript comparators/analyse-cl2-arms.R <cl2_dir> <comp_dir> <ref_syn> <ref_plas> <out_dir>
# Each replicate's record = the CL2-arms task record + the comparator task record's
# columns it lacks + the study record's NCV columns (truth and est_REML must agree).
# Two replicate windows: 1-200 (all arms but pcre + CL2) and 1-100 (all arms).
# Writes cl2-arms-{cells,contrasts,fits,crn,failures}.csv and cl2-arms-SUMMARY.md;
# no number in the summary is typed by hand.

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")

# Scoring helpers of analyse-comparators.R (function definitions only).
local({
  want <- c(
    "arm_matrices",
    "interval_points",
    "score_arm",
    "md_table",
    "cell_label"
  )
  for (e in parse("comparators/analyse-comparators.R")) {
    if (
      is.call(e) &&
        identical(e[[1]], as.name("<-")) &&
        is.name(e[[2]]) &&
        as.character(e[[2]]) %in% want
    ) {
      eval(e, globalenv())
    }
  }
})

ARMS <- tibble::tribble(
  ~arm,
  ~est,
  ~se,
  ~lo,
  ~hi,
  ~label,
  "reml_cl2",
  "est_REML",
  "se_cl2_REML",
  NA,
  NA,
  "REML + CL2",
  "ncv_cl2",
  "est_NCV",
  "se_cl2_NCV",
  NA,
  NA,
  "NCV + CL2",
  "gls_fpc",
  "est_GLSfpc",
  "se_GLSfpc",
  NA,
  NA,
  "GLS-FPCA",
  "gls_fpc_cl2",
  "est_GLSfpc",
  "se_cl2_GLSfpc",
  NA,
  NA,
  "GLS-FPCA + CL2",
  "gls_raw",
  "est_GLSraw",
  "se_GLSraw",
  NA,
  NA,
  "GLS-raw",
  "gls_raw_cl2",
  "est_GLSraw",
  "se_cl2_GLSraw",
  NA,
  NA,
  "GLS-raw + CL2",
  "pcre",
  "est_PCRE",
  "se_PCRE",
  NA,
  NA,
  "pcre",
  "pcre_cl2",
  "est_PCRE",
  "se_cl2_PCRE",
  NA,
  NA,
  "pcre + CL2"
)
ESTIMANDS <- c("beta", "mean", "alpha", "gamma")
WINDOWS <- list(`1-200` = 1:200, `1-100` = 1:100)
DETECT_SHARE <- 0.25 # detection: grid points with |truth| > this share of max |truth|

# Merging ------------------------------------------------------------------------------

merge_record <- function(r, comp_path, ref_path) {
  add_from <- function(records, path, columns = NULL) {
    if (!file.exists(path)) return(records)
    other <- readRDS(path)
    if (!identical(other$status %||% "success", "success")) return(records)
    for (e in names(records)) {
      o <- other$records[[e]]
      if (is.null(o)) next
      if (max(abs(o[, "truth"] - records[[e]][, "truth"])) > 1e-10) {
        stop("truth differs from ", path)
      }
      if (max(abs(o[, "est_REML"] - records[[e]][, "est_REML"])) > 1e-6) {
        stop("est_REML differs from ", path)
      }
      new <- setdiff(columns %||% colnames(o), colnames(records[[e]]))
      new <- intersect(new, colnames(o))
      records[[e]] <- cbind(records[[e]], o[, new, drop = FALSE])
    }
    records
  }
  r$records <- add_from(r$records, comp_path)
  r$records <- add_from(r$records, ref_path, c("est_NCV", "se_cl2_NCV"))
  r
}

# Scoring --------------------------------------------------------------------------------

detection_rate <- function(m) {
  z <- qnorm(1 - (1 - LEVEL) / 2)
  lower <- m$est - z * m$se
  upper <- m$est + z * m$se
  truth <- m$truth
  big <- abs(truth) > DETECT_SHARE * apply(abs(truth), 2, max)[col(truth)]
  hit <- (truth > 0 & lower > 0) | (truth < 0 & upper < 0)
  sum(hit & big) / sum(big)
}

paired <- function(per_rep, keys, reference) {
  ref <- dplyr::filter(per_rep, arm == reference)
  others <- setdiff(unique(per_rep$arm), reference)
  dplyr::bind_rows(lapply(others, \(a) {
    d <- dplyr::inner_join(
      dplyr::filter(per_rep, arm == a),
      ref,
      by = "rep",
      suffix = c("", "_ref")
    )
    diff <- d$coverage - d$coverage_ref
    dplyr::bind_cols(
      keys,
      tibble::tibble(
        reference = reference,
        arm = a,
        n_pairs = nrow(d),
        coverage_diff_se = sd(diff) / sqrt(length(diff)),
        coverage_diff = mean(diff),
        width_ratio = sum(d$width) / sum(d$width_ref),
        score_ratio = sum(d$score) / sum(d$score_ref),
        mse_ratio = sum(d$sq_error) / sum(d$sq_error_ref)
      )
    )
  }))
}

score_window <- function(results, base, window_name, reps) {
  results <- Filter(\(r) r$rep %in% reps, results)
  rows <- list()
  contrasts <- list()
  for (estimand in ESTIMANDS) {
    if (!length(results) || is.null(results[[1]]$records[[estimand]])) next
    keys <- dplyr::mutate(base, window = window_name, estimand = estimand)
    per_rep <- list()
    for (i in seq_len(nrow(ARMS))) {
      if (window_name == "1-200" && ARMS$arm[i] == "pcre_cl2") next
      m <- arm_matrices(results, estimand, ARMS[i, ])
      if (length(m$reps) == 0) next
      out <- score_arm(m, keys, ARMS[i, ], length(reps))
      rows[[length(rows) + 1]] <- dplyr::mutate(
        out$row,
        detection = detection_rate(m)
      )
      per_rep[[length(per_rep) + 1]] <- out$per_rep
    }
    per_rep <- dplyr::bind_rows(per_rep)
    contrasts[[estimand]] <- dplyr::bind_rows(
      paired(per_rep, keys, "reml_cl2"),
      if ("ncv_cl2" %in% per_rep$arm) paired(per_rep, keys, "ncv_cl2")
    )
  }
  list(cells = dplyr::bind_rows(rows), contrasts = dplyr::bind_rows(contrasts))
}

parse_info <- function(info, field) {
  pat <- sprintf(".*%s=([^;]+).*", field)
  ifelse(grepl(paste0(field, "="), info), sub(pat, "\\1", info), NA)
}

summarise_cell <- function(cell, cl2_dir, comp_dir, ref_syn, ref_plas) {
  files <- list.files(
    cl2_dir,
    sprintf("^task-%s-r\\d+\\.rds$", cell$key),
    full.names = TRUE
  )
  if (!length(files)) return(NULL)
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  results <- lapply(all[ok], \(r) {
    task <- list(key = r$key, rep = r$rep)
    merge_record(
      r,
      file.path(comp_dir, sprintf("task-%s-r%03d.rds", r$key, r$rep)),
      reference_path(task, cell, ref_syn, ref_plas)
    )
  })
  # locals first: tibble() evaluates sequentially, a `cell` column would shadow the row
  label <- cell_label(cell)
  base <- tibble::tibble(
    key = cell$key,
    source = cell$source,
    label = label,
    error = cell$error,
    app = cell$app,
    G = cell$G,
    cell = cell$cell
  )
  scored <- lapply(names(WINDOWS), \(w) {
    score_window(results, base, w, WINDOWS[[w]])
  })
  failures <- dplyr::bind_rows(lapply(all, \(r) {
    if (is.null(r$failures) || !NROW(r$failures)) return(NULL)
    dplyr::mutate(r$failures, key = r$key, rep = r$rep, task_status = r$status)
  }))
  fits <- dplyr::bind_rows(lapply(all[ok], \(r) {
    dplyr::mutate(r$fits, key = r$key, rep = r$rep)
  }))
  crn <- dplyr::bind_rows(lapply(all[ok], \(r) {
    if (is.null(r$crn)) return(NULL)
    dplyr::mutate(r$crn, key = r$key, rep = r$rep)
  }))
  list(
    cells = dplyr::bind_rows(lapply(scored, `[[`, "cells")),
    contrasts = dplyr::bind_rows(lapply(scored, `[[`, "contrasts")),
    failures = failures,
    fits = fits,
    crn = crn,
    tasks = tibble::tibble(
      key = cell$key,
      n_files = length(all),
      n_success = sum(ok),
      n_merged_ncv = sum(vapply(
        results,
        \(r) "est_NCV" %in% colnames(r$records$beta),
        TRUE
      )),
      task_time = sum(vapply(all, \(r) r$time %||% NA_real_, 0), na.rm = TRUE),
      peak_mb = max(vapply(all, \(r) r$peak_mb %||% NA_real_, 0), na.rm = TRUE)
    )
  )
}

# Summary -----------------------------------------------------------------------------------

order_rows <- function(df) {
  dplyr::arrange(
    df,
    factor(source, c("synthetic", "plasmode")),
    factor(error, c("iid", "ou", "fpc")),
    factor(app, COMP_APPS),
    dplyr::desc(G)
  )
}

# One row per cell x arm: coverage (MC SE), q05, width / score vs REML + CL2, MSE vs REML
# and NCV, SE/SD, detection.
long_table <- function(cells, contrasts, est, win, src, arms) {
  vs <- \(ref, metric) {
    contrasts |>
      dplyr::filter(reference == ref) |>
      dplyr::select(key, window, estimand, arm, value = dplyr::all_of(metric))
  }
  cells |>
    dplyr::filter(
      estimand == est,
      window == win,
      source == src,
      arm %in% arms
    ) |>
    dplyr::left_join(
      dplyr::rename(vs("reml_cl2", "width_ratio"), width_rel = value),
      by = c("key", "window", "estimand", "arm")
    ) |>
    dplyr::left_join(
      dplyr::rename(vs("reml_cl2", "score_ratio"), score_rel = value),
      by = c("key", "window", "estimand", "arm")
    ) |>
    dplyr::left_join(
      dplyr::rename(vs("reml_cl2", "mse_ratio"), mse_vs_reml = value),
      by = c("key", "window", "estimand", "arm")
    ) |>
    dplyr::left_join(
      dplyr::rename(vs("ncv_cl2", "mse_ratio"), mse_vs_ncv = value),
      by = c("key", "window", "estimand", "arm")
    ) |>
    dplyr::mutate(
      width_rel = ifelse(arm == "reml_cl2", 1, width_rel),
      score_rel = ifelse(arm == "reml_cl2", 1, score_rel),
      mse_vs_reml = ifelse(arm == "reml_cl2", 1, mse_vs_reml),
      mse_vs_ncv = ifelse(arm == "ncv_cl2", 1, mse_vs_ncv)
    ) |>
    order_rows() |>
    dplyr::mutate(arm = factor(ARMS$label[match(arm, ARMS$arm)], ARMS$label)) |>
    dplyr::arrange(factor(label, unique(label)), arm)
}

format_long <- function(df) {
  df |>
    dplyr::transmute(
      cell = label,
      arm = as.character(arm),
      n = n_rep,
      `cov %` = sprintf("%.1f (%.1f)", 100 * coverage, 100 * coverage_se),
      `q05 %` = sprintf("%.1f", 100 * pointwise_q05),
      `width/REML` = sprintf("%.2f", width_rel),
      `score/REML` = sprintf("%.2f", score_rel),
      `MSE/REML` = sprintf("%.2f", mse_vs_reml),
      `MSE/NCV` = sprintf("%.2f", mse_vs_ncv),
      `SE/SD` = sprintf("%.2f", se_sd_ratio),
      `detect %` = sprintf("%.1f", 100 * detection)
    )
}

# Means over cells per source and arm (synthetic: 6 cells; plasmode: 12 cells).
mean_table <- function(long) {
  long |>
    dplyr::group_by(arm) |>
    dplyr::summarise(
      cells = dplyr::n(),
      `cov %` = sprintf("%.1f", 100 * mean(coverage)),
      `min cell cov %` = sprintf("%.1f", 100 * min(coverage)),
      `q05 %` = sprintf("%.1f", 100 * mean(pointwise_q05)),
      `width/REML` = sprintf("%.2f", mean(width_rel)),
      `score/REML` = sprintf("%.2f", mean(score_rel)),
      `MSE/REML` = sprintf("%.2f", mean(mse_vs_reml)),
      `MSE/NCV` = sprintf("%.2f", mean(mse_vs_ncv)),
      `SE/SD` = sprintf("%.2f", mean(se_sd_ratio)),
      `detect %` = sprintf("%.1f", 100 * mean(detection)),
      .groups = "drop"
    ) |>
    dplyr::mutate(arm = as.character(arm))
}

fit_diagnostics <- function(fits, cells) {
  fits |>
    dplyr::filter(arm %in% c("GLSraw", "GLSfpc", "PCRE")) |>
    dplyr::left_join(
      dplyr::distinct(cells, key, label, source, error, app, G),
      by = "key"
    ) |>
    dplyr::mutate(
      hat = parse_info(info, "hat_eig"),
      hat_lo = as.numeric(sub("\\[([^,]+),.*", "\\1", hat)),
      hat_hi = as.numeric(sub(".*,([^]]+)\\]", "\\1", hat)),
      n_adj = as.numeric(parse_info(info, "n_adjusted")),
      sig2 = as.numeric(parse_info(info, "sig2")),
      cl2_time = as.numeric(parse_info(info, "cl2_time"))
    ) |>
    dplyr::group_by(source, error, app, G, label, arm) |>
    dplyr::summarise(
      reps = dplyr::n(),
      `min eig H_gg` = min(hat_lo),
      `max eig H_gg` = max(hat_hi),
      `reps with floored block` = sum(n_adj > 0),
      `median sig2` = median(sig2),
      `fit s` = mean(time),
      `CL2 s` = mean(cl2_time),
      .groups = "drop"
    ) |>
    order_rows() |>
    dplyr::select(-source, -error, -app, -G)
}

write_summary <- function(out, path) {
  crn_by <- out$crn |>
    dplyr::group_by(reference, column) |>
    dplyr::summarise(
      records = dplyr::n_distinct(paste(key, rep)),
      max_abs = formatC(max(max_abs), format = "e", digits = 1),
      max_rel = formatC(max(max_rel), format = "e", digits = 1),
      .groups = "drop"
    )
  cl2_check <- as.numeric(parse_info(
    out$fits$info[out$fits$arm == "REML"],
    "cl2_check"
  ))
  fail_tab <- if (nrow(out$failures)) {
    out$failures |>
      dplyr::mutate(message = substr(gsub("[\n|]", " ", message), 1, 90)) |>
      dplyr::count(key, stage, task_status, message, name = "n")
  } else {
    tibble::tibble(note = "no failures")
  }
  lines <- c(
    "# CL2 on GLS and pcre fits (comparator cells)",
    "",
    "Generated by `comparators/analyse-cl2-arms.R`; no number typed by hand.",
    "Nominal 95% pointwise intervals, Gaussian model M1. CL2 = exact CL2, Bayesian form (+ Vp - Ve), curves as clusters.",
    "GLS + CL2: CL2 of the whitened penalized fit (design, response and bread of the whitened data);",
    "pcre + CL2: refund's CL2 of the pcre fit. GLS-FPCA / GLS-raw / pcre without CL2: model-based Vc.",
    "Width, score: paired ratio of means to REML + CL2; MSE: paired ratio to REML (same as REML + CL2) and to NCV.",
    sprintf(
      "Detection: share of grid points with |truth| > %.0f%% of max |truth| whose interval excludes zero on the truth's side.",
      100 * DETECT_SHARE
    ),
    "",
    "## Run and checks",
    "",
    sprintf(
      "- Task files: %d (successful %d); arm failure rows: %d; records with NCV merged from the study: %d.",
      sum(out$tasks$n_files),
      sum(out$tasks$n_success),
      nrow(out$failures),
      sum(out$tasks$n_merged_ncv)
    ),
    sprintf(
      "- pcre (refit) on replicates %s; GLS arms on 1-%d.",
      paste(range(out$fits$rep[out$fits$arm == "PCRE"]), collapse = "-"),
      max(out$fits$rep)
    ),
    sprintf(
      "- Worker time: %.1f core-hours; peak memory per task %d MB.",
      sum(out$tasks$task_time) / 3600,
      as.integer(max(out$tasks$peak_mb))
    ),
    sprintf(
      "- CL2 implementation check (own exact CL2 on the REML design vs refund::pffr_vcov): max relative difference %.1e over %d fits.",
      max(cl2_check, na.rm = TRUE),
      sum(!is.na(cl2_check))
    ),
    "- CRN against the study records and the comparator task records:",
    "",
    md_table(crn_by),
    "Whitened (GLS) and pcre hat: eigenvalue range of the curve blocks H_gg, replicates with a floored residual block,",
    "GLS scale estimate, mean fit and CL2 time (s):",
    "",
    md_table(fit_diagnostics(out$fits, out$cells), digits = 3),
    "Failures:",
    "",
    md_table(fail_tab)
  )
  sections <- list(
    list(win = "1-200", arms = setdiff(ARMS$arm, "pcre_cl2")),
    list(win = "1-100", arms = ARMS$arm)
  )
  for (sec in sections) {
    lines <- c(lines, sprintf("## Replicates %s", sec$win), "")
    for (est in ESTIMANDS) {
      for (src in c("synthetic", "plasmode")) {
        long <- long_table(
          out$cells,
          out$contrasts,
          est,
          sec$win,
          src,
          sec$arms
        )
        if (!nrow(long)) next
        lines <- c(
          lines,
          sprintf("### %s, %s: mean over cells", est, src),
          "",
          md_table(mean_table(long)),
          if (est == "beta" || sec$win == "1-100") {
            c(
              sprintf("### %s, %s: per cell", est, src),
              "",
              md_table(format_long(long))
            )
          }
        )
      }
    }
  }
  writeLines(lines, path)
}

# Main -----------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: analyse-cl2-arms.R cl2_dir comp_dir ref_syn ref_plas out_dir")
}
cells <- dplyr::bind_rows(
  comparator_syn_cells() |> dplyr::mutate(app = NA_character_),
  comparator_plas_cells()
)
per_cell <- lapply(seq_len(nrow(cells)), \(i) {
  summarise_cell(cells[i, ], args[1], args[2], args[3], args[4])
})
per_cell <- Filter(Negate(is.null), per_cell)
parts <- c("cells", "contrasts", "failures", "fits", "crn", "tasks")
out <- lapply(
  setNames(nm = parts),
  \(p) dplyr::bind_rows(lapply(per_cell, `[[`, p))
)
if (!nrow(out$failures)) {
  out$failures <- tibble::tibble(
    stage = character(),
    message = character(),
    key = character(),
    rep = integer(),
    task_status = character()
  )
}
out_dir <- args[5]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (p in c("cells", "contrasts", "failures", "crn")) {
  write.csv(
    out[[p]],
    file.path(out_dir, sprintf("cl2-arms-%s.csv", p)),
    row.names = FALSE
  )
}
write.csv(
  out$fits[, c("key", "rep", "arm", "time", "converged", "n_warnings", "info")],
  file.path(out_dir, "cl2-arms-fits.csv"),
  row.names = FALSE
)
write_summary(out, file.path(out_dir, "cl2-arms-SUMMARY.md"))
cat("done:", nrow(out$cells), "rows\n")
