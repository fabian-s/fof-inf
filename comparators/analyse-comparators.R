# Comparator study: cell summaries, paired contrasts against REML + CL2, timing,
# failures, CRN check and SUMMARY.md from the task records of run-comparators.R.
# Run from the repository root (code copy; the cardiology task files never leave LRZ):
#   Rscript comparators/analyse-comparators.R <results_dir> <out_dir>
# No number in SUMMARY.md is typed by hand.

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")

# Interval arms: centre column, and either an SE column (Wald, z_0.975) or limit columns.
ARMS <- tibble::tribble(
  ~arm,
  ~est,
  ~se,
  ~lo,
  ~hi,
  ~label,
  "reml_model",
  "est_REML",
  "se_model_REML",
  NA,
  NA,
  "REML model",
  "reml_cl2",
  "est_REML",
  "se_cl2_REML",
  NA,
  NA,
  "REML CL2",
  "pcre",
  "est_PCRE",
  "se_PCRE",
  NA,
  NA,
  "pcre",
  "gls_raw",
  "est_GLSraw",
  "se_GLSraw",
  NA,
  NA,
  "GLS raw",
  "gls_fpc",
  "est_GLSfpc",
  "se_GLSfpc",
  NA,
  NA,
  "GLS fpc",
  "boot_pct",
  "est_REML",
  NA,
  "lo_BOOT",
  "hi_BOOT",
  "boot pct",
  "boot_wald",
  "est_REML",
  "se_BOOT",
  NA,
  NA,
  "boot Wald"
)
REFERENCE <- "reml_cl2"
ESTIMANDS <- c("beta", "mean", "alpha", "gamma")
FIT_ARMS <- c(
  REML = "reml",
  PCRE = "pcre",
  GLSraw = "gls_raw",
  GLSfpc = "gls_fpc",
  BOOT = "boot"
)

# Scoring ---------------------------------------------------------------------------------

arm_matrices <- function(results, estimand, spec) {
  need <- na.omit(c(spec$est, spec$se, spec$lo, spec$hi))
  ok <- vapply(
    results,
    \(r) all(need %in% colnames(r$records[[estimand]])),
    TRUE
  )
  pick <- \(col) {
    if (is.na(col)) return(NULL)
    do.call(cbind, lapply(results[ok], \(r) r$records[[estimand]][, col]))
  }
  list(
    reps = vapply(results[ok], \(r) r$rep, 0),
    est = pick(spec$est),
    se = pick(spec$se),
    lo = pick(spec$lo),
    hi = pick(spec$hi),
    truth = pick("truth")
  )
}

interval_points <- function(m) {
  z <- qnorm(1 - (1 - LEVEL) / 2)
  lower <- if (is.null(m$se)) m$lo else m$est - z * m$se
  upper <- if (is.null(m$se)) m$hi else m$est + z * m$se
  alpha <- 1 - LEVEL
  truth <- m$truth
  list(
    covered = lower <= truth & truth <= upper,
    width = upper - lower,
    score = (upper - lower) +
      2 /
        alpha *
        ((lower - truth) * (truth < lower) + (truth - upper) * (truth > upper)),
    sq_error = (m$est - truth)^2
  )
}

score_arm <- function(m, keys, spec, n_total) {
  s <- interval_points(m)
  n <- ncol(m$est)
  per_rep <- tibble::tibble(
    rep = m$reps,
    coverage = colMeans(s$covered),
    width = colMeans(s$width),
    score = colMeans(s$score),
    sq_error = colMeans(s$sq_error)
  )
  se_sd <- if (is.null(m$se) || n < 3) NA_real_ else {
    mean(sqrt(rowMeans(m$se^2)) / apply(m$est, 1, sd))
  }
  row <- dplyr::bind_cols(
    keys,
    tibble::tibble(
      arm = spec$arm,
      n_rep = n,
      n_missing = n_total - n,
      coverage_se = sd(per_rep$coverage) / sqrt(n),
      coverage = mean(per_rep$coverage),
      pointwise_q05 = unname(quantile(rowMeans(s$covered), 0.05)),
      width = mean(per_rep$width),
      score = mean(per_rep$score),
      se_sd_ratio = se_sd,
      mse_se = sd(per_rep$sq_error) / sqrt(n),
      mse = mean(per_rep$sq_error)
    )
  )
  list(row = row, per_rep = dplyr::mutate(per_rep, arm = spec$arm))
}

# Paired contrasts against REML + CL2 on the replicates both arms have.
paired_vs_reference <- function(per_rep, keys) {
  ref <- dplyr::filter(per_rep, arm == REFERENCE)
  others <- setdiff(unique(per_rep$arm), REFERENCE)
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

cell_label <- function(cell) {
  if (cell$source == "synthetic") sprintf("%s G=%d", cell$error, cell$G) else
    sprintf("%s G=%d", cell$app, cell$G)
}

summarise_cell <- function(files, cell, n_total) {
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  results <- all[ok]
  # locals first: tibble() evaluates sequentially, a `cell` column would shadow the row
  spec <- list(
    key = cell$key,
    source = cell$source,
    cell = cell$cell,
    label = cell_label(cell),
    error = cell$error,
    app = cell$app,
    G = cell$G
  )
  base <- tibble::as_tibble(spec)
  failures <- dplyr::bind_rows(lapply(all, \(r) {
    if (is.null(r$failures) || !NROW(r$failures)) return(NULL)
    dplyr::mutate(r$failures, key = r$key, rep = r$rep, task_status = r$status)
  }))
  crn <- dplyr::bind_rows(lapply(results, \(r) {
    if (is.null(r$crn)) return(NULL)
    dplyr::mutate(r$crn, key = r$key, rep = r$rep)
  }))
  fits <- dplyr::bind_rows(lapply(results, \(r) {
    dplyr::mutate(r$fits, key = r$key, rep = r$rep, boot_run = r$boot, B = r$B)
  }))
  rows <- list()
  contrasts <- list()
  for (estimand in ESTIMANDS) {
    if (!length(results) || is.null(results[[1]]$records[[estimand]])) next
    keys <- dplyr::mutate(base, estimand = estimand)
    per_rep <- list()
    for (i in seq_len(nrow(ARMS))) {
      m <- arm_matrices(results, estimand, ARMS[i, ])
      if (length(m$reps) == 0) next
      out <- score_arm(m, keys, ARMS[i, ], n_total)
      rows[[length(rows) + 1]] <- out$row
      per_rep[[length(per_rep) + 1]] <- out$per_rep
    }
    contrasts[[estimand]] <- paired_vs_reference(
      dplyr::bind_rows(per_rep),
      keys
    )
  }
  list(
    cells = dplyr::bind_rows(rows),
    contrasts = dplyr::bind_rows(contrasts),
    failures = failures,
    crn = crn,
    fits = fits,
    tasks = tibble::tibble(
      key = cell$key,
      n_files = length(all),
      n_success = sum(ok),
      task_time = sum(vapply(all, \(r) r$time %||% NA_real_, 0), na.rm = TRUE),
      peak_mb = max(vapply(all, \(r) r$peak_mb %||% NA_real_, 0), na.rm = TRUE)
    ),
    software = dplyr::distinct(dplyr::bind_rows(lapply(all, \(r) {
      tibble::as_tibble(r$software)
    })))
  )
}

# Markdown ---------------------------------------------------------------------------------

md_table <- function(df, digits = 1) {
  fmt <- \(x) {
    if (is.na(x)) return("")
    if (is.double(x)) formatC(x, format = "f", digits = digits) else
      as.character(x)
  }
  body <- vapply(
    seq_len(nrow(df)),
    \(i) {
      paste0(
        "| ",
        paste(vapply(df[i, ], \(x) fmt(x[[1]]), ""), collapse = " | "),
        " |"
      )
    },
    ""
  )
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body,
    ""
  )
}

order_cells <- function(df) {
  dplyr::arrange(
    df,
    factor(source, c("synthetic", "plasmode")),
    factor(error, c("iid", "ou", "fpc")),
    factor(app, COMP_APPS),
    dplyr::desc(G)
  )
}

# Wide table of one metric (cells x arms) for one estimand and source.
wide_metric <- function(df, est, src, value, arms = ARMS$arm) {
  df |>
    dplyr::filter(estimand == est, source == src, arm %in% arms) |>
    order_cells() |>
    dplyr::mutate(
      arm = factor(
        ARMS$label[match(arm, ARMS$arm)],
        ARMS$label[ARMS$arm %in% arms]
      )
    ) |>
    dplyr::select(label, arm, value = dplyr::all_of(value)) |>
    tidyr::pivot_wider(names_from = arm, values_from = value)
}

with_mean_row <- function(wide) {
  avg <- dplyr::summarise(
    wide,
    dplyr::across(where(is.numeric), \(x) mean(x, na.rm = TRUE))
  )
  dplyr::bind_rows(wide, dplyr::mutate(avg, label = "mean"))
}

write_summary <- function(out, path, settings) {
  cells_df <- out$cells
  contr <- out$contrasts
  crn_by <- out$crn |>
    dplyr::group_by(column) |>
    dplyr::summarise(
      records = dplyr::n_distinct(paste(key, rep)),
      max_abs = formatC(max(max_abs), format = "e", digits = 1),
      max_rel = formatC(max(max_rel), format = "e", digits = 1),
      .groups = "drop"
    )
  timing <- out$fits |>
    dplyr::left_join(
      dplyr::distinct(cells_df, key, label, source, error, app, G),
      by = "key"
    ) |>
    dplyr::group_by(source, error, app, G, label, arm) |>
    dplyr::summarise(time = mean(time), .groups = "drop") |>
    order_cells() |>
    tidyr::pivot_wider(names_from = arm, values_from = time) |>
    dplyr::select(-source, -error, -app, -G)
  info <- out$fits |>
    dplyr::left_join(
      dplyr::distinct(cells_df, key, label, source, error, app, G),
      by = "key"
    ) |>
    dplyr::mutate(
      npc = as.numeric(sub(
        ".*npc=(\\d+).*",
        "\\1",
        ifelse(grepl("npc=", info), info, NA)
      )),
      nearpd = grepl("nearPD=TRUE", info),
      boot_failed = as.numeric(sub(
        ".*failed=(\\d+).*",
        "\\1",
        ifelse(grepl("failed=", info), info, NA)
      )),
      L_check = as.numeric(sub(
        ".*L_check=([^;]+).*",
        "\\1",
        ifelse(grepl("L_check=", info), info, NA)
      ))
    ) |>
    dplyr::group_by(source, error, app, G, label) |>
    dplyr::summarise(
      `pcre npc median` = median(npc[arm == "PCRE"]),
      `pcre npc range` = if (any(arm == "PCRE"))
        paste(range(npc[arm == "PCRE"]), collapse = "-") else "",
      `GLS raw nearPD (%)` = 100 * mean(nearpd[arm == "GLSraw"]),
      `boot reps` = sum(arm == "BOOT"),
      `boot failed refits` = sum(boot_failed[arm == "BOOT"], na.rm = TRUE),
      `max L check` = suppressWarnings(max(L_check, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    order_cells() |>
    dplyr::select(-source, -error, -app, -G)
  fail_tab <- if (nrow(out$failures)) {
    out$failures |>
      dplyr::mutate(message = substr(gsub("[\n|]", " ", message), 1, 90)) |>
      dplyr::count(key, stage, task_status, message, name = "n")
  } else {
    tibble::tibble(note = "no failures")
  }
  syn_check <- out$study_check
  lines <- c(
    "# Comparator study: refund's remedies for dependent residuals vs REML + CL2",
    "",
    "Generated by `comparators/analyse-comparators.R` from the task records; no number typed by hand.",
    "Nominal 95% pointwise intervals, grid-averaged coverage (%), Gaussian model M1.",
    "Arms: REML model = Vc; REML CL2 = exact CL2, Bayesian (the study's reml_cl2);",
    sprintf(
      "pcre = REML + pcre() curve effect, eigenfunctions from fpca.sc of the REML residual curves (pve %.2f), Vc;",
      PCRE_PVE
    ),
    "GLS raw / fpc = pffr_gls() algorithm with hatSigma = cov() of the REML residual curves / the fpca.sc",
    "covariance + white noise, Vc; boot = curve bootstrap of the REML fit (refit with REML λ),",
    "percentile (pct) or bootstrap-SE Wald interval around the REML estimate.",
    "",
    "## Run and checks",
    "",
    sprintf(
      "- Task files: %d (successful %d); arm failure rows: %d.",
      sum(out$tasks$n_files),
      sum(out$tasks$n_success),
      nrow(out$failures)
    ),
    sprintf("- Settings: %s.", settings),
    sprintf(
      "- Software: refund %s, mgcv %s.",
      paste(unique(substr(out$software$refund_sha, 1, 8)), collapse = ", "),
      paste(unique(out$software$mgcv_version), collapse = ", ")
    ),
    sprintf(
      "- Worker time: %.1f core-hours in total; peak memory per task %d MB.",
      sum(out$tasks$task_time) / 3600,
      as.integer(max(out$tasks$peak_mb))
    ),
    "- CRN: REML columns against the study's stored records (synthetic run 373f6fd65cac, plasmode production results):",
    "",
    md_table(crn_by),
    if (!is.null(syn_check)) {
      sprintf(
        "- REML + CL2 coverage vs `summaries/final/cells.csv` (reml_cl2, synthetic): max |difference| %.2e over %d cell × estimand rows.",
        max(abs(syn_check$coverage - syn_check$study_coverage)),
        nrow(syn_check)
      )
    },
    "",
    "Per-cell diagnostics (pcre: number of FPCs; GLS raw: share of replicates whose hatSigma needed",
    "refund's nearPD repair, condition number > 500; boot: replicates with the arm, failed refits, and",
    "max |difference| between the first refit's own evaluation matrices and the reused REML ones):",
    "",
    md_table(info, digits = 2),
    "Mean fit time per replicate (s; pcre includes the FPCA, REML excludes CL2, boot = all refits):",
    "",
    md_table(timing, digits = 1),
    "Failures (task-level and per arm):",
    "",
    md_table(fail_tab)
  )
  for (src in c("synthetic", "plasmode")) {
    lines <- c(lines, sprintf("## %s cells", tools::toTitleCase(src)), "")
    for (est in ESTIMANDS) {
      cov <- wide_metric(cells_df, est, src, "coverage")
      if (!nrow(cov)) next
      cov <- with_mean_row(dplyr::mutate(
        cov,
        dplyr::across(where(is.numeric), \(x) 100 * x)
      ))
      q05 <- with_mean_row(dplyr::mutate(
        wide_metric(cells_df, est, src, "pointwise_q05"),
        dplyr::across(where(is.numeric), \(x) 100 * x)
      ))
      diffs <- contr |>
        dplyr::filter(estimand == est, source == src) |>
        order_cells() |>
        dplyr::mutate(
          value = sprintf(
            "%+.1f (%.1f)",
            100 * coverage_diff,
            100 * coverage_diff_se
          ),
          arm = factor(ARMS$label[match(arm, ARMS$arm)], ARMS$label)
        ) |>
        dplyr::select(label, arm, value) |>
        tidyr::pivot_wider(names_from = arm, values_from = value)
      ratio <- \(metric) {
        contr |>
          dplyr::filter(estimand == est, source == src) |>
          order_cells() |>
          dplyr::mutate(
            arm = factor(ARMS$label[match(arm, ARMS$arm)], ARMS$label)
          ) |>
          dplyr::select(label, arm, value = dplyr::all_of(metric)) |>
          tidyr::pivot_wider(names_from = arm, values_from = value) |>
          with_mean_row()
      }
      sesd <- with_mean_row(wide_metric(
        cells_df,
        est,
        src,
        "se_sd_ratio",
        setdiff(ARMS$arm, "boot_pct")
      ))
      lines <- c(
        lines,
        sprintf("### %s", est),
        "",
        "Coverage (%):",
        "",
        md_table(cov),
        "5% quantile of pointwise coverage over the grid (%):",
        "",
        md_table(q05),
        "Coverage minus REML + CL2 coverage, paired, percentage points (MC SE):",
        "",
        md_table(diffs),
        "Width relative to REML + CL2 (ratio of mean widths, paired replicates):",
        "",
        md_table(ratio("width_ratio"), digits = 2),
        "Interval score relative to REML + CL2:",
        "",
        md_table(ratio("score_ratio"), digits = 2),
        "MSE of the arm's estimate relative to REML (boot and REML arms share the REML estimate):",
        "",
        md_table(ratio("mse_ratio"), digits = 2),
        "SE/SD (grid mean of sqrt(mean SE²) / SD of the estimate over replicates):",
        "",
        md_table(sesd, digits = 2)
      )
    }
  }
  writeLines(lines, path)
}

# Main ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: analyse-comparators.R results_dir out_dir")
results_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cells <- dplyr::bind_rows(
  comparator_syn_cells() |> dplyr::mutate(app = NA_character_),
  comparator_plas_cells()
)
files <- list.files(
  results_dir,
  "^task-[sp]\\d+-r\\d+\\.rds$",
  full.names = TRUE
)
file_key <- sub("^task-([sp]\\d+)-r.*", "\\1", basename(files))
n_total <- as.integer(Sys.getenv("COMP_MAX_REP", N_REP))
per_cell <- lapply(seq_len(nrow(cells)), \(i) {
  f <- files[file_key == cells$key[i]]
  if (!length(f)) return(NULL)
  summarise_cell(f, cells[i, ], n_total)
})
per_cell <- Filter(Negate(is.null), per_cell)
parts <- c("cells", "contrasts", "failures", "crn", "fits", "tasks", "software")
out <- lapply(
  setNames(nm = parts),
  \(p) dplyr::bind_rows(lapply(per_cell, `[[`, p))
)
out$software <- dplyr::distinct(out$software)
if (!nrow(out$crn)) {
  out$crn <- tibble::tibble(
    estimand = character(),
    column = character(),
    max_abs = numeric(),
    max_rel = numeric(),
    key = character(),
    rep = integer()
  )
}
if (!nrow(out$failures)) {
  out$failures <- tibble::tibble(
    stage = character(),
    message = character(),
    key = character(),
    rep = integer(),
    task_status = character()
  )
}
study_file <- Sys.getenv("STUDY_CELLS", "summaries/final/cells.csv")
if (file.exists(study_file)) {
  study <- read.csv(study_file)
  out$study_check <- out$cells |>
    dplyr::filter(source == "synthetic", arm == REFERENCE, n_rep == N_REP) |>
    dplyr::inner_join(
      dplyr::filter(study, arm == "reml_cl2") |>
        dplyr::select(cell, estimand, study_coverage = coverage_mean),
      by = c("cell", "estimand")
    )
  if (!nrow(out$study_check)) out$study_check <- NULL
}
B_used <- unique(out$fits$B[out$fits$arm == "BOOT"])
boot_reps <- out$fits |>
  dplyr::filter(arm == "BOOT") |>
  dplyr::left_join(dplyr::select(cells, key, source), by = "key") |>
  dplyr::group_by(source) |>
  dplyr::summarise(max_rep = max(rep), .groups = "drop")
settings <- sprintf(
  "replicates per cell up to %d; bootstrap B = %s, run on replicates %s",
  n_total,
  paste(B_used, collapse = "/"),
  paste(
    sprintf("1-%d (%s)", boot_reps$max_rep, boot_reps$source),
    collapse = ", "
  )
)
write.csv(
  out$cells,
  file.path(out_dir, "comparators-cells.csv"),
  row.names = FALSE
)
write.csv(
  out$contrasts,
  file.path(out_dir, "comparators-contrasts.csv"),
  row.names = FALSE
)
write.csv(
  out$failures,
  file.path(out_dir, "comparators-failures.csv"),
  row.names = FALSE
)
write.csv(out$crn, file.path(out_dir, "comparators-crn.csv"), row.names = FALSE)
fit_cols <- c(
  "key",
  "rep",
  "arm",
  "time",
  "converged",
  "n_warnings",
  "info",
  "boot_run",
  "B"
)
write.csv(
  out$fits[, fit_cols],
  file.path(out_dir, "comparators-fits.csv"),
  row.names = FALSE
)
write_summary(out, file.path(out_dir, "SUMMARY.md"), settings)
cat("done:", nrow(out$cells), "rows\n")
