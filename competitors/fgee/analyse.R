# Summaries of the fastFGEE comparison (run.R output).
#   Rscript competitors/fgee/analyse.R <results dir with A/ and B/> <out dir>
# Per replicate: grid-average coverage, mean width, squared error per (method, estimand);
# then means over replicates with Monte Carlo SEs (replicates are independent data
# sets; points within a replicate are not). Also the SE calibration ratio
# mean(SE) / sd(estimate) per point, averaged over the grid.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
args <- commandArgs(TRUE)
in_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summarise_records <- function(records) {
  per_rep <- records |>
    mutate(covered = truth >= lower & truth <= upper, width = upper - lower) |>
    group_by(cell, method, estimand, rep) |>
    summarise(
      coverage = mean(covered),
      width = mean(width),
      sq_err = mean((est - truth)^2),
      .groups = "drop"
    )
  cover <- per_rep |>
    group_by(cell, method, estimand) |>
    summarise(
      R = n(),
      coverage_mc_se = sd(coverage) / sqrt(n()),
      coverage = mean(coverage),
      width = mean(width),
      rmse = sqrt(mean(sq_err)),
      .groups = "drop"
    )
  calib <- records |>
    group_by(cell, method, estimand, point) |>
    summarise(
      mean_se = mean(se),
      sd_est = sd(est),
      bias = mean(est - truth),
      .groups = "drop"
    ) |>
    group_by(cell, method, estimand) |>
    summarise(
      se_sd_ratio = mean(mean_se / sd_est),
      rel_bias = sqrt(mean(bias^2)) / sqrt(mean(sd_est^2)),
      .groups = "drop"
    )
  list(
    per_rep = per_rep,
    table = left_join(cover, calib, by = c("cell", "method", "estimand"))
  )
}

# Files are read and summarised one cell at a time (all records at once do not fit in
# memory for B).
read_cell <- function(files) {
  out <- lapply(files, readRDS)
  failed <- vapply(out, \(x) !is.null(x[["error", exact = TRUE]]), TRUE)
  errors <- lapply(out[failed], \(x) {
    data.frame(
      cell = x$cell$cell,
      rep = x$rep,
      error = x[["error", exact = TRUE]]
    )
  })
  method_errors <- lapply(out[!failed], \(x) {
    e <- x$errors[nzchar(x$errors)]
    if (length(e))
      data.frame(cell = x$cell$cell, rep = x$rep, method = names(e), error = e)
  })
  recs <- lapply(out[!failed], \(x) {
    r <- x$records
    r$cell <- x$cell$cell
    r$rep <- x$rep
    r
  })
  cells <- bind_rows(lapply(out, \(x) {
    cl <- x$cell
    as.data.frame(cl[vapply(cl, length, 1L) == 1])
  })) |>
    distinct()
  timing <- bind_rows(lapply(out[!failed], \(x) {
    data.frame(
      cell = x$cell$cell,
      rep = x$rep,
      part = names(x$timing),
      time = x$timing
    )
  }))
  list(
    records = bind_rows(recs),
    cells = cells,
    errors = bind_rows(errors),
    method_errors = bind_rows(method_errors),
    timing = timing
  )
}

read_experiment <- function(experiment) {
  files <- list.files(
    file.path(in_dir, experiment),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  cell_of <- as.integer(sub(".*-c(\\d+)-r\\d+\\.rds$", "\\1", files))
  parts <- lapply(split(files, cell_of), \(f) {
    x <- read_cell(f)
    x$summary <- if (nrow(x$records)) summarise_records(x$records)
    x$records <- NULL
    x
  })
  pick <- \(name) bind_rows(lapply(parts, `[[`, name))
  list(
    per_rep = bind_rows(lapply(parts, \(x) x$summary$per_rep)),
    table = bind_rows(lapply(parts, \(x) x$summary$table)),
    cells = pick("cells"),
    errors = pick("errors"),
    method_errors = pick("method_errors"),
    timing = pick("timing"),
    n_tasks = length(files)
  )
}

for (experiment in c("A", "B", "P")) {
  if (!dir.exists(file.path(in_dir, experiment))) next
  x <- read_experiment(experiment)
  if (!nrow(x$table)) next
  s <- list(per_rep = x$per_rep, table = x$table)
  tab <- left_join(s$table, x$cells, by = "cell")
  write.csv(
    tab,
    file.path(out_dir, paste0(experiment, "-summary.csv")),
    row.names = FALSE
  )
  if (experiment == "A") {
    # Their Table 3 statistic: coverage pooled over beta0, beta1, beta2 per replicate.
    pooled <- s$per_rep |>
      group_by(cell, method, rep) |>
      summarise(coverage = mean(coverage), .groups = "drop") |>
      group_by(cell, method) |>
      summarise(
        R = n(),
        coverage_mc_se = sd(coverage) / sqrt(n()),
        coverage = mean(coverage),
        .groups = "drop"
      ) |>
      left_join(x$cells, by = "cell")
    write.csv(pooled, file.path(out_dir, "A-pooled.csv"), row.names = FALSE)
  }
  write.csv(
    x$errors,
    file.path(out_dir, paste0(experiment, "-errors.csv")),
    row.names = FALSE
  )
  write.csv(
    x$method_errors,
    file.path(out_dir, paste0(experiment, "-method-errors.csv")),
    row.names = FALSE
  )
  timing <- x$timing |>
    group_by(cell, part) |>
    summarise(median_s = median(time), .groups = "drop")
  write.csv(
    timing,
    file.path(out_dir, paste0(experiment, "-timing.csv")),
    row.names = FALSE
  )
  cat(
    experiment,
    ":",
    x$n_tasks,
    "tasks,",
    nrow(x$errors),
    "task errors,",
    nrow(x$method_errors),
    "method errors\n"
  )
}
