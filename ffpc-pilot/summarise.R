# Summary of the ffpc pilot: Rscript ffpc-pilot/summarise.R <dir> [<dir> ...] <out.csv>
# Per cell and pve: retained PCs, relative RMSE of beta (sqrt(sum sq. error / sum truth^2) over
# replicates and grid), pointwise 95% coverage of the true beta and of its projection onto the FPC
# span, half-width relative to ff REML + CL2, and detection (share of points with |truth| > 25% of
# max |truth| whose interval excludes zero on the truth's side). ff arms from the study records.
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out_csv <- args[length(args)]
files <- unlist(lapply(
  args[-length(args)],
  list.files,
  pattern = "^ffpc-[pk][0-9]+.*\\.rds$",
  full.names = TRUE
))
z <- qnorm(0.975)
rows <- lapply(files, \(f) {
  r <- readRDS(f)
  # target pve (or fixed K) from the file name: ffpc-p<100 pve>-... or ffpc-k<K>-...
  tag <- regmatches(basename(f), regexpr("^ffpc-[pk][0-9]+", basename(f)))
  r$pve_target <- sub("^ffpc-", "", tag)
  if (!is.null(r$error))
    return(tibble(cell = r$cell, rep = r$rep, error = r$error))
  m <- r$records
  tr <- m[, "truth"]
  big <- abs(tr) > 0.25 * max(abs(tr))
  arm <- \(name, est, se, target = tr) {
    lo <- est - z * se
    hi <- est + z * se
    tibble(
      arm = name,
      sq_err = sum((est - tr)^2),
      sq_err_proj = sum((est - m[, "truth_proj"])^2),
      sq_truth = sum(tr^2),
      cover = mean(lo <= target & target <= hi),
      cover_proj = mean(lo <= m[, "truth_proj"] & m[, "truth_proj"] <= hi),
      halfwidth = mean(z * se),
      detect = mean(((lo > 0 & tr > 0) | (hi < 0 & tr < 0))[big])
    )
  }
  bind_rows(
    arm("ffpc_cl2", m[, "est_ffpc"], m[, "se_cl2_ffpc"]),
    arm("ffpc_model", m[, "est_ffpc"], m[, "se_model_ffpc"]),
    if ("est_REML" %in% colnames(m))
      arm("ff_reml_cl2", m[, "est_REML"], m[, "se_cl2_REML"]),
    if ("est_NCV" %in% colnames(m))
      arm("ff_ncv_cl2", m[, "est_NCV"], m[, "se_cl2_NCV"])
  ) |>
    mutate(
      cell = r$cell,
      rep = r$rep,
      pve = r$pve_target,
      method = r$method %||% "REML",
      K = r$K,
      proj_err = sqrt(mean((m[, "truth_proj"] - tr)^2) / mean(tr^2))
    )
})
d <- bind_rows(rows)
if ("error" %in% names(d)) print(filter(d, !is.na(error)) |> count(cell, error))
d <- filter(d, !is.na(arm))
s <- d |>
  group_by(cell, pve, method, arm) |>
  summarise(
    n_rep = n(),
    K_median = median(K),
    proj_rel_err = mean(proj_err),
    rel_rmse = sqrt(sum(sq_err) / sum(sq_truth)),
    rel_rmse_vs_proj = sqrt(sum(sq_err_proj) / sum(sq_truth)),
    coverage = mean(cover),
    coverage_mcse = sd(cover) / sqrt(n()),
    coverage_proj = mean(cover_proj),
    halfwidth = mean(halfwidth),
    detect = mean(detect),
    .groups = "drop"
  ) |>
  group_by(cell, pve, method) |>
  mutate(
    hw_vs_reml = if (any(arm == "ff_reml_cl2")) {
      halfwidth / halfwidth[arm == "ff_reml_cl2"]
    } else {
      NA_real_
    }
  ) |>
  ungroup()
readr::write_csv(s, out_csv)
print(as.data.frame(mutate(s, across(where(is.double), \(x) signif(x, 3)))))
