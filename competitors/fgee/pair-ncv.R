# Paired comparison of fastFGEE (experiment B) with the study's own REML and
# curve-blocked NCV fits on the SAME data sets (simulate(cell, rep) is deterministic).
#   Rscript competitors/fgee/pair-ncv.R <study results dir> <fgee B dir> <out csv>
# Also checks common random numbers: the REML estimates of both runs must agree.
args <- commandArgs(TRUE)
study_dir <- args[1]
b_dir <- args[2]
prefix <- if (length(args) >= 4) args[4] else "B"
out_csv <- args[3]

study_files <- list.files(
  study_dir,
  pattern = "^task-.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
b_files <- list.files(
  b_dir,
  pattern = paste0("^", prefix, "-c\\d+-r\\d+\\.rds$"),
  full.names = TRUE
)
key <- \(cell, rep) paste(cell, rep)
b_keys <- vapply(
  b_files,
  \(f) {
    m <- regmatches(
      basename(f),
      regexec(paste0(prefix, "-c(\\d+)-r(\\d+)"), basename(f))
    )[[1]]
    key(as.integer(m[2]), as.integer(m[3]))
  },
  ""
)
wanted <- setNames(b_files, b_keys)

rows <- list()
seen <- character()
for (f in study_files) {
  x <- tryCatch(readRDS(f), error = \(e) NULL)
  if (is.null(x) || !identical(x$status, "success")) next
  k <- key(x$cell, x$rep)
  if (!k %in% names(wanted) || k %in% seen) next
  seen <- c(seen, k)
  b <- readRDS(wanted[[k]])
  if (!is.null(b[["error", exact = TRUE]])) next
  for (est in c("beta", "alpha", "gamma")) {
    s <- x$records[[est]]
    r <- b$records[b$records$estimand == est, ]
    truth <- s[, "truth"]
    get <- \(m) r$est[r$method == m]
    reml_b <- get("reml_cl2")
    rows[[length(rows) + 1]] <- data.frame(
      cell = x$cell,
      rep = x$rep,
      estimand = est,
      crn_max_diff = max(abs(reml_b - s[, "est_REML"])),
      truth_max_diff = max(abs(r$truth[r$method == "reml_cl2"] - truth)),
      se_reml = mean((s[, "est_REML"] - truth)^2),
      se_ncv = mean((s[, "est_NCV"] - truth)^2),
      se_fgee_ind = mean((get("fgee_independent") - truth)^2),
      se_fgee_ar1 = mean((get("fgee_ar1") - truth)^2),
      se_fgee_fpca = mean((get("fgee_fpca") - truth)^2)
    )
  }
}
out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
cat(
  "paired data sets:",
  length(seen),
  "; max REML difference (CRN check):",
  signif(max(out$crn_max_diff), 3),
  "; max truth difference:",
  signif(max(out$truth_max_diff), 3),
  "\n"
)
