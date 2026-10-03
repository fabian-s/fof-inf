# Markdown tables for the fastFGEE comparison, from analyse.R's CSVs.
#   Rscript competitors/fgee/report.R <summary dir>   (writes <summary dir>/REPORT.md)
# All numbers in REPORT.md are computed here; prose elsewhere cites this file.
suppressPackageStartupMessages(library(dplyr))
dir <- commandArgs(TRUE)[1]
out <- character()
say <- \(...) out <<- c(out, paste0(...))
fmt <- \(x, d = 3) formatC(x, format = "f", digits = d)
range_txt <- \(x, d = 2) paste0(fmt(min(x), d), "–", fmt(max(x), d))

# kable-free markdown table: rows = row labels, columns = methods
md_table <- function(df, row, col, value, digits = 3) {
  wide <- tidyr::pivot_wider(
    df[, c(row, col, value)],
    names_from = all_of(col),
    values_from = all_of(value)
  )
  wide <- as.data.frame(wide)
  num <- vapply(wide, is.numeric, TRUE)
  wide[num] <- lapply(wide[num], fmt, d = digits)
  c(
    paste0("| ", paste(names(wide), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(wide)), collapse = "|"), "|"),
    apply(wide, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
}

# Experiment A ------------------------------------------------------------------------
if (file.exists(file.path(dir, "A-pooled.csv"))) {
  a <- read.csv(file.path(dir, "A-pooled.csv"))
  s <- read.csv(file.path(dir, "A-summary.csv"))
  a$design <- sprintf("%s N=%d n_i=%d", a$family, a$N, a$n_i)
  s$design <- sprintf("%s N=%d n_i=%d", s$family, s$N, s$n_i)
  methods <- c(
    "model",
    "fgee_base",
    "fgee_base_wild",
    "fgee_base_fixed",
    "refund_cr_freq_centred",
    "refund_cr",
    "refund_cl2",
    "reml_cl2",
    "fgee_independent_independent",
    "fgee_independent_independent_wild",
    "fgee_ar1_exchangeable",
    "fgee_ar1_exchangeable_wild"
  )
  a <- a[a$method %in% methods, ]
  a$method <- factor(a$method, methods)
  a <- a[order(a$family, a$n_i, a$N, a$method), ]
  say("# fastFGEE comparison: results\n")
  say(
    "## A. fastFGEE's simulation design (Loewinger et al., arXiv 2506.20437v3, §4)\n"
  )
  say(sprintf(
    "Replicates per cell: %d–%d. Pointwise 95%% coverage pooled over beta0, beta1, beta2 (their Table 3 statistic).\n",
    min(a$R),
    max(a$R)
  ))
  out <- c(out, md_table(a, "design", "method", "coverage"), "")
  say(sprintf(
    "Largest Monte Carlo SE of a pooled coverage: %s.\n",
    fmt(max(a$coverage_mc_se))
  ))
  rng <- \(m, fam = NULL, n_i = NULL) {
    x <- a[a$method == m, ]
    if (!is.null(fam)) x <- x[x$family %in% fam, ]
    if (!is.null(n_i)) x <- x[x$n_i %in% n_i, ]
    range_txt(x$coverage)
  }
  say("Ranges of pooled coverage over cells (n_i = 5 | n_i = 1):\n")
  for (m in methods) {
    say(sprintf("- %s: %s | %s", m, rng(m, n_i = 5), rng(m, n_i = 1)))
  }
  say("")
  say(
    "Per family, n_i = 5 (their design): fgee_base Gaussian ",
    rng("fgee_base", "gaussian", 5),
    ", Poisson ",
    rng("fgee_base", "poisson", 5),
    ", binomial ",
    rng("fgee_base", "binomial", 5),
    ".\n"
  )
  s <- s[s$method %in% methods, ]
  s$method <- factor(s$method, methods)
  for (est in c("beta1", "beta2")) {
    say(sprintf("### SE calibration mean(SE)/sd(estimate), %s\n", est))
    x <- s[s$estimand == est, ]
    x <- x[order(x$family, x$n_i, x$N, x$method), ]
    out <- c(out, md_table(x, "design", "method", "se_sd_ratio", 2), "")
  }
  say("### RMSE relative to the pffr fit (their Table 1 statistic), beta1\n")
  x <- s[s$estimand == "beta1", ]
  base <- x[x$method == "refund_cl2", c("design", "rmse")]
  x <- merge(x, setNames(base, c("design", "rmse_pffr")), by = "design")
  x$rmse_ratio <- x$rmse / x$rmse_pffr
  x <- x[
    x$method %in%
      c("reml_cl2", "fgee_independent_independent", "fgee_ar1_exchangeable"),
  ]
  out <- c(out, md_table(x, "design", "method", "rmse_ratio", 2), "")
}

# Experiments B and P -------------------------------------------------------------------
bp_report <- function(file, label, design_fun) {
  if (!file.exists(file)) return(invisible())
  x <- read.csv(file)
  x$design <- design_fun(x)
  methods <- c(
    "reml_model",
    "reml_cl2",
    "fgee_independent",
    "fgee_independent_wild",
    "fgee_ar1",
    "fgee_ar1_wild",
    "fgee_fpca",
    "fgee_fpca_wild"
  )
  x <- x[x$method %in% methods, ]
  x$method <- factor(x$method, methods)
  x <- x[order(x$design, x$method), ]
  say(label, "\n")
  say(sprintf("Replicates per cell: %d–%d.\n", min(x$R), max(x$R)))
  for (est in c("beta", "mean", "alpha", "gamma")) {
    y <- x[x$estimand == est, ]
    if (!nrow(y)) next
    say(sprintf("### Coverage, %s\n", est))
    out <<- c(out, md_table(y, "design", "method", "coverage"), "")
  }
  for (est in c("beta", "mean")) {
    y <- x[x$estimand == est, ]
    ref <- y[y$method == "reml_cl2", c("design", "rmse", "width")]
    y <- merge(
      y,
      setNames(ref, c("design", "rmse_reml", "width_cl2")),
      by = "design"
    )
    y$rmse_ratio <- y$rmse / y$rmse_reml
    y$width_ratio <- y$width / y$width_cl2
    y <- y[y$method %in% c("fgee_independent", "fgee_ar1", "fgee_fpca"), ]
    say(sprintf("### %s: RMSE ratio to REML\n", est))
    out <<- c(out, md_table(y, "design", "method", "rmse_ratio", 2), "")
    say(sprintf("### %s: interval width ratio to REML + CL2\n", est))
    out <<- c(out, md_table(y, "design", "method", "width_ratio", 2), "")
    z <- x[x$estimand == est, ]
    say(sprintf("### %s: mean(SE)/sd(estimate)\n", est))
    out <<- c(out, md_table(z, "design", "method", "se_sd_ratio", 2), "")
  }
  say(sprintf(
    "Largest Monte Carlo SE of a coverage: %s.\n",
    fmt(max(x$coverage_mc_se))
  ))
}
bp_report(
  file.path(dir, "B-summary.csv"),
  "## B. Restart study data sets (M1: ff + z, one curve per subject, core cells at mid signal)",
  \(x) sprintf("%s %s G=%d", x$family, x$error, x$G)
)
bp_report(
  file.path(dir, "P-summary.csv"),
  "## P. Plasmode (public datasets, MID truth, REML residual curves attached)",
  \(x) sprintf("%s G=%d", x$app, x$G)
)
writeLines(out, file.path(dir, "REPORT.md"))
cat("wrote", file.path(dir, "REPORT.md"), "\n")
