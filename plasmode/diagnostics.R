# Per-dataset diagnostics of the plasmode datasets, from
# the truth caches (built on demand): the survey metrics regenerated on the fitted model,
# subject frame and grid; per truth source EDF, log sp, signal and residual variance and
# the MID distances; per residual source DE, variance profile and leading eigenvalues;
# the test-cohort check; block partitions. Aggregates only (no per-subject values).
# Run from the repository root:  Rscript plasmode/diagnostics.R [<app> ...]
# Writes PLASMODE_SUMMARIES/diagnostics.csv (default summaries/plasmode).
source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
suppressPackageStartupMessages(library(refund))

N_EIGEN <- 3 # leading covariance eigenvalue shares reported
N_SCORES <- 3 # residual FPC scores for the between-block correlation
FPC_SHARE <- 0.9
AEMET_NEAR_KM <- 250

# Residual-curve metrics (those of the dataset screening) of an error matrix E (G x D).
error_metrics <- function(E, yind) {
  S <- crossprod(E) / nrow(E) # second moments: sign-flipped errors have mean zero
  R <- stats::cov2cor(S)
  off <- R[upper.tri(R)]
  D <- length(yind)
  lag_at <- \(share) {
    k <- max(1, round(share * (D - 1)))
    mean(R[cbind(seq_len(D - k), seq_len(D - k) + k)])
  }
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0) / sum(pmax(ev, 0))
  sds <- sqrt(diag(S))
  c(
    residual_variance = mean(diag(S)),
    DE = sum(S) / sum(diag(S)),
    lag1 = mean(R[cbind(seq_len(D - 1), seq_len(D - 1) + 1)]),
    lag10 = lag_at(0.1),
    lag25 = lag_at(0.25),
    neg_share = mean(off < 0),
    min_r = min(off),
    sd_ratio = max(sds) / min(sds),
    fpc90 = which(cumsum(ev) >= FPC_SHARE)[1],
    setNames(ev[seq_len(N_EIGEN)], paste0("eigen_share", seq_len(N_EIGEN))),
    variance_q10 = unname(stats::quantile(diag(S), 0.1)),
    variance_q90 = unname(stats::quantile(diag(S), 0.9))
  )
}

# Within- and between-block correlation of the residual FPC scores: mean normalised
# score products over pairs (ordered data: adjacent curves; aemet: stations < 250 km).
block_metrics <- function(E, d, app) {
  if (is.null(d$block)) return(NULL)
  scores <- stats::prcomp(E, center = FALSE)$x[, seq_len(N_SCORES)]
  scores <- sweep(scores, 2, sqrt(colMeans(scores^2)), "/")
  G <- nrow(E)
  pairs <- if (app == "aemet") {
    aemet_near_pairs(d)
  } else {
    cbind(seq_len(G - 1), seq_len(G - 1) + 1)
  }
  same <- d$block[pairs[, 1]] == d$block[pairs[, 2]]
  product <- rowMeans(scores[pairs[, 1], ] * scores[pairs[, 2], ])
  sizes <- table(d$block)
  c(
    n_blocks = length(sizes),
    block_size_min = min(sizes),
    block_size_max = max(sizes),
    score_corr_within = mean(product[same]),
    score_corr_between = mean(product[!same]),
    n_pairs_between = sum(!same),
    n_blocks_alt = if (!is.null(d$block_alt)) length(unique(d$block_alt))
  )
}

aemet_near_pairs <- function(d) {
  utils::data("aemet", package = "fda.usc", envir = environment())
  stopifnot(identical(make.unique(aemet$df$ind), d$ids))
  lat <- aemet$df$latitude * pi / 180
  lon <- aemet$df$longitude * pi / 180
  a <- outer(lat, lat, \(p, q) sin((q - p) / 2)^2) +
    outer(cos(lat), cos(lat)) * outer(lon, lon, \(p, q) sin((q - p) / 2)^2)
  km <- 2 * 6371 * asin(pmin(sqrt(a), 1))
  which(upper.tri(km) & km < AEMET_NEAR_KM, arr.ind = TRUE)
}

rms <- \(x) sqrt(mean(x^2))

app_diagnostics <- function(app) {
  cells <- plasmode_cells()
  mine <- cells[cells$app == app, ]
  truths <- lapply(
    setNames(nm = intersect(c("REML", "NCV", "MID"), mine$truth)),
    \(s) {
      plasmode_truth(app, s)
    }
  )
  d <- truths$REML$data
  row <- \(source, values) {
    tibble::tibble(
      app = app,
      source = source,
      metric = names(values),
      value = unname(values)
    )
  }
  out <- list()
  check <- test_cohort_check(d$X, d$z)
  # the fitted model (REML truth) on the plasmode frame and grid
  Y_centred <- sweep(d$Y, 2, colMeans(d$Y))
  out$model <- row(
    "model (REML)",
    c(
      G = length(d$z),
      D = length(d$yind),
      n_eval_t = length(d$eval_t),
      R2 = 1 - sum(truths$REML$residuals^2) / sum(Y_centred^2),
      rank_X = centred_rank(d$X),
      test_cohort_p = check$p,
      test_cohort_one_sided = check$one_sided,
      test_cohort_first50 = identical(
        d$test_rows %||% seq_len(N_TEST),
        seq_len(N_TEST)
      ),
      setNames(check$smd, paste0("test_smd_", names(check$smd))),
      mu_diff_share = if (!is.null(truths$NCV)) {
        sum((truths$REML$mu - truths$NCV$mu)^2) / sum(truths$REML$residuals^2)
      }
    )
  )
  for (s in names(truths)) {
    tr <- truths[[s]]
    terms <- tr$terms %||% truth_fit_info(app, s)$terms
    sp <- tr$meta$sp %||% truth_fit_info(app, s)$sp
    mu_centred <- sweep(tr$mu, 2, colMeans(tr$mu))
    values <- c(
      edf_total = sum(terms$edf),
      setNames(terms$edf, paste0("edf:", terms$term)),
      setNames(log(sp), paste0("log_sp:", names(sp))),
      signal_variance = mean(mu_centred^2),
      mid_interior = if (s == "MID") as.numeric(tr$mid$interior)
    )
    if (s == "MID") {
      for (other in intersect(c("REML", "NCV"), names(truths))) {
        values <- c(
          values,
          vapply(
            setNames(
              names(tr$targets),
              paste0("rms_to_", other, ":", names(tr$targets))
            ),
            \(e) rms(tr$targets[[e]] - truths[[other]]$targets[[e]]),
            numeric(1)
          )
        )
      }
    }
    out[[paste("truth", s)]] <- row(paste("truth", s), values)
  }
  for (r in unique(mine$residual)) {
    E <- plasmode_errors(app, r, truths$REML)
    values <- error_metrics(E, d$yind)
    if (startsWith(r, "beat"))
      values <- c(
        values,
        c_scale = if (r == "beat_m") plasmode_beats(app)$c else 1
      )
    out[[paste("residual", r)]] <- row(paste("residual", r), values)
  }
  blocks <- block_metrics(truths$REML$residuals, d, app)
  if (!is.null(blocks)) out$blocks <- row("blocks (REML residuals)", blocks)
  dplyr::bind_rows(out)
}

apps <- commandArgs(trailingOnly = TRUE)
if (!length(apps)) apps <- names(PLASMODE_EXTENSION)
diagnostics <- dplyr::bind_rows(lapply(apps, app_diagnostics))
out_dir <- Sys.getenv("PLASMODE_SUMMARIES", "summaries/plasmode")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(diagnostics, file.path(out_dir, "diagnostics.csv"))
cat("diagnostics for", toString(apps), "written to", out_dir, "\n")
