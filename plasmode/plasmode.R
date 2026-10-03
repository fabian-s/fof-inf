# Plasmode applications: real covariates and real residual curves around a known truth.
# The model is fitted once to the real data (by REML and by curve-blocked NCV, and at the
# log-midpoint of their smoothing parameters, "MID"); a fit's means are frozen as the
# truth, and replicate data are
#   Y*_i(t) = mu_i(t) + w_i e_{pi(i)}(t),  w_i = +-1,
# with pi the identity ("attached": each subject keeps its own residual curve) or a random
# permutation ("detached"), e the residual curves of a named fit (or, for the ECG study,
# its beat-to-beat differences), and w drawn per curve or per block of curves.
# Everything else (fits, intervals, scoring) is the study's code.
# Needs scripts/setup.R and plasmode/sources.R sourced first.

PLASMODE_APPS <- c("ecg_gls4", "running") # the original 16 cells
PLASMODE_G_SMALL <- 40L # subsample of subjects per replicate
PLASMODE_CACHE_SCHEMA <- 1L # version of the truth-cache metadata
MID_EDF_EDGE <- 0.01 # EDF share of a term's penalized range that counts as its boundary

# Main design: per app its G (all subjects), fitted D, residual
# sources (crossed with the truths REML, NCV, MID) or explicit truth:residual pairs, the
# (truth:residual) cells at G = all with an AR(1) fit, and the block arm.
PLASMODE_EXTENSION <- list(
  ecg_gls4b = list(
    G = 78L,
    D = 61L,
    residuals = c("REML", "NCV", "beat_m", "beat_1"),
    ar1 = c("MID:REML", "MID:beat_m")
  ),
  af = list(G = 118L, D = 30L),
  running = list(
    G = 90L,
    D = 61L,
    pairs = c("REML:NCV", "NCV:REML", "MID:REML", "MID:NCV")
  ),
  dti = list(G = 92L, D = 55L),
  gait_vc = list(G = 138L, D = 101L),
  ecg_roahd = list(G = 100L, D = 101L, basis = "ecg8"),
  ocean = list(G = 116L, D = 101L, block = TRUE),
  aemet = list(G = 73L, D = 73L, block = TRUE, block_alt = TRUE),
  adelaide = list(G = 102L, D = 48L, block = TRUE)
)

# Design ------------------------------------------------------------------------------------

# The plasmode design: the original 16 cells first and unchanged (ids 1-16), then the
# extension's cells appended.
plasmode_cells <- function(extension = PLASMODE_EXTENSION) {
  legacy <- legacy_plasmode_cells()
  new <- purrr::imap(extension, extension_cells) |> dplyr::bind_rows()
  new$cell <- nrow(legacy) + seq_len(nrow(new))
  cells <- dplyr::bind_rows(legacy, new)
  dplyr::relocate(cells, cell)
}

# The original design: app x truth source x error scheme x G (all subjects, 40); residuals
# from the truth's own fit, signs per curve.
legacy_plasmode_cells <- function() {
  sizes <- c(ecg_gls4 = 79L, running = 90L) # subjects after cleaning (plasmode_data)
  cells <- purrr::map(PLASMODE_APPS, \(app) {
    design_block(
      "plasmode",
      truth = c("REML", "NCV"), # which fit's means are the truth
      error = c("attached", "detached"),
      G = c(sizes[[app]], PLASMODE_G_SMALL),
      signal = "real"
    ) |>
      dplyr::mutate(app = app)
  }) |>
    dplyr::bind_rows()
  # AR(1) comparison on real errors: all subjects, own residuals.
  ar1 <- cells$error == "attached" & cells$G > PLASMODE_G_SMALL
  cells$extra[ar1] <- "ar1"
  cells$fit_family <- cells$family
  cells$residual <- cells$truth
  cells$flip <- "curve"
  cells$cell <- seq_len(nrow(cells))
  dplyr::relocate(cells, cell)
}

# One app's extension cells: truth x residual pairs x G {all, 40} with curve flips, plus
# (block arm) the diagonal pairs at G = all with block flips. Extra fits: AR(1) in the
# listed G = all cells; the block-NCV fit in every G = all cell of a block app; the
# alternative partition's NCV fit in its block-flip cells.
extension_cells <- function(spec, app) {
  truths <- spec$truths %||% c("REML", "NCV", "MID")
  pairs <- spec$pairs %||%
    with(
      tidyr::expand_grid(
        truth = truths,
        residual = spec$residuals %||% c("REML", "NCV")
      ),
      paste(truth, residual, sep = ":")
    )
  rows <- tidyr::expand_grid(pair = pairs, G = c(spec$G, PLASMODE_G_SMALL))
  rows$flip <- "curve"
  if (isTRUE(spec$block)) {
    rows <- dplyr::bind_rows(
      rows,
      tibble::tibble(
        pair = c("REML:REML", "NCV:NCV"),
        G = spec$G,
        flip = "block"
      )
    )
  }
  parts <- strsplit(rows$pair, ":", fixed = TRUE)
  rows$truth <- vapply(parts, `[`, "", 1)
  rows$residual <- vapply(parts, `[`, "", 2)
  all_subjects <- rows$G == spec$G
  extra <- \(name, on) ifelse(on, name, "")
  rows$extra <- trimws(paste(
    extra("ar1", all_subjects & rows$pair %in% (spec$ar1 %||% "MID:REML")),
    extra("block", all_subjects & isTRUE(spec$block)),
    extra("block_alt", rows$flip == "block" & isTRUE(spec$block_alt))
  ))
  rows$extra <- gsub(" +", " ", rows$extra)
  defaults <- design_block(
    "plasmode",
    error = "attached",
    signal = "real",
    D = spec$D
  )
  defaults <- defaults[setdiff(names(defaults), names(rows))]
  cells <- dplyr::bind_cols(
    defaults[rep(1, nrow(rows)), ],
    dplyr::select(rows, -pair)
  )
  cells$app <- app
  cells$basis <- spec$basis %||% "default"
  cells$fit_family <- cells$family
  cells
}

# Truths ------------------------------------------------------------------------------------

truth_file <- function(app, source) {
  file.path(PLASMODE_CACHE, paste0("truth-", app, "-", source, ".rds"))
}

# Loaded caches, per R process: a worker checks each cache once.
.plasmode_memo <- new.env()

# The frozen truth of one app and truth source: fitted means, residual curves, the
# estimand targets on the test cohort, per-term EDF and sp, and metadata. Cached outside
# the repository (per-subject curves). A cache with metadata is checked against the
# current data, grids, model and software and refused if stale; caches without metadata
# (the original 16 cells') are accepted only for the original apps' REML/NCV truths.
plasmode_truth <- function(app, source) {
  file <- truth_file(app, source)
  if (!is.null(.plasmode_memo[[file]])) return(.plasmode_memo[[file]])
  if (file.exists(file)) {
    truth <- readRDS(file)
    legacy <- app %in% PLASMODE_LEGACY_APPS && source %in% c("REML", "NCV")
    if (is.null(truth$meta) && !legacy) {
      stop("truth cache without metadata for ", app, "/", source, ": ", file)
    }
    if (!is.null(truth$meta)) check_truth_meta(truth, app, source)
  } else {
    truth <- build_truth(app, source)
    dir.create(PLASMODE_CACHE, recursive = TRUE, showWarnings = FALSE)
    Sys.chmod(PLASMODE_CACHE, "0700")
    atomic_save(truth, file) # concurrent workers may build the same (identical) truth
  }
  assign(file, truth, envir = .plasmode_memo)
  truth
}

PLASMODE_FIT_CELL <- list(model = "M1", fit_family = "gaussian")
# The model an app's truths are fitted with: its cells' bases (ecg_roahd: wider intercept).
plasmode_fit_cell <- function(app) {
  c(
    PLASMODE_FIT_CELL,
    list(basis = PLASMODE_EXTENSION[[app]]$basis %||% "default")
  )
}

build_truth <- function(app, source) {
  d <- plasmode_data(app)
  sim <- plasmode_sim(d, d$Y)
  mid <- NULL
  fit <- if (source == "MID") {
    mid <- mid_fit(app, sim)
    mid$fit
  } else {
    fit_model(sim, plasmode_fit_cell(app), source)
  }
  truth <- truth_from_fit(fit, d, sim)
  truth$meta <- truth_meta(app, source, d, fit)
  truth$mid <- mid[c("interior", "terms", "sp_reml", "sp_ncv")]
  truth
}

# Fitted means, residual curves, targets and per-term EDF of one fit of the real data.
truth_from_fit <- function(fit, d, sim) {
  lp <- predict(fit, type = "lpmatrix", reformat = FALSE) # rows curve-major
  mu <- matrix(drop(lp %*% fit$coefficients), nrow(d$Y), byrow = TRUE)
  targets <- lapply(setNames(nm = model_estimands("M1")), \(e) {
    drop(evaluation_matrix(fit, sim, e) %*% fit$coefficients)
  })
  list(
    data = d,
    mu = mu,
    residuals = d$Y - mu,
    targets = targets,
    terms = term_edf(fit),
    sp = fit_sp(fit)
  )
}

truth_meta <- function(app, source, d, fit) {
  c(
    list(schema = PLASMODE_CACHE_SCHEMA, app = app, source = source),
    frame_meta(d, app),
    list(sp = fit_sp(fit))
  )
}

# Everything a cache must agree with the current code and data on (sp is the fit's).
frame_meta <- function(d, app) {
  list(
    refund_sha = substr(loaded_refund_sha(), 1, 8),
    frame = digest::digest(d$ids),
    data = digest::digest(list(d$Y, d$X, d$z)), # content, not only the frame
    G = length(d$ids),
    yind = d$yind,
    eval_t = d$eval_t,
    s_grid = S_GRID,
    test_rows = d$test_rows,
    model = digest::digest(list(
      "Y ~ ff(X) + z",
      BASES[[plasmode_fit_cell(app)$basis]]
    )),
    blocks = digest::digest(list(d$block, d$block_alt))
  )
}

check_truth_meta <- function(truth, app, source) {
  d <- plasmode_data(app)
  expected <- c(
    list(schema = PLASMODE_CACHE_SCHEMA, app = app, source = source),
    frame_meta(d, app)
  )
  stored <- truth$meta[names(expected)]
  stale <- names(expected)[!mapply(identical, stored, expected)]
  if (length(stale)) {
    stop(
      "stale truth cache ",
      truth_file(app, source),
      " (differs in: ",
      toString(stale),
      "); delete it and the app's results"
    )
  }
  if (source == "MID") {
    sp <- mid_sp(
      truth_fit_info(app, "REML")$sp,
      truth_fit_info(app, "NCV")$sp
    )
    if (!isTRUE(all.equal(sp, truth$meta$sp))) {
      stop("MID cache does not match the current REML/NCV truths of ", app)
    }
  }
  invisible(TRUE)
}

# Per-smooth EDF with the smooth's null-space and full (penalized-range) dimensions.
term_edf <- function(fit) {
  data.frame(
    term = vapply(fit$smooth, \(sm) sm$label, ""),
    edf = vapply(
      fit$smooth,
      \(sm) sum(fit$edf[sm$first.para:sm$last.para]),
      numeric(1)
    ),
    null = vapply(fit$smooth, \(sm) as.numeric(sm$null.space.dim), numeric(1)),
    full = vapply(
      fit$smooth,
      \(sm) as.numeric(sm$last.para - sm$first.para + 1),
      numeric(1)
    )
  )
}

# Smoothing parameters and per-term EDF of an app's REML or NCV truth fit. Caches of the
# original apps have neither: their fit is repeated on the cached data and must reproduce
# the cached means.
truth_fit_info <- function(app, source) {
  key <- paste("fit-info", app, source)
  if (!is.null(.plasmode_memo[[key]])) return(.plasmode_memo[[key]])
  truth <- plasmode_truth(app, source)
  info <- if (!is.null(truth$meta)) {
    list(sp = truth$meta$sp, terms = truth$terms)
  } else {
    sim <- plasmode_sim(truth$data, truth$data$Y)
    fit <- fit_model(sim, plasmode_fit_cell(app), source)
    refit <- truth_from_fit(fit, truth$data, sim)
    gap <- max(abs(refit$mu - truth$mu)) / stats::sd(truth$mu)
    if (gap > 1e-5) {
      # numerical noise across machines is ~1e-7
      stop("refit of the cached ", app, "/", source, " truth differs: ", gap)
    }
    list(sp = fit_sp(fit), terms = refit$terms)
  }
  assign(key, info, envir = .plasmode_memo)
  info
}

# A fit's smoothing parameters; mgcv reports fixed ones in full.sp only.
fit_sp <- function(fit) if (length(fit$sp)) fit$sp else fit$full.sp

# Elementwise log-midpoint of two named sp vectors (aligned by name).
mid_sp <- function(sp_reml, sp_ncv) {
  if (
    !setequal(names(sp_reml), names(sp_ncv)) || anyDuplicated(names(sp_reml))
  ) {
    stop("REML and NCV smoothing parameters do not align by name")
  }
  sp_ncv <- sp_ncv[names(sp_reml)]
  exp((log(sp_reml) + log(sp_ncv)) / 2)
}

# The MID fit (all penalties fixed at the log-midpoint) and its interior check: every
# term's EDF strictly between the REML and NCV fits', and neither of those at a term's
# null-space or full EDF (a selector at an optimiser boundary).
mid_fit <- function(app, sim) {
  sp_reml <- truth_fit_info(app, "REML")$sp
  sp_ncv <- truth_fit_info(app, "NCV")$sp
  sp <- mid_sp(sp_reml, sp_ncv)
  fit <- collect_warnings(fit_pffr(
    sim,
    plasmode_fit_cell(app),
    "REML",
    sp = sp
  ))$value
  if (!isTRUE(all.equal(unname(fit_sp(fit)), unname(sp)))) {
    stop("pffr did not keep the fixed MID smoothing parameters")
  }
  terms <- term_edf(fit)
  edf_of <- \(source) {
    other <- truth_fit_info(app, source)$terms
    other$edf[match(terms$term, other$term)]
  }
  terms$edf_reml <- edf_of("REML")
  terms$edf_ncv <- edf_of("NCV")
  list(
    fit = fit,
    sp_reml = sp_reml,
    sp_ncv = sp_ncv,
    terms = terms,
    interior = mid_interior(terms)
  )
}

mid_interior <- function(terms) {
  lo <- pmin(terms$edf_reml, terms$edf_ncv)
  hi <- pmax(terms$edf_reml, terms$edf_ncv)
  share <- \(edf) (edf - terms$null) / (terms$full - terms$null)
  edge <- \(edf) share(edf) < MID_EDF_EDGE | share(edf) > 1 - MID_EDF_EDGE
  between <- terms$edf > lo & terms$edf < hi
  isTRUE(all(between & !edge(terms$edf_reml) & !edge(terms$edf_ncv)))
}

# Errors -------------------------------------------------------------------------------------

# Error curves of a residual source on the truth's subject frame: the residuals of the
# named fit of the real data, or (ECG) the beat-to-beat differences, scaled by c (beat_m:
# average variance of the REML residuals) or 1 (beat_1).
plasmode_errors <- function(app, residual, truth) {
  if (residual %in% c("REML", "NCV", "MID")) {
    source <- plasmode_truth(app, residual)
    # the original apps' caches hold no ids: then Y, X and z must agree
    ids_differ <- !is.null(source$data$ids) &&
      !is.null(truth$data$ids) &&
      !identical(source$data$ids, truth$data$ids)
    if (
      !identical(source$data$Y, truth$data$Y) ||
        !identical(source$data$X, truth$data$X) ||
        !identical(source$data$z, truth$data$z) ||
        ids_differ
    ) {
      stop("residual source ", residual, " is on another subject frame")
    }
    return(source$residuals)
  }
  beats <- plasmode_beats(app)
  if (
    !identical(beats$ids, truth$data$ids) ||
      !identical(beats$yind, truth$data$yind)
  ) {
    stop("beat frame differs from the truth frame (ids, order or grid)")
  }
  switch(
    residual,
    beat_m = beats$c * beats$diff,
    beat_1 = beats$diff,
    stop("unknown residual source: ", residual)
  )
}

# Beat differences of an app (ECG only) with their scale c, cached with the truths.
plasmode_beats <- function(app) {
  if (app != "ecg_gls4b") stop("beat-to-beat errors exist for ecg_gls4b only")
  file <- file.path(PLASMODE_CACHE, paste0("beats-", app, ".rds"))
  if (!is.null(.plasmode_memo[[file]])) return(.plasmode_memo[[file]])
  reml <- plasmode_truth(app, "REML")
  if (file.exists(file)) {
    beats <- readRDS(file)
    if (
      !identical(beats$reml_frame, reml$meta$frame) ||
        !identical(beats$reml_data, reml$meta$data) ||
        !identical(beats$reml_sp, reml$meta$sp)
    ) {
      stop("stale beat cache ", file, ": delete it and the app's results")
    }
  } else {
    beats <- ecg_gls4b_beats(reml$data$ids)
    beats$yind <- reml$data$yind
    beats$c <- beat_scale(reml$residuals, beats$diff)
    beats$reml_frame <- reml$meta$frame
    beats$reml_data <- reml$meta$data
    beats$reml_sp <- reml$meta$sp
    atomic_save(beats, file)
  }
  assign(file, beats, envir = .plasmode_memo)
  beats
}

# c such that c * beat differences have the REML residuals' average variance over t; the
# variance of sign-flipped errors at t is the mean of squares over subjects.
beat_scale <- function(residuals, diff) sqrt(mean(residuals^2) / mean(diff^2))

# Replicates ---------------------------------------------------------------------------------

# A study-format sim object from plasmode data rows and responses Y; `blocks` names the
# block-id columns (block, block_alt) to attach to the data for the block arm.
plasmode_sim <- function(
  d,
  Y,
  rows = seq_len(nrow(Y)),
  targets = NULL,
  blocks = character()
) {
  test <- d$test_rows %||% seq_len(N_TEST) # the test cohort's covariates
  data <- data.frame(X = I(d$X[rows, ]), z = d$z[rows], Y = I(Y))
  for (b in blocks) {
    if (is.null(d[[b]])) stop("no ", b, " ids in this dataset")
    data[[b]] <- d[[b]][rows]
  }
  list(
    data = data,
    yind = d$yind,
    xind = d$xind,
    eval_t = d$eval_t %||% EVAL_T,
    test = data.frame(X = I(d$X[test, ]), z = d$z[test]),
    truth = targets,
    diag = list()
  )
}

# All random draws of replicate `rep` for n subjects, in a fixed order: the first three
# are the original ones, the block signs are appended (one per possible block id).
plasmode_draws <- function(n, rep) {
  withr::with_seed(
    SEED_BASE + rep,
    list(
      subset = sample.int(n),
      sign = sample(c(-1, 1), n, replace = TRUE),
      donor = sample.int(n),
      block_sign = sample(c(-1, 1), n, replace = TRUE)
    )
  )
}

# Replicate `rep` of a plasmode cell. The draws are made for all subjects in a fixed order,
# so cells of one app share them (common random numbers across truths, residual sources
# and G).
plasmode_simulate <- function(cell, rep) {
  truth <- plasmode_truth(cell$app, cell$truth)
  errors <- plasmode_errors(cell$app, cell$residual %||% cell$truth, truth)
  n <- nrow(truth$mu)
  draws <- plasmode_draws(n, rep)
  rows <- if (cell$G < n) sort(draws$subset[seq_len(cell$G)]) else seq_len(n)
  donor <- if (cell$error == "attached") rows else draws$donor[rows]
  extras <- strsplit(trimws(cell$extra), " +")[[1]]
  blocks <- intersect(c("block", "block_alt"), extras)
  if (length(blocks) && cell$G < n) stop("block ids exist at G = all only")
  signs <- switch(
    cell$flip %||% "curve",
    curve = draws$sign[rows],
    block = {
      if (cell$G < n || is.null(truth$data$block)) {
        stop("block flips need all subjects and block ids")
      }
      draws$block_sign[truth$data$block[rows]]
    }
  )
  Y <- truth$mu[rows, ] + signs * errors[donor, ]
  plasmode_sim(truth$data, Y, rows, truth$targets, blocks)
}
