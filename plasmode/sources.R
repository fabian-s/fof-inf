# Plasmode data sources. Each source returns one curve per subject on its fitted grids:
# list(Y (G x D, on `yind`), X (G x length(xind)), z, yind, xind, eval_t, ids) plus, for
# the datasets with dependent curves, `block` and `block_alt` (one id per curve).
# PRIVACY: the cardiology data (ecg_gls4, ecg_gls4b, af) are patient data. They are read
# from their local location only; never copy them, anything derived from them per subject
# (fits, residual or beat curves, task files) or these paths' contents into the
# repository or to other machines.
# The public datasets are read from plain .rds files written once by
# plasmode/export-public.R, so a cluster needs no package beyond the study library.

# Data locations are set by environment variables; the defaults lie beside the repository.
PLASMODE_CACHE <- path.expand(
  Sys.getenv("PLASMODE_CACHE", "../plasmode-cache") # outside the repository
)
PLASMODE_PUBLIC_DATA <- path.expand(
  Sys.getenv("PLASMODE_PUBLIC_DATA", "../plasmode-public-data")
)

# Cardiology data (not public): the ECG study's project directory and the AF trial's file.
ECG_STUDY <- Sys.getenv("PLASMODE_ECG_DIR")
AF_FILE <- Sys.getenv("PLASMODE_AF_FILE")

# Apps whose data must never leave this machine (checked by run-plasmode.R).
PLASMODE_LOCAL_ONLY <- c("ecg_gls4", "ecg_gls4b", "af")
# Apps of the original 16 cells: fixed test cohort (first N_TEST) and legacy caches.
# (running is not one: its cells use the checked cohort)
PLASMODE_LEGACY_APPS <- "ecg_gls4"

# Cardiology input: "original" (the local files above; this machine only) or "deid" (the
# de-identified numeric arrays under PLASMODE_DEID_DATA/<app>/data.rds: no ids, names,
# dates or other covariates; the only input allowed on a cluster, and only for the apps
# in PLASMODE_DEID_APPS).
PLASMODE_CARDIO <- Sys.getenv("PLASMODE_CARDIO", "original")
PLASMODE_DEID_DATA <- path.expand(
  Sys.getenv("PLASMODE_DEID_DATA", file.path(PLASMODE_CACHE, "deid"))
)
PLASMODE_DEID_APPS <- c("ecg_gls4b", "af")

# Test cohort: the first N_TEST subjects, unless they are one-sided in z or in the
# covariate's leading scores given the row order; then a fixed
# seeded sample of N_TEST subjects.
TEST_COHORT_P <- 0.05 # permutation p-value of max |SMD| below which the cohort is one-sided
TEST_COHORT_N_PERM <- 2000
TEST_COHORT_SCORES <- 3 # leading covariate FPC scores checked besides z

standardize <- function(x) as.numeric(scale(x))

# tf's tfd() looks its evaluator up by name, so tf must be attached where tf curves are
# built or evaluated (the cardiology loaders; export-public.R).
attach_tf <- function() suppressPackageStartupMessages(library(tf))

# Curves (tf vector, regular or irregular, any domain) mapped to [0, 1] and linearly
# interpolated on `grid`: a matrix with one row per curve.
on_grid <- function(f, grid) {
  attach_tf()
  f <- tf::tfd(f)
  domain <- tf::tf_domain(f)
  args <- tf::tf_arg(f)
  values <- tf::tf_evaluations(f)
  if (!is.list(args)) args <- rep(list(args), length(values))
  rescaled <- \(a) (a - domain[1]) / diff(domain)
  t(mapply(
    \(a, v) stats::approx(rescaled(a), v, xout = grid, rule = 2)$y,
    args,
    values
  ))
}

# Mean-estimand evaluation times on a native grid: all grid points
# for D <= 73, every second one for D = 101; EVAL_T on the study's D = 61 grid.
native_eval_t <- function(yind) {
  D <- length(yind)
  if (D == 61 && isTRUE(all.equal(yind, t_grid(61)))) return(EVAL_T)
  if (D <= 73) return(yind)
  if (D == 101) return(yind[seq(1, D, by = 2)])
  stop("no mean-evaluation rule for D = ", D)
}

# Cardiology loaders (local only) --------------------------------------------------------

# GLS4 beats of the ECG study with age: one row per beat, beats in file order per subject.
ecg_gls4_beats <- function() {
  base <- path.expand(ECG_STUDY)
  beats <- readRDS(file.path(
    base,
    "Data and Code for FOFR This time it should be correct!/data/data_f.rds"
  ))
  scalars <- readRDS(file.path(
    base,
    "Strain curves and ECG/data/scalar.rds"
  )) |>
    dplyr::ungroup() |>
    dplyr::mutate(id = as.character(id))
  beats |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::filter(strain_curve == "GLS4") |>
    dplyr::group_by(id) |>
    dplyr::mutate(beat = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::inner_join(scalars, by = "id") |>
    dplyr::filter(!is.na(age))
}

# ECG -> GLS4 strain, first beat per subject; z = age. `min_beats` = 2 gives the
# 78-subject frame of the extension (beat-difference errors need two beats).
ecg_gls4_source <- function(min_beats = 1L) {
  b <- ecg_gls4_beats()
  n_beats <- table(b$id)
  keep <- names(n_beats)[n_beats >= min_beats]
  first <- b[b$beat == 1 & b$id %in% keep, ]
  # sorted by id in the C locale, as the original loader's group_by() + slice(1)
  first <- first[order(first$id, method = "radix"), ]
  list(
    Y = on_grid(first$strain_val, t_grid(61)),
    X = on_grid(first$ecg_val, S_GRID),
    z = standardize(first$age),
    ids = first$id
  )
}

# Beat-to-beat differences (Y_beat1 - Y_beat2) / sqrt(2) on t_grid(61), for `ids` in order.
ecg_gls4b_beats <- function(ids) {
  if (identical(PLASMODE_CARDIO, "deid")) {
    d <- read_public("ecg_gls4b", PLASMODE_DEID_DATA)
    if (!identical(d$ids, ids)) stop("beat frame differs from the truth frame")
    return(list(diff = d$beat_diff, ids = ids))
  }
  b <- ecg_gls4_beats()
  beat_curves <- \(k) {
    rows <- b[b$beat == k, ]
    rows <- rows[match(ids, rows$id), ]
    if (anyNA(rows$id)) stop("subjects without beat ", k)
    if (!identical(rows$id, ids))
      stop("beat frame differs from the truth frame")
    on_grid(rows$strain_val, t_grid(61))
  }
  list(diff = (beat_curves(1) - beat_curves(2)) / sqrt(2), ids = ids)
}

# AF trial: LA strain (4-chamber) at day 30 ~ day 0 + randomised treatment (digoxin = 1).
af_source <- function() {
  attach_tf()
  a <- readRDS(path.expand(AF_FILE))
  grid <- tf::tf_arg(a$ch4d30)
  yind <- (grid - min(grid)) / diff(range(grid))
  list(
    Y = on_grid(a$ch4d30, yind),
    X = on_grid(a$ch4d0, S_GRID),
    z = as.numeric(a$case),
    yind = yind,
    ids = as.character(a$id)
  )
}

# The cardiology apps of the extension, from the original or the de-identified input.
cardio_source <- function(app) {
  if (identical(PLASMODE_CARDIO, "deid")) {
    return(read_public(app, PLASMODE_DEID_DATA))
  }
  if (!identical(PLASMODE_CARDIO, "original")) {
    stop("PLASMODE_CARDIO must be 'original' or 'deid'")
  }
  switch(app, ecg_gls4b = ecg_gls4_source(min_beats = 2L), af = af_source())
}

# Public sets: the exported .rds files ------------------------------------------------------

# (also reads the de-identified cardiology files, which carry beat_diff for ECG)
read_public <- function(set, dir = PLASMODE_PUBLIC_DATA) {
  file <- file.path(dir, set, "data.rds")
  if (!file.exists(file)) {
    stop("data file missing: ", file, " (run plasmode/export-public.R)")
  }
  d <- readRDS(file)
  hash <- digest::digest(d[intersect(
    c(PUBLIC_HASHED_FIELDS, "beat_diff"),
    names(d)
  )])
  if (!identical(hash, d$hash)) stop("content hash mismatch in ", file)
  d
}
PUBLIC_HASHED_FIELDS <- c(
  "Y",
  "X",
  "z",
  "yind",
  "xind",
  "eval_t",
  "ids",
  "block",
  "block_alt"
)

PLASMODE_SOURCES <- list(
  ecg_gls4 = function() ecg_gls4_source(min_beats = 1L),
  ecg_gls4b = function() cardio_source("ecg_gls4b"),
  af = function() cardio_source("af"),
  running = function() read_public("running"),
  dti = function() read_public("dti"),
  gait_vc = function() read_public("gait_vc"),
  ecg_roahd = function() read_public("ecg_roahd"),
  ocean = function() read_public("ocean"),
  aemet = function() read_public("aemet"),
  adelaide = function() read_public("adelaide")
)

# One source, checked on ingestion: incomplete curves dropped (legacy behaviour), Y on
# its own strictly increasing yind, eval_t a subset of yind, X on xind, one z and id per
# curve; test_rows from the one-sidedness check (legacy apps: the first N_TEST).
plasmode_data <- function(app) {
  raw <- PLASMODE_SOURCES[[app]]()
  if (is.null(raw)) stop("unknown plasmode app: ", app)
  yind <- raw$yind %||% t_grid(61)
  xind <- raw$xind %||% S_GRID
  ok <- stats::complete.cases(raw$Y) &
    stats::complete.cases(raw$X) &
    !is.na(raw$z)
  d <- list(
    Y = raw$Y[ok, , drop = FALSE],
    X = raw$X[ok, , drop = FALSE],
    z = raw$z[ok],
    yind = yind,
    xind = xind,
    eval_t = raw$eval_t %||% native_eval_t(yind),
    ids = as.character(raw$ids %||% seq_along(raw$z))[ok]
  )
  if (!is.null(raw$block)) d$block <- raw$block[ok]
  if (!is.null(raw$block_alt)) d$block_alt <- raw$block_alt[ok]
  check_plasmode_data(d)
  d$test_check <- test_cohort_check(d$X, d$z)
  d$test_rows <- if (app %in% PLASMODE_LEGACY_APPS || !d$test_check$one_sided) {
    seq_len(N_TEST)
  } else {
    withr::with_seed(SEED_TEST, sort(sample.int(length(d$z), N_TEST)))
  }
  d
}

check_plasmode_data <- function(d) {
  G <- length(d$z)
  stopifnot(
    "yind must be strictly increasing" = all(diff(d$yind) > 0),
    "Y must have one column per yind point" = ncol(d$Y) == length(d$yind),
    "X must have one column per xind point" = ncol(d$X) == length(d$xind),
    "one Y row and id per z" = nrow(d$Y) == G && length(d$ids) == G,
    "ids must be unique" = !anyDuplicated(d$ids),
    "eval_t must be a subset of yind" = all(
      vapply(d$eval_t, \(t) any(abs(d$yind - t) < 1e-9), logical(1))
    ),
    "G must be at least N_TEST" = G >= N_TEST
  )
  for (b in c("block", "block_alt")) {
    if (!is.null(d[[b]]) && (length(d[[b]]) != G || anyNA(d[[b]]))) {
      stop(b, " must hold one id per curve")
    }
  }
  invisible(d)
}

# Is the first-N_TEST cohort one-sided? Standardised mean differences (first N_TEST vs the
# rest) of z and of the covariate's leading centred FPC scores; the maximum |SMD| is
# compared with its distribution over random cohorts of the same size.
test_cohort_check <- function(X, z) {
  scores <- stats::prcomp(X)$x[, seq_len(TEST_COHORT_SCORES), drop = FALSE]
  vars <- cbind(z = z, scores)
  smd <- \(rows) {
    apply(vars, 2, \(v) {
      (mean(v[rows]) - mean(v[-rows])) / stats::sd(v)
    })
  }
  observed <- smd(seq_len(N_TEST))
  null <- withr::with_seed(SEED_TEST + 1, {
    replicate(TEST_COHORT_N_PERM, max(abs(smd(sample.int(length(z), N_TEST)))))
  })
  p <- mean(null >= max(abs(observed)))
  list(smd = observed, p = p, one_sided = p < TEST_COHORT_P)
}
