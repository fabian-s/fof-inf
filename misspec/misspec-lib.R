# Misspecified-truth stress test: shared definitions. Sourced after scripts/setup.R
# (the study's constants, data generation and fit code), from the repository root.
#
# Truths T1 ("bump") and T2 ("weak") are NOT projected onto the fitted bases: beta is the
# study's (projected, scaled) smooth beta plus a feature outside the default spline space,
# rescaled to the study's signal (main-term variance SIGNAL_SHARE_MAIN of S = 1). T0 is the
# study's own truth ("smooth"). alpha and gamma are the study's.

# Constants -------------------------------------------------------------------------

FEATURE_SHARE <- 0.2 # share of the beta term's signal variance carried by the feature
MS_BUMP_CENTER <- c(s = 0.3, t = 0.6)
MS_BUMP_SD <- 0.05 # in s and in t
MS_BUMP_REGION_SD <- 2 # beta points within 2 sd (elliptical) of the centre = bump region
WEAK_J <- 12:16 # cos(j pi s), covariate variance (j + 1)^-2
WEAK_T <- \(t) sin(pi * t) # smooth function of t multiplying the weak component
MS_TRUTHS <- c(T0 = "smooth", T1 = "bump", T2 = "weak")
MS_FIT_GRID <- 101 # grid for relative L2 projection errors
N_CHECK_REUSED <- 4 # reps recomputed in cells whose records are reused from the study
STUDY_RESULTS <- "../results" # task-cCCC-rRRR.rds of cells 96-123 (LRZ, from code copy)
STUDY_CORE_RESULTS <- "../results/373f6fd65cac" # task-<id>.rds of the 95-cell run
STUDY_N_CELLS <- 95

# Feature surfaces (raw, unscaled): functions (s, t) -> length(s) x length(t) matrix.
MS_FEATURES <- list(
  bump = \(s, t) {
    outer(
      exp(-(s - MS_BUMP_CENTER[["s"]])^2 / (2 * MS_BUMP_SD^2)),
      exp(-(t - MS_BUMP_CENTER[["t"]])^2 / (2 * MS_BUMP_SD^2))
    )
  },
  weak = \(s, t) {
    outer(rowSums(sqrt(2) * cos(outer(s, WEAK_J) * pi)), WEAK_T(t))
  }
)

# Truths ----------------------------------------------------------------------------

# The study's truth constructor, kept under another name: make_truth() is masked below so
# that the study's simulate() builds the misspecified truths through its own code path
# (T0 then runs exactly the study's code).
if (!exists("study_make_truth")) study_make_truth <- make_truth

make_truth <- function(model, roughness, x_sd = X_SD) {
  if (!roughness %in% names(MS_FEATURES)) {
    return(study_make_truth(model, roughness, x_sd))
  }
  if (model != "M1") stop("misspecified truths are defined for M1 only")
  misspec_truth(roughness, x_sd)$truth
}

# beta = a (beta0 + c feature): c gives the feature FEATURE_SHARE of the beta term's
# variance V(f) = int Var_X(int X f ds) dt, and a restores V = SIGNAL_SHARE_MAIN.
misspec_truth <- function(feature, x_sd = X_SD) {
  base <- study_make_truth("M1", "smooth", x_sd)
  f <- MS_FEATURES[[feature]]
  V <- \(fun) main_term_variance(fun, "M1", x_sd)
  share <- \(c)
    V(\(s, t) c * f(s, t)) / V(\(s, t) base$beta(s, t) + c * f(s, t))
  c_feat <- uniroot(
    \(c) share(c) - FEATURE_SHARE,
    c(1e-6, 1e3),
    tol = 1e-12
  )$root
  raw <- \(s, t) base$beta(s, t) + c_feat * f(s, t)
  a <- sqrt(SIGNAL_SHARE_MAIN / V(raw))
  truth <- base
  truth$beta <- \(s, t) a * raw(s, t)
  list(
    truth = truth,
    smooth_part = \(s, t) a * base$beta(s, t),
    feature_part = \(s, t) a * c_feat * f(s, t),
    amplitude = a * c_feat,
    feature_share = share(c_feat)
  )
}

# Relative L2 error (on a fine [0, 1]^2 grid) of the least-squares projection of g onto a
# basis's ff tensor space, and the share of the beta term's variance (int Var int X r) that
# the projection residual r carries.
projection_error <- function(g, basis_name, x_sd = X_SD) {
  ff <- BASES[[basis_name]]$ff
  margin <- \(i) ps_basis(list(bs = ff$bs, k = ff$k[i], m = ff$m[[i]]))
  bu <- margin(1)
  bv <- margin(2)
  u <- seq(0, 1, length.out = MS_FIT_GRID)
  # as project_tensor(), for g returning the length(s) x length(t) matrix itself
  Bu <- ps_eval(bu, u)
  Bv <- ps_eval(bv, u)
  theta <- solve(crossprod(Bu), t(Bu)) %*%
    g(u, u) %*%
    t(solve(crossprod(Bv), t(Bv)))
  proj <- \(s, t) ps_eval(bu, s) %*% theta %*% t(ps_eval(bv, t))
  resid <- \(s, t) g(s, t) - proj(s, t)
  c(
    rel_l2 = sqrt(sum(resid(u, u)^2) / sum(g(u, u)^2)),
    resid_signal_share = main_term_variance(resid, "M1", x_sd) /
      main_term_variance(g, "M1", x_sd)
  )
}

# Description of the truths: amplitudes, signal shares and projection errors.
describe_truths <- function() {
  rows <- lapply(names(MS_TRUTHS), \(label) {
    name <- MS_TRUTHS[[label]]
    beta <- make_truth("M1", name)$beta
    info <- if (name %in% names(MS_FEATURES)) misspec_truth(name) else NULL
    u <- seq(0, 1, length.out = MS_FIT_GRID)
    l2_share <- if (is.null(info)) 0 else {
      sum(info$feature_part(u, u)^2) / sum(beta(u, u)^2)
    }
    errs <- lapply(c("default", "xlarge"), \(b) projection_error(beta, b))
    tibble::tibble(
      truth = label,
      name = name,
      feature_amplitude = info$amplitude %||% 0,
      feature_signal_share = info$feature_share %||% 0,
      feature_l2_share = l2_share,
      beta_signal_var = main_term_variance(beta, "M1", X_SD),
      rel_l2_default = errs[[1]][["rel_l2"]],
      resid_signal_share_default = errs[[1]][["resid_signal_share"]],
      rel_l2_xlarge = errs[[2]][["rel_l2"]],
      resid_signal_share_xlarge = errs[[2]][["resid_signal_share"]]
    )
  })
  dplyr::bind_rows(rows)
}

# Cells and tasks ---------------------------------------------------------------------

# Gaussian, default signal, D = 61, M1: truth x {iid G100, fpc G100, fpc G40} x basis.
# study_cell: the study cell with identical data and fits (T0 only), else NA.
misspec_cells <- function() {
  combos <- tibble::tribble(
    ~error,
    ~G,
    "iid",
    100L,
    "fpc",
    100L,
    "fpc",
    40L
  )
  cells <- tidyr::expand_grid(
    truth_label = names(MS_TRUTHS),
    combos,
    basis = c("default", "xlarge")
  )
  defaults <- cell_defaults()
  for (col in setdiff(names(defaults), c("error", "G", "basis"))) {
    cells[[col]] <- defaults[[col]]
  }
  cells$truth <- unname(MS_TRUTHS[cells$truth_label])
  cells$fit_family <- cells$family
  cells$block <- "misspec"
  study <- design_table()
  key <- \(d) paste(d$error, d$G, d$basis, d$truth, d$signal, d$family, d$model)
  hit <- study[
    study$D == 61 &
      study$covariate == "rich" &
      study$block %in% c("core", "basis_size"),
  ]
  cells$study_cell <- hit$cell[match(key(cells), key(hit))]
  cells$cell <- seq_len(nrow(cells))
  stopifnot(sum(!is.na(cells$study_cell)) == 5)
  dplyr::relocate(cells, cell)
}

# Cells whose full study records are reused (expensive xlarge T0 cells): only the first
# N_CHECK_REUSED replicates are recomputed there, as a common-random-numbers check.
reused_cell <- function(cells)
  !is.na(cells$study_cell) & cells$basis == "xlarge"

# Rep-major task table; tasks are looked up by id. `ncv_reps`: replicates in which the
# xlarge cells also get the NCV fit (all default-basis tasks do).
misspec_tasks <- function(cells = misspec_cells(), reps = seq_len(N_REP)) {
  tasks <- tidyr::expand_grid(rep = reps, cell = cells$cell)
  reused <- cells$cell[reused_cell(cells)]
  tasks <- tasks[!(tasks$cell %in% reused & tasks$rep > N_CHECK_REUSED), ]
  tasks$task <- seq_len(nrow(tasks))
  dplyr::relocate(tasks, task)
}

# Path of the study's record of (study cell, rep).
study_record_path <- function(study_cell, rep) {
  if (study_cell <= STUDY_N_CELLS) {
    return(file.path(
      STUDY_CORE_RESULTS,
      sprintf("task-%d.rds", (rep - 1) * STUDY_N_CELLS + study_cell)
    ))
  }
  file.path(STUDY_RESULTS, sprintf("task-c%03d-r%03d.rds", study_cell, rep))
}
