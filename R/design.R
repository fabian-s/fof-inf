# The one design table: a core factorial plus extensions, one row per cell.
# Every column not varied by a block keeps its default. A block is only rows.

cell_defaults <- function() {
  list(
    family = "gaussian", # generating family
    fit_family = NA_character_, # fitted family; NA = same as `family`
    error = "fpc",
    G = 100L,
    signal = "mid",
    truth = "smooth",
    D = 61L,
    model = "M1",
    covariate = "rich", # spectrum of the functional covariate: rich or lowrank
    basis = "default", # fitted bases, a name in BASES (R/constants.R)
    extra = "" # extra fits: "pointwise", "ar1"
  )
}

# Build one block from the factors it varies; all other columns get defaults.
design_block <- function(block, ...) {
  varied <- tidyr::expand_grid(...)
  defaults <- cell_defaults()
  fixed <- defaults[setdiff(names(defaults), names(varied))]
  dplyr::bind_cols(
    tibble::tibble(block = block),
    varied,
    tibble::as_tibble(fixed)
  )
}

# Signal levels of a family, in the order of SIGNAL_LEVELS.
signal_levels <- function(family) names(SIGNAL_LEVELS[[family]])

# A block crossed separately per family, because the signal levels differ by family.
per_family_block <- function(block, families, ...) {
  purrr::map(families, \(fam) {
    design_block(block, family = fam, signal = signal_levels(fam), ...)
  }) |>
    dplyr::bind_rows()
}

design_table <- function() {
  main_families <- c("gaussian", "poisson", "binomial")
  dependent <- c("iid", "ou", "fpc")
  cells <- dplyr::bind_rows(
    per_family_block(
      "core",
      main_families,
      error = dependent,
      G = c(40L, 100L)
    ),
    design_block(
      "warp",
      family = main_families,
      error = "warp",
      G = c(40L, 100L),
      truth = c("smooth", "wiggly"),
      signal = "high"
    ),
    design_block(
      "dense_grid",
      family = main_families,
      error = dependent,
      D = 241L
    ),
    per_family_block(
      "rough_truth",
      main_families,
      error = c("iid", "fpc"),
      truth = "wiggly"
    ),
    design_block(
      "term_type",
      family = main_families,
      error = dependent,
      model = "M2"
    ),
    design_block(
      "families",
      family = c("scat", "betar", "nb"),
      error = c("iid", "fpc")
    ),
    design_block(
      "families",
      family = "nb",
      fit_family = "poisson",
      error = c("iid", "fpc")
    ),
    design_block("ar1_home", error = "ou0", extra = "ar1"),
    # Extension (2026-09-25), appended so that the ids of cells 1-95 are unchanged:
    # an error process whose correlations change sign, and AR(1) on the wiggly warp.
    design_block("oscillating", family = main_families, error = "osc"),
    design_block(
      "warp_ar1",
      error = "warp",
      truth = "wiggly",
      signal = "high",
      extra = "ar1"
    ),
    # Extension (2026-09-26): heteroskedastic Gaussian errors, along t (real DTI variance
    # profile) and across subjects (scale depending on z), appended after cell 99.
    design_block(
      "heteroskedastic",
      error = c("fpc_raw", "fpc_het", "fpc_raw_het"),
      G = c(40L, 100L)
    ),
    # Extension (2026-09-26): a low-rank covariate like the running data's (rank 9 < 1.5
    # k_s), to test whether NCV's advantage depends on how well beta is identified.
    design_block(
      "lowrank_covariate",
      family = main_families,
      error = c("iid", "fpc"),
      covariate = "lowrank"
    ),
    # Extension (2026-09-28): larger nested bases (same truth, more over-parameterized
    # fit), compared with the default-basis core cells at G = 100 by common random numbers.
    design_block(
      "basis_size",
      family = main_families,
      error = c("iid", "fpc"),
      basis = c("large", "xlarge")
    ),
    # Extension (2026-10-03): the core factorial with 20 curves (first 20 of each
    # replicate's draws), and a null effect: beta = 0 in the smooth-truth setting
    # (gamma, intercept and Gaussian error SD unchanged), Gaussian and binary.
    per_family_block("small_G", main_families, error = dependent, G = 20L),
    design_block(
      "null_beta",
      family = c("gaussian", "binomial"),
      error = c("iid", "fpc"),
      G = c(40L, 100L),
      truth = "null"
    )
  )
  cells <- mark_extra_fits(cells)
  osc_gaussian <- cells$block == "oscillating" & cells$family == "gaussian"
  cells$extra[osc_gaussian] <- "pointwise ar1"
  het_ar1 <- cells$block == "heteroskedastic" & cells$G == 100L
  cells$extra[het_ar1] <- "ar1"
  cells$fit_family <- dplyr::coalesce(cells$fit_family, cells$family)
  cells$cell <- seq_len(nrow(cells))
  dplyr::relocate(cells, cell)
}

# Extra fits run inside existing cells (Gaussian, G = 100; core at mid signal, warp with
# the smooth truth):
# pointwise NCV under iid and fpc; AR(1) under ou, fpc and warp.
mark_extra_fits <- function(cells) {
  host <- cells$family == "gaussian" &
    cells$G == 100L &
    ((cells$block == "core" & cells$signal == "mid") |
      (cells$block == "warp" & cells$truth == "smooth"))
  add <- function(extra, errors) {
    ifelse(host & cells$error %in% errors, extra, "")
  }
  extras <- paste(
    add("pointwise", c("iid", "fpc")),
    add("ar1", c("ou", "fpc", "warp"))
  )
  cells$extra <- ifelse(cells$extra == "", trimws(extras), cells$extra)
  cells
}

# Tasks: every cell x replicate, enumerated rep-major so a partial run leaves
# equal replicate counts per cell. Tasks are looked up by id, never by position.
task_table <- function(cells = design_table(), reps = seq_len(N_REP)) {
  tasks <- tidyr::expand_grid(rep = reps, cell = cells$cell)
  tasks$task <- seq_len(nrow(tasks))
  dplyr::relocate(tasks, task)
}
