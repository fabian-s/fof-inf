# Plasmode main design: design, truths, draws, native grids, block arm, validation keys
# and the decision rules.
source(testthat::test_path("..", "analysis", "contrasts.R"))
source(testthat::test_path("..", "plasmode", "sources.R"))
source(testthat::test_path("..", "plasmode", "plasmode.R"))
source(testthat::test_path("..", "plasmode", "rule.R"))

# A stand-in public dataset: G = 56 curves on a native D = 48 grid, blocks of two.
standin_source <- function(G = 56L, D = 48L) {
  withr::with_seed(7, {
    t <- seq(0, 1, length.out = D)
    X <- matrix(rnorm(G * 6), G) %*% (cos(outer(0:5, S_GRID * pi)) / (1:6))
    z <- rnorm(G)
    signal <- outer(rowMeans(X), sin(pi * t)) + outer(z, 0.5 * cos(pi * t))
    errors <- matrix(rnorm(G * D), G) %*% chol(0.7^abs(outer(1:D, 1:D, "-")))
    list(
      Y = 1 + signal + 0.5 * errors,
      X = X,
      z = z,
      yind = t,
      ids = sprintf("s%02d", seq_len(G)),
      block = as.integer(ceiling(seq_len(G) / 2))
    )
  })
}
STANDIN <- list(standin = list(G = 56L, D = 48L, block = TRUE))

local_standin <- function(env = parent.frame()) {
  cache <- withr::local_tempdir(.local_envir = env)
  old <- list(cache = PLASMODE_CACHE, sources = PLASMODE_SOURCES)
  assign("PLASMODE_CACHE", cache, envir = globalenv())
  sources <- PLASMODE_SOURCES
  sources$standin <- standin_source
  assign("PLASMODE_SOURCES", sources, envir = globalenv())
  rm(list = ls(.plasmode_memo), envir = .plasmode_memo)
  withr::local_envvar(
    REFUND_SHA = substr(REFUND_SHA_EXPECTED, 1, 8),
    .local_envir = env
  )
  withr::defer(
    {
      assign("PLASMODE_CACHE", old$cache, envir = globalenv())
      assign("PLASMODE_SOURCES", old$sources, envir = globalenv())
      rm(list = ls(.plasmode_memo), envir = .plasmode_memo)
    },
    envir = env
  )
  plasmode_cells(STANDIN)
}

test_that("the original 16 cells come first and unchanged; the extension is appended", {
  cells <- plasmode_cells()
  legacy <- cells[1:16, ]
  expect_equal(legacy$cell, 1:16)
  expect_equal(legacy$app, rep(c("ecg_gls4", "running"), each = 8))
  expect_equal(legacy$truth, rep(rep(c("REML", "NCV"), each = 4), 2))
  expect_equal(legacy$error, rep(rep(c("attached", "detached"), each = 2), 4))
  expect_equal(legacy$G, c(rep(c(79L, 40L), 4), rep(c(90L, 40L), 4)))
  expect_equal(legacy$extra, rep(c("ar1", "", "", ""), 4))
  expect_true(all(legacy$residual == legacy$truth & legacy$flip == "curve"))
  new <- cells[-(1:16), ]
  expect_equal(nrow(new), 122)
  expect_equal(new$cell, 16 + seq_len(122))
  expect_equal(
    c(table(new$app))[names(PLASMODE_EXTENSION)],
    c(
      ecg_gls4b = 24,
      af = 12,
      running = 8,
      dti = 12,
      gait_vc = 12,
      ecg_roahd = 12,
      ocean = 14,
      aemet = 14,
      adelaide = 14
    )
  )
  expect_false(
    anyDuplicated(cells[c("app", "truth", "residual", "flip", "error", "G")]) >
      0
  )
  expect_true(all(new$error == "attached"))
  # ECG 8-lead only: the wider intercept basis, in its cells and its truth fits
  expect_true(all(new$basis[new$app == "ecg_roahd"] == "ecg8"))
  expect_true(all(new$basis[new$app != "ecg_roahd"] == "default"))
  expect_equal(cell_bases(plasmode_fit_cell("ecg_roahd"))$int$k, 40)
  expect_null(cell_bases(plasmode_fit_cell("dti"))$int)
  # block arm: extra fit in every G = all cell of the block apps, never at G = 40
  has_block <- grepl("\\bblock\\b", new$extra)
  expect_equal(sum(has_block), 3 * 8)
  expect_true(all(new$G[has_block] > PLASMODE_G_SMALL))
  expect_true(all(
    new$app[new$flip == "block"] %in% c("ocean", "aemet", "adelaide")
  ))
  expect_equal(
    which(grepl("block_alt", new$extra)),
    which(new$app == "aemet" & new$flip == "block")
  )
  ar1 <- new[grepl("ar1", new$extra), c("app", "truth", "residual")]
  expect_equal(nrow(ar1), 10)
  expect_true(all(ar1$truth == "MID"))
  expect_setequal(ar1$residual[ar1$app == "ecg_gls4b"], c("REML", "beat_m"))
  # the running pool: 3 truths x 2 residual sources, old and new cells together
  expect_setequal(
    with(cells[plasmode_pool(cells, "running"), ], paste(truth, residual)),
    c("REML REML", "NCV NCV", "REML NCV", "NCV REML", "MID REML", "MID NCV")
  )
  expect_length(plasmode_pool(cells, "ecg_gls4b"), 12)
  # the pool excludes block-flip and G = 40 cells
  pool <- cells[plasmode_pool(cells, "ocean"), ]
  expect_equal(nrow(pool), 6)
  expect_true(all(pool$flip == "curve" & pool$G == 116))
})

test_that("draws: the original three are unchanged and block signs are appended", {
  n <- 78
  old <- withr::with_seed(
    SEED_BASE + 3,
    list(
      subset = sample.int(n),
      sign = sample(c(-1, 1), n, replace = TRUE),
      donor = sample.int(n)
    )
  )
  new <- plasmode_draws(n, 3)
  expect_identical(new[1:3], old)
  expect_named(new, c("subset", "sign", "donor", "block_sign"))
})

test_that("beat scaling matches the residuals' average variance", {
  res <- matrix(rnorm(200, sd = 2), 20)
  diff <- matrix(rnorm(200, sd = 0.5), 20)
  c <- beat_scale(res, diff)
  expect_equal(mean((c * diff)^2), mean(res^2))
})

test_that("MID: log-midpoint by name, interior check, and equal sp reproduce the fit", {
  sp_r <- c(a = 1, b = 100)
  sp_n <- c(b = 1, a = 100)
  expect_equal(mid_sp(sp_r, sp_n), c(a = 10, b = 10))
  expect_error(mid_sp(sp_r, c(a = 1, c = 2)), "align")
  terms <- data.frame(
    term = c("t1", "t2"),
    edf = c(5, 8),
    null = 2,
    full = 20,
    edf_reml = c(6, 10),
    edf_ncv = c(4, 6)
  )
  expect_true(mid_interior(terms))
  terms$edf[2] <- 11
  expect_false(mid_interior(terms))
  terms$edf[2] <- 8
  terms$edf_reml[1] <- 2.05 # at the null space: a selector at its boundary
  terms$edf[1] <- 3
  expect_false(mid_interior(terms))

  cells <- local_standin()
  d <- plasmode_data("standin")
  sim <- plasmode_sim(d, d$Y)
  reml <- fit_model(sim, PLASMODE_FIT_CELL, "REML")
  fixed <- fit_pffr(sim, PLASMODE_FIT_CELL, "REML", sp = reml$sp)
  expect_equal(unname(fit_sp(fixed)), unname(reml$sp))
  expect_equal(fixed$coefficients, reml$coefficients, tolerance = 1e-6)
})

test_that("stand-in app: truths, native grid, common random numbers, blocks", {
  cells <- local_standin()
  mine <- cells[cells$app == "standin", ]
  expect_equal(nrow(mine), 14)
  cell <- \(truth, residual, G = 56L, flip = "curve") {
    mine[
      mine$truth == truth &
        mine$residual == residual &
        mine$G == G &
        mine$flip == flip,
    ]
  }
  reml <- plasmode_truth("standin", "REML")
  expect_equal(reml$meta$G, 56L)
  expect_equal(reml$data$eval_t, reml$data$yind) # native D = 48: all points
  expect_length(reml$targets$mean, N_TEST * 48)
  mid <- plasmode_truth("standin", "MID")
  expect_true(is.logical(mid$mid$interior))
  expect_equal(
    unname(mid$meta$sp),
    unname(exp(
      (log(reml$meta$sp) +
        log(plasmode_truth("standin", "NCV")$meta$sp[names(reml$meta$sp)])) /
        2
    ))
  )
  # a stale cache is refused
  bad <- reml
  bad$meta$frame <- "other"
  saveRDS(bad, truth_file("standin", "REML"))
  rm(list = ls(.plasmode_memo), envir = .plasmode_memo)
  expect_error(plasmode_truth("standin", "REML"), "stale")
  saveRDS(reml, truth_file("standin", "REML"))
  rm(list = ls(.plasmode_memo), envir = .plasmode_memo)

  # common random numbers across residual sources: same signs, same rows
  a <- plasmode_simulate(cell("MID", "REML"), 5)
  b <- plasmode_simulate(cell("MID", "NCV"), 5)
  mu <- mid$mu
  e_r <- plasmode_truth("standin", "REML")$residuals
  e_n <- plasmode_truth("standin", "NCV")$residuals
  w_a <- rowMeans((unclass(a$data$Y) - mu) * e_r) / rowMeans(e_r^2)
  w_b <- rowMeans((unclass(b$data$Y) - mu) * e_n) / rowMeans(e_n^2)
  expect_equal(w_a, w_b)
  expect_true(all(abs(abs(w_a) - 1) < 1e-8))
  expect_equal(a$eval_t, mid$data$yind)
  # block ids only where the block arm runs (G = all), none at G = 40
  expect_identical(a$data$block, mid$data$block)
  small <- plasmode_simulate(cell("MID", "REML", G = 40L), 5)
  expect_null(small$data$block)
  expect_equal(nrow(small$data), 40)
  no_arm <- cell("MID", "REML")
  no_arm$extra <- ""
  expect_null(plasmode_simulate(no_arm, 5)$data$block)
  # block flips: one sign per block
  blk <- plasmode_simulate(cell("REML", "REML", flip = "block"), 5)
  w <- rowMeans((unclass(blk$data$Y) - reml$mu) * e_r) / rowMeans(e_r^2)
  expect_true(all(
    tapply(w, reml$data$block, \(x) length(unique(round(x, 8)))) == 1
  ))
})

test_that("a D = 48 block-arm task runs end to end through scoring", {
  cells <- local_standin()
  target <- cells[
    cells$app == "standin" &
      cells$truth == "REML" &
      cells$residual == "REML" &
      cells$flip == "block",
  ]
  target$extra <- "block" # skip nothing but AR(1)
  tasks <- data.frame(task = 1, cell = target$cell, rep = 1)
  out <- withr::local_tempdir()
  path <- run_task(1, out, target, tasks, simulate_fun = plasmode_simulate)
  record <- readRDS(path)
  expect_identical(record$status, "success")
  expect_length(record$extra_failures, 0)
  expect_true(all(
    c("se_cl2b_REML", "se_cl2_NCV_block", "se_cl2b_NCV_block") %in%
      colnames(record$records$mean)
  ))
  expect_equal(nrow(record$records$mean), N_TEST * 48)
  scores <- task_scores(list(record), target)
  expect_true(all(
    c("ncvb_cl2b_bias", "ncv_cl2b_bias", "reml_cl2b") %in% scores$arm
  ))
  expect_false("ncva_cl2a" %in% scores$arm)
  expect_true(all(is.finite(scores$coverage)))
})

test_that("a D = 101 stand-in is fitted and scored on every second grid point", {
  G <- 20L
  D <- 101L
  t <- seq(0, 1, length.out = D)
  d <- withr::with_seed(3, {
    X <- matrix(rnorm((G + N_TEST) * 5), G + N_TEST) %*%
      cos(outer(0:4, S_GRID * pi))
    z <- rnorm(G + N_TEST)
    list(X = X, z = z)
  })
  eta <- outer(rowMeans(d$X), sin(pi * t)) + outer(d$z, cos(pi * t))
  eval_t <- native_eval_t(t)
  expect_length(eval_t, 51)
  sim <- list(
    data = data.frame(
      X = I(d$X[1:G, ]),
      z = d$z[1:G],
      Y = I(eta[1:G, ] + withr::with_seed(4, matrix(rnorm(G * D, sd = 0.3), G)))
    ),
    yind = t,
    xind = S_GRID,
    eval_t = eval_t,
    test = data.frame(X = I(d$X[G + 1:N_TEST, ]), z = d$z[G + 1:N_TEST]),
    diag = list()
  )
  sim$truth <- list(
    mean = as.vector(eta[G + 1:N_TEST, seq(1, D, by = 2)]),
    beta = rep(0, nrow(eval_points("beta"))),
    alpha = rep(0, length(EVAL_T)),
    gamma = cos(pi * EVAL_T)
  )
  cell <- fixture_cell()
  cell$D <- D
  fits <- list(
    REML = fit_model(sim, cell, "REML"),
    NCV = fit_model(sim, cell, "NCV")
  )
  scored <- score_task(
    interval_inputs(fits, sim, cell),
    sim,
    cell,
    data.frame(task = 1, cell = cell$cell, rep = 1)
  )
  expect_equal(nrow(scored$records$mean), N_TEST * 51)
  s <- score_task_records(scored, cell)
  expect_gt(mean(s$coverage[s$estimand == "mean" & s$arm == "reml_cl2"]), 0.5)
  sim$eval_t <- eval_t + 0.001
  expect_error(evaluation_matrix(fits$REML, sim, "mean"), "subset")
})

test_that("validation keys: original and extension records aggregate together", {
  cells <- plasmode_cells()
  sim <- fixture_sim()
  inputs <- \(cell) fixture_inputs(sim, cell)
  record <- \(cell, spec) {
    c(
      list(status = "success", cell_spec = spec),
      score_task(
        inputs(cell),
        sim,
        cell,
        data.frame(task = 1, cell = cell$cell, rep = 1)
      )
    )
  }
  old_cell <- cells[cells$cell == 9, ]
  # the original records: no basis, residual or flip in the stored specification
  old_spec <- old_cell[setdiff(
    names(old_cell),
    c("basis", "covariate", "residual", "flip")
  )]
  new_cell <- cells[cells$app == "running" & cells$residual != cells$truth, ][
    1,
  ]
  mixed <- list(record(old_cell, old_spec), record(new_cell, new_cell))
  expect_setequal(unique(task_scores(mixed, cells)$cell), c(9, new_cell$cell))
  wrong <- new_cell
  wrong$residual <- new_cell$truth
  expect_error(task_scores(list(record(new_cell, wrong)), cells), "differs")
  wrong <- new_cell
  wrong$flip <- "block"
  expect_error(task_scores(list(record(new_cell, wrong)), cells), "differs")
  # synthetic records carry no app
  syn <- design_table()[1, ]
  expect_gt(nrow(task_scores(list(record(syn, syn)), design_table())), 0)
})

test_that("§4.1 decision table: one case per row, both U triggers, degenerate CE", {
  ok <- c(0.001, 0.015) # resolved adequate interval
  bad <- c(0.03, 0.06) # resolved inadequate interval
  mixed <- c(0.01, 0.03) # straddles 2 pp
  r <- c(0.8, 0.9)
  cls <- \(ce_p, ce_f, ratio, ci_p, ci_f, ci_r = r) {
    rule_class(ce_p, ce_f, ratio, list(ce_p = ci_p, ce_f = ci_f, ratio = ci_r))
  }
  # row 1: only P adequate
  expect_equal(cls(0.01, 0.04, 1, ok, bad), list(decision = "P", class = "P"))
  expect_equal(cls(0.01, 0.04, 1, mixed, bad)$class, "U")
  expect_equal(cls(0.01, 0.04, 1, ok, mixed)$class, "U")
  # row 2: only F adequate
  expect_equal(cls(0.04, 0.01, 1, bad, ok), list(decision = "F", class = "F"))
  expect_equal(cls(0.04, 0.01, 1, mixed, ok)$class, "U")
  expect_equal(cls(0.04, 0.01, 1, bad, mixed)$class, "U")
  # row 3: both adequate, ratio <= 0.95
  expect_equal(cls(0.01, 0.01, 0.9, ok, ok, c(0.85, 0.93))$class, "P")
  expect_equal(cls(0.01, 0.01, 0.9, ok, ok, c(0.85, 0.97))$class, "U")
  expect_equal(cls(0.01, 0.01, 0.9, mixed, ok, c(0.85, 0.93))$class, "U")
  # row 4: both adequate, ratio >= 1/0.95
  expect_equal(cls(0.01, 0.01, 1.1, ok, ok, c(1.06, 1.15))$class, "F")
  expect_equal(cls(0.01, 0.01, 1.1, ok, ok, c(1.04, 1.15))$class, "U")
  # row 5: both adequate, ratio in between
  expect_equal(cls(0.01, 0.01, 1, ok, ok, c(0.96, 1.04))$class, "E")
  expect_equal(cls(0.01, 0.01, 1, ok, ok, c(0.94, 1.04))$class, "U")
  expect_equal(cls(0.01, 0.01, 1, ok, ok, c(0.96, 1.06))$class, "U")
  # row 6: neither adequate: N, never U
  expect_equal(
    cls(0.04, 0.05, 1, mixed, mixed),
    list(decision = "N", class = "N")
  )
  # degenerate CE = 0: interval [0, upper] straddles 2 pp iff upper > 2 pp
  expect_equal(cls(0, 0.04, 1, c(0, 0.015), bad)$class, "P")
  expect_equal(cls(0, 0.04, 1, c(0, 0.025), bad)$class, "U")
})

test_that("§4.2 classifier: every label, 3 P + 3 E, N before U, question 1 classes", {
  expect_equal(consistency_label(c("P", "P")), "no consistency claim (k < 3)")
  expect_equal(
    consistency_label(c("P", "P", "P", "P", "U", "U")),
    "consistent for P"
  )
  expect_equal(
    consistency_label(c("F", "F", "F", "F", "E", "U")),
    "consistent for F"
  )
  expect_equal(
    consistency_label(c("P", "P", "P", "E", "E", "E")),
    "consistent for P-or-equivalent"
  )
  expect_equal(
    consistency_label(c("F", "E", "E", "E", "U", "U")),
    "consistent for F-or-equivalent"
  )
  expect_equal(
    consistency_label(c("E", "E", "E", "E", "U", "U")),
    "consistent for E"
  )
  expect_equal(consistency_label(c("P", "P", "P", "U", "U", "U")), "unresolved")
  expect_equal(
    consistency_label(c("N", "U", "U", "U", "P", "P")),
    "inconsistent (N)"
  )
  expect_equal(
    consistency_label(c("P", "P", "P", "P", "F", "E")),
    "inconsistent (mixed)"
  )
  expect_equal(
    consistency_label(c("F", "F", "P", "E", "E", "E")),
    "inconsistent (mixed)"
  )
  # m = ceiling(2k/3), u = floor(k/3): k = 7 needs 5
  expect_equal(
    consistency_label(c(rep("P", 4), "E", "E", "E")),
    "consistent for P-or-equivalent"
  )
  expect_equal(consistency_label(c(rep("P", 5), "U", "U")), "consistent for P")
  expect_equal(
    consistency_label(c("C", "C", "C", "C", "E1", "U1")),
    "consistent for P"
  )
  expect_equal(
    consistency_label(c("M", "M", "M", "M", "U1", "U1")),
    "consistent for F"
  )
  expect_match(consistency_label(c("P", "P", NA, "P")), "incomplete")
  expect_equal(q1_class(0.025, 0.05), "C")
  expect_equal(q1_class(-0.05, -0.025), "M")
  expect_equal(q1_class(-0.01, 0.015), "E1")
  expect_equal(q1_class(0.01, 0.03), "U1")
})
