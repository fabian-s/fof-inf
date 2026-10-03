# Basis-size extension (2026-09-28): design, nesting of the bases, fits and scoring.

test_that("cells 1-111 are unchanged and the basis block is appended", {
  d <- design_table()
  frozen <- readr::read_csv(
    "../summaries/final/cells.csv",
    show_col_types = FALSE,
    col_types = readr::cols(extra = "c", .default = readr::col_guess())
  ) |>
    dplyr::distinct(
      cell,
      block,
      family,
      fit_family,
      error,
      G,
      signal,
      truth,
      D,
      model,
      covariate,
      extra
    ) |>
    dplyr::filter(cell <= 111) |>
    dplyr::arrange(cell)
  old <- d[d$cell <= 111, names(frozen)]
  old$extra[old$extra == ""] <- NA # read_csv reads empty strings as NA
  expect_equal(nrow(frozen), 111)
  expect_equal(as.data.frame(old), as.data.frame(frozen), ignore_attr = TRUE)
  expect_true(all(d$basis[d$cell <= 111] == "default"))
  # cells 112-123; later blocks (small_G, null_beta, 2026-10-03) are appended after them
  new <- d[d$cell > 111 & d$cell <= 123, ]
  expect_equal(nrow(new), 12)
  expect_true(all(new$block == "basis_size"))
  expect_setequal(new$basis, c("large", "xlarge"))
  expect_true(all(new$G == 100L & new$signal == "mid" & new$model == "M1"))
  expect_true(all(new$extra == ""))
  # every basis cell has exactly one default-basis core cell to be compared with
  source("../analysis/contrasts.R", local = TRUE)
  for (id in new$cell) {
    core <- matching_core_cell(d, id)
    expect_length(core, 1)
    expect_identical(d$basis[d$cell == core], "default")
    same <- c("family", "error", "G", "signal", "truth", "D", "covariate")
    expect_equal(d[d$cell == core, same], new[new$cell == id, same])
  }
})

test_that("default bases are the study's and the larger ones are nested", {
  expect_identical(BASES$default, list(t = BASIS_T, ff = BASIS_FF))
  margins <- \(b) {
    list(
      t = b$t,
      s = list(bs = b$ff$bs, k = b$ff$k[1], m = b$ff$m[[1]]),
      ff_t = list(bs = b$ff$bs, k = b$ff$k[2], m = b$ff$m[[2]])
    )
  }
  u <- seq(0, 1, length.out = 1001)
  residual <- \(small, big) {
    A <- ps_eval(ps_basis(small), u)
    B <- ps_eval(ps_basis(big), u)
    max(abs(A - B %*% qr.solve(B, A)))
  }
  pairs <- list(c("default", "large"), c("large", "xlarge"))
  for (p in pairs) {
    small <- margins(BASES[[p[1]]])
    big <- margins(BASES[[p[2]]])
    for (m in names(small)) {
      expect_lt(residual(small[[m]], big[[m]]), 1e-10)
      # the interval count doubles: the smaller basis' knots on [0, 1] (widened by
      # mgcv by 0.1% on each side) are a subset of the larger basis' knots
      k_small <- ps_basis(small[[m]])$knots
      k_small <- k_small[k_small > -0.01 & k_small < 1.01]
      nearest <- apply(
        abs(outer(k_small, ps_basis(big[[m]])$knots, "-")),
        1,
        min
      )
      expect_lt(max(nearest), 1e-12)
    }
  }
  # control: non-nested dimensions are not representable
  expect_gt(residual(BASIS_T, modifyList(BASIS_T, list(k = 20))), 1e-3)
})

test_that("a large-basis fit uses the cell's bases with the same knots", {
  sim <- fixture_sim()
  cell <- fixture_cell()
  cell$basis <- "large"
  fit <- fit_model(sim, cell, "REML")
  alpha <- estimand_smooth(fit, "alpha")
  gamma <- estimand_smooth(fit, "gamma")
  beta <- estimand_smooth(fit, "beta")
  expect_equal(alpha$knots, ps_basis(BASES$large$t)$knots)
  expect_equal(gamma$knots, ps_basis(BASES$large$t)$knots)
  expect_equal(
    beta$margin[[1]]$knots,
    ps_basis(list(bs = "ps", k = 13, m = c(2, 1)))$knots
  )
  expect_equal(
    beta$margin[[2]]$knots,
    ps_basis(list(bs = "ps", k = 17, m = c(2, 1)))$knots
  )
  expect_equal(beta$last.para - beta$first.para + 1, 13 * 17)
  input <- interval_inputs(list(REML = fit), sim, cell)$REML
  expect_true(all(is.finite(input$estimands$beta$se_cl2)))
  # cells without a basis column (plasmode) get the default; M2 has no larger bases
  expect_identical(
    cell_bases(fixture_cell()[setdiff(names(cell), "basis")]),
    BASES$default
  )
  m2 <- fixture_cell("M2")
  m2$basis <- "large"
  expect_error(cell_bases(m2), "M1 only")
  cell$basis <- "huge"
  expect_error(cell_bases(cell), "unknown basis")
})

test_that("scoring checks the basis and accepts records from before the column", {
  cells <- design_table()
  sim <- fixture_sim()
  # the fixture inputs plus the frequentist CL2 SEs that scoring expects
  inputs <- \(cell) {
    lapply(fixture_inputs(sim, cell), \(input) {
      input$estimands <- lapply(input$estimands, \(e) {
        c(e, list(se_cl2f = e$se_cl2))
      })
      input
    })
  }
  cell <- cells[cells$block == "basis_size", ][1, ]
  scored <- score_task(
    inputs(cell),
    sim,
    cell,
    data.frame(task = 1, cell = cell$cell, rep = 1)
  )
  record <- c(list(status = "success", cell_spec = cell), scored)
  expect_gt(nrow(task_scores(list(record), cells)), 0)
  record$cell_spec$basis <- "default"
  expect_error(task_scores(list(record), cells), "differs")
  old_cell <- cells[1, ]
  old <- c(
    list(
      status = "success",
      cell_spec = old_cell[setdiff(names(old_cell), "basis")]
    ),
    score_task(
      inputs(old_cell),
      sim,
      old_cell,
      data.frame(task = 1, cell = 1, rep = 1)
    )
  )
  expect_gt(nrow(task_scores(list(old), cells)), 0)
})
