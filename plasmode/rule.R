# The plasmode analysis: the §5 rule per dataset,
# the cross-dataset consistency classifier, the truth x residual attribution and the
# block arm. Needs scripts/setup.R, analysis/contrasts.R and plasmode/plasmode.R.
# All uncertainty is the replicate bootstrap of boot_stat() (replicates resampled jointly
# across the cells and arms of a statistic).

# Datasets counted in the consistency criterion (independent curves) and stress tests.
PLASMODE_COUNTED <- c(
  "ecg_gls4b",
  "running",
  "dti",
  "gait_vc",
  "ecg_roahd",
  "af"
)
PLASMODE_STRESS <- c("ocean", "aemet", "adelaide")
Q1_ARMS <- c("reml_model", "reml_cl2", "ncv_model", "ncv_cl2")
ATTRIBUTION_ARMS <- c(Q1_ARMS, "ncv_cl2_bias")
SMALL_COVERAGE <- 0.02 # |coverage effect| that counts as small (§4.3)
SMALL_LOG_MSE <- log(c(0.95, 1.05)) # band of a small log MSE-ratio effect
BLOCK_GAIN <- 0.9 # MSE ratio NCV_block / NCV that counts as help (§2.5)
BLOCK_COST_MSE <- 1.05 # MSE ratio of acceptable cost under curve flips
BLOCK_COST_COVERAGE <- 0.02 # coverage loss of acceptable cost under curve flips
# Block-arm recipes (§2.5): standard, block, the two cross combinations, REML + CL2.
BLOCK_RECIPES <- c(
  standard = "ncv_cl2_bias",
  block = "ncvb_cl2b_bias",
  ncv_blockcl2 = "ncv_cl2b_bias",
  blockncv_curvecl2 = "ncvb_cl2_bias",
  reml_curvecl2 = "reml_cl2",
  reml_blockcl2 = "reml_cl2b"
)
MISMATCH_RECIPES <- c(mismatch = "ncva_cl2a_bias", reml_altcl2 = "reml_cl2a")

# Pools ---------------------------------------------------------------------------------------

# An app's §4.1 pool: its attached, curve-flip cells with all subjects, optionally for
# some residual sources only (block-flip, G = 40 and detached cells are reported beside).
plasmode_pool <- function(cells, app, residuals = NULL) {
  mine <- cells[cells$app == app, ]
  keep <- mine$G == max(mine$G) &
    mine$error == "attached" &
    mine$flip == "curve"
  if (!is.null(residuals)) keep <- keep & mine$residual %in% residuals
  mine$cell[keep]
}

# Task scores of a pool restricted to replicates with every needed arm scored (complete
# sets, as paired_contrast()).
complete_sets <- function(tasks, pool, estimand, arms) {
  sub <- dplyr::filter(
    tasks,
    cell %in% pool,
    estimand == !!estimand,
    arm %in% arms
  )
  complete <- sub |>
    dplyr::count(cell, rep) |>
    dplyr::filter(n == length(arms))
  dplyr::semi_join(sub, complete, by = c("cell", "rep"))
}

# Decision table (§4.1) -------------------------------------------------------------------------

# Point decision from the calibration errors and the interval-score ratio.
rule_decision <- function(ce_p, ce_f, ratio) {
  adequate <- c(ce_p, ce_f) <= ADEQUATE
  if (adequate[1] && !adequate[2]) return("P")
  if (!adequate[1] && adequate[2]) return("F")
  if (!any(adequate)) return("N")
  if (ratio <= PRACTICAL_RATIO) return("P")
  if (ratio >= 1 / PRACTICAL_RATIO) return("F")
  "E"
}

# Class: the decision, or U when a criterion that decides it is unresolved (its interval
# straddles its threshold). `ci` holds c(lower, upper) for ce_p, ce_f and ratio. N is never
# resolved into P/F and never U.
rule_class <- function(ce_p, ce_f, ratio, ci) {
  decision <- rule_decision(ce_p, ce_f, ratio)
  if (decision == "N") return(list(decision = "N", class = "N"))
  st <- \(x, threshold) straddles(x[1], x[2], threshold)
  ce_unresolved <- st(ci$ce_p, ADEQUATE) || st(ci$ce_f, ADEQUATE)
  both <- ce_p <= ADEQUATE && ce_f <= ADEQUATE
  ratio_unresolved <- both &&
    switch(
      decision,
      P = st(ci$ratio, PRACTICAL_RATIO),
      F = st(ci$ratio, 1 / PRACTICAL_RATIO),
      E = st(ci$ratio, PRACTICAL_RATIO) || st(ci$ratio, 1 / PRACTICAL_RATIO)
    )
  unresolved <- ce_unresolved || isTRUE(ratio_unresolved)
  list(decision = decision, class = if (unresolved) "U" else decision)
}

# Question 1 per fit (§4.1): d = CE(model-based) - CE(CL2).
q1_class <- function(lower, upper) {
  if (lower > ADEQUATE) return("C")
  if (upper < -ADEQUATE) return("M")
  if (lower > -ADEQUATE && upper < ADEQUATE) return("E1")
  "U1"
}

# Cross-dataset classifier (§4.2), first matching label wins. `classes`: one of P, F, E,
# N, U per counted dataset (question 1: C, M, E1, U1, mapped to P, F, E, U).
consistency_label <- function(classes) {
  if (anyNA(classes)) return("incomplete (a dataset's pool is incomplete)")
  classes <- unname(
    c(C = "P", M = "F", E1 = "E", U1 = "U")[classes] %|% classes
  )
  k <- length(classes)
  if (k < 3) return("no consistency claim (k < 3)")
  m <- ceiling(2 * k / 3)
  u <- floor(k / 3)
  n <- \(label) sum(classes == label)
  if (n("N") > 0) return("inconsistent (N)")
  if (n("U") > u) return("unresolved")
  if (n("P") >= m && n("F") == 0) return("consistent for P")
  if (n("F") >= m && n("P") == 0) return("consistent for F")
  if (n("P") + n("E") >= m && n("P") >= 1 && n("F") == 0) {
    return("consistent for P-or-equivalent")
  }
  if (n("F") + n("E") >= m && n("F") >= 1 && n("P") == 0) {
    return("consistent for F-or-equivalent")
  }
  if (n("E") >= m && n("P") == 0 && n("F") == 0) return("consistent for E")
  "inconsistent (mixed)"
}
`%|%` <- function(x, y) ifelse(is.na(x), y, x)

# §4.1 per dataset -------------------------------------------------------------------------------

# The rule and question 1 on one pool of one app and estimand, with bootstrap intervals.
plasmode_rule <- function(tasks, pool, estimand) {
  arms <- unique(c(PROPOSAL, FALLBACK, Q1_ARMS))
  sub <- complete_sets(tasks, pool, estimand, arms)
  if (!setequal(unique(sub$cell), pool)) {
    return(tibble::tibble(estimand = estimand, decision = "incomplete pool"))
  }
  stat <- \(m) {
    m <- m[order(m$cell), ] # base indexing: this runs once per bootstrap draw
    ce <- \(a) {
      i <- m$arm == a
      calibration_error(m$coverage[i], m$coverage_se[i])
    }
    score <- \(a) m$score[m$arm == a]
    c(
      ce_p = ce(PROPOSAL),
      ce_f = ce(FALLBACK),
      ratio = exp(mean(log(score(PROPOSAL) / score(FALLBACK)))),
      q1_reml = ce("reml_model") - ce("reml_cl2"),
      q1_ncv = ce("ncv_model") - ce("ncv_cl2")
    )
  }
  b <- boot_stat(sub, stat)
  get <- \(name, col) b[[col]][b$name == name]
  ci <- \(name) c(get(name, "lower"), get(name, "upper"))
  cls <- rule_class(
    get("ce_p", "estimate"),
    get("ce_f", "estimate"),
    get("ratio", "estimate"),
    list(ce_p = ci("ce_p"), ce_f = ci("ce_f"), ratio = ci("ratio"))
  )
  tibble::tibble(
    estimand = estimand,
    n_cells = length(pool),
    n_rep_min = min(dplyr::count(dplyr::distinct(sub, cell, rep), cell)$n),
    ce_proposal = get("ce_p", "estimate"),
    ce_proposal_lower = get("ce_p", "lower"),
    ce_proposal_upper = get("ce_p", "upper"),
    ce_fallback = get("ce_f", "estimate"),
    ce_fallback_lower = get("ce_f", "lower"),
    ce_fallback_upper = get("ce_f", "upper"),
    score_ratio = get("ratio", "estimate"),
    score_ratio_lower = get("ratio", "lower"),
    score_ratio_upper = get("ratio", "upper"),
    decision = cls$decision,
    class = cls$class,
    q1_reml = get("q1_reml", "estimate"),
    q1_reml_lower = get("q1_reml", "lower"),
    q1_reml_upper = get("q1_reml", "upper"),
    q1_class_reml = q1_class(get("q1_reml", "lower"), get("q1_reml", "upper")),
    q1_ncv = get("q1_ncv", "estimate"),
    q1_ncv_lower = get("q1_ncv", "lower"),
    q1_ncv_upper = get("q1_ncv", "upper"),
    q1_class_ncv = q1_class(get("q1_ncv", "lower"), get("q1_ncv", "upper"))
  )
}

# The rule for every app and primary estimand (plus gamma where z is binary), on the full
# pool and, for the attribution's flip check, on the REML- and NCV-residual pools alone.
rule_table <- function(tasks, cells, apps = unique(cells$app)) {
  gamma_apps <- c("af", "dti", "ecg_roahd")
  rows <- list()
  for (app in intersect(apps, c(PLASMODE_COUNTED, PLASMODE_STRESS))) {
    estimands <- c("mean", "beta", if (app %in% gamma_apps) "gamma")
    pools <- list(
      all = plasmode_pool(cells, app),
      REML = plasmode_pool(cells, app, "REML"),
      NCV = plasmode_pool(cells, app, "NCV")
    )
    for (e in estimands) {
      for (p in names(pools)) {
        rows[[length(rows) + 1]] <- plasmode_rule(tasks, pools[[p]], e) |>
          dplyr::mutate(app = app, pool = p, .before = 1)
      }
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      role = dplyr::if_else(
        app %in% PLASMODE_COUNTED,
        "counted",
        "stress test"
      ),
      .after = app
    )
}

# §4.2 per primary estimand over the counted datasets present, for the rule and for
# question 1 per fit.
consistency_table <- function(rules) {
  counted <- dplyr::filter(
    rules,
    pool == "all",
    role == "counted",
    estimand %in% c("mean", "beta")
  )
  if (!nrow(counted) || !"class" %in% names(counted)) return(tibble::tibble())
  counted |>
    dplyr::group_by(estimand) |>
    dplyr::summarise(
      k = dplyr::n(),
      classes = paste(app, class, sep = ":", collapse = " "),
      rule = consistency_label(class),
      q1_reml = consistency_label(q1_class_reml),
      q1_ncv = consistency_label(q1_class_ncv),
      .groups = "drop"
    )
}

# §4.3 attribution -------------------------------------------------------------------------------

# Paired truth x residual effects on the 2 x 2 of REML/NCV truths and residuals at G = all
# (curve flips): per-cell values v, truth effect = mean over residuals of v(NCV truth) -
# v(REML truth), residual effect likewise, interaction = difference of differences.
attribution_2x2 <- function(tasks, cells, app, estimand) {
  pool <- plasmode_pool(cells, app, c("REML", "NCV"))
  pool <- pool[cells$truth[match(pool, cells$cell)] %in% c("REML", "NCV")]
  if (length(pool) != 4) return(NULL)
  at <- \(truth, residual) {
    pool[
      cells$truth[match(pool, cells$cell)] == truth &
        cells$residual[match(pool, cells$cell)] == residual
    ]
  }
  cell_of <- c(
    RR = at("REML", "REML"),
    RN = at("REML", "NCV"),
    NR = at("NCV", "REML"),
    NN = at("NCV", "NCV")
  )
  arms <- unique(c(ATTRIBUTION_ARMS, PROPOSAL, FALLBACK))
  sub <- complete_sets(tasks, pool, estimand, arms)
  if (!setequal(unique(sub$cell), pool)) return(NULL)
  effects <- \(v) {
    c(
      truth = ((v[["NR"]] + v[["NN"]]) - (v[["RR"]] + v[["RN"]])) / 2,
      residual = ((v[["RN"]] + v[["NN"]]) - (v[["RR"]] + v[["NR"]])) / 2,
      interaction = (v[["NN"]] - v[["NR"]]) - (v[["RN"]] - v[["RR"]])
    )
  }
  stat <- \(m) {
    value <- \(arm, metric) {
      vapply(
        cell_of,
        \(c) m[[metric]][m$cell == c & m$arm == arm],
        numeric(1)
      )
    }
    ce <- \(arm) {
      vapply(
        cell_of,
        \(c) {
          x <- m[m$cell == c & m$arm == arm, ]
          calibration_error(x$coverage, x$coverage_se)
        },
        numeric(1)
      )
    }
    out <- c(
      unlist(lapply(
        setNames(ATTRIBUTION_ARMS, paste0("coverage:", ATTRIBUTION_ARMS)),
        \(a) effects(value(a, "coverage"))
      )),
      log_mse_ratio = effects(log(
        value("ncv_cl2", "sq_error") / value("reml_cl2", "sq_error")
      )),
      ce_diff = effects(ce(PROPOSAL) - ce(FALLBACK)),
      log_is_ratio = effects(log(
        value(PROPOSAL, "score") / value(FALLBACK, "score")
      ))
    )
    out
  }
  b <- boot_stat(sub, stat)
  parts <- strsplit(b$name, ".", fixed = TRUE)
  b$quantity <- vapply(parts, \(p) paste(head(p, -1), collapse = "."), "")
  b$effect <- vapply(parts, \(p) tail(p, 1), "")
  b |>
    dplyr::mutate(
      app = app,
      estimand = estimand,
      verdict = smallness(quantity, effect, lower, upper)
    ) |>
    dplyr::select(
      app,
      estimand,
      quantity,
      effect,
      estimate,
      se,
      lower,
      upper,
      verdict
    )
}

# "small" / "not small" / "undetermined" for residual effects and interactions (§4.3);
# NA for truth effects and for quantities without a threshold.
smallness <- function(quantity, effect, lower, upper) {
  band <- \(lo, hi) {
    dplyr::case_when(
      lower > lo & upper < hi ~ "small",
      upper < lo | lower > hi ~ "not small",
      TRUE ~ "undetermined"
    )
  }
  coverage <- band(-SMALL_COVERAGE, SMALL_COVERAGE)
  mse <- band(SMALL_LOG_MSE[1], SMALL_LOG_MSE[2])
  dplyr::case_when(
    effect == "truth" ~ NA_character_,
    startsWith(quantity, "coverage:") ~ coverage,
    quantity == "log_mse_ratio" ~ mse,
    TRUE ~ NA_character_
  )
}

# Beat cells (ECG): per truth, coverage of every attribution arm under (truth, beat_m) and
# (truth, beat_1) minus (truth, REML), paired by replicate.
beat_contrasts <- function(tasks, cells, app, estimand) {
  mine <- cells[cells$app == app & cells$cell %in% plasmode_pool(cells, app), ]
  if (!any(grepl("^beat", mine$residual))) return(NULL)
  rows <- list()
  for (truth in unique(mine$truth)) {
    ref <- mine$cell[mine$truth == truth & mine$residual == "REML"]
    for (beat in c("beat_m", "beat_1")) {
      other <- mine$cell[mine$truth == truth & mine$residual == beat]
      sub <- complete_sets(tasks, c(ref, other), estimand, ATTRIBUTION_ARMS)
      if (!setequal(unique(sub$cell), c(ref, other))) next
      stat <- \(m) {
        setNames(
          vapply(
            ATTRIBUTION_ARMS,
            \(a) {
              m$coverage[m$cell == other & m$arm == a] -
                m$coverage[m$cell == ref & m$arm == a]
            },
            numeric(1)
          ),
          paste0("coverage:", ATTRIBUTION_ARMS)
        )
      }
      rows[[length(rows) + 1]] <- boot_stat(sub, stat) |>
        dplyr::transmute(
          app = app,
          estimand = estimand,
          quantity = name,
          effect = paste0(beat, " - REML residuals (truth ", truth, ")"),
          estimate,
          se,
          lower,
          upper,
          verdict = NA_character_
        )
    }
  }
  dplyr::bind_rows(rows)
}

attribution_table <- function(tasks, cells, apps = unique(cells$app)) {
  apps <- intersect(apps, c(PLASMODE_COUNTED, PLASMODE_STRESS))
  dplyr::bind_rows(lapply(apps, \(app) {
    dplyr::bind_rows(lapply(c("mean", "beta"), \(e) {
      dplyr::bind_rows(
        attribution_2x2(tasks, cells, app, e),
        beat_contrasts(tasks, cells, app, e)
      )
    }))
  }))
}

# Flips (§4.3): the class on the REML-residual pool differs from the NCV-residual pool's.
flip_table <- function(rules) {
  if (!"class" %in% names(rules)) return(tibble::tibble())
  rules |>
    dplyr::filter(pool %in% c("REML", "NCV")) |>
    dplyr::select(app, estimand, pool, class) |>
    tidyr::pivot_wider(names_from = pool, values_from = class) |>
    dplyr::mutate(flips = REML != NCV)
}

# §4.5 block arm -------------------------------------------------------------------------------

# Per block app and estimand: under block flips (per truth and pooled) the paired MSE
# ratio NCV_block / NCV and the calibration error and interval score of each recipe; under
# curve flips (pooled over the G = all curve-flip cells) the cost: MSE ratio and the
# coverage loss of the block recipe against the standard one; and the verdict.
block_arm <- function(tasks, cells, apps = PLASMODE_STRESS) {
  rows <- list()
  for (app in intersect(intersect(apps, PLASMODE_STRESS), unique(cells$app))) {
    mine <- cells[
      cells$app == app & cells$G == max(cells$G[cells$app == app]),
    ]
    block_cells <- mine$cell[mine$flip == "block"]
    if (!length(block_cells)) next
    curve_cells <- mine$cell[mine$flip == "curve" & grepl("block", mine$extra)]
    mismatch <- any(grepl("block_alt", mine$extra))
    recipes <- c(BLOCK_RECIPES, if (mismatch) MISMATCH_RECIPES)
    for (e in c("mean", "beta")) {
      sets <- c(
        list(block_flips = block_cells, curve_flips = curve_cells),
        setNames(
          as.list(block_cells),
          paste0("block_flips:", mine$truth[match(block_cells, mine$cell)])
        )
      )
      # the alternative partition's arms exist in the block-flip cells only
      out <- lapply(names(sets), \(s) {
        alt <- mismatch && startsWith(s, "block_flips")
        rec <- if (alt) recipes else BLOCK_RECIPES
        block_arm_set(tasks, sets[[s]], e, rec, alt)
      })
      names(out) <- names(sets)
      out <- out[!vapply(out, is.null, logical(1))]
      if (!length(out)) next
      stats <- dplyr::bind_rows(out, .id = "cells")
      rows[[length(rows) + 1]] <- stats |>
        dplyr::mutate(
          app = app,
          estimand = e,
          verdict = block_verdict(stats),
          .before = 1
        )
    }
  }
  dplyr::bind_rows(rows)
}

block_arm_set <- function(tasks, pool, estimand, recipes, mismatch) {
  mse_arms <- c("ncvb_cl2", "ncv_cl2", if (mismatch) "ncva_cl2")
  sub <- complete_sets(tasks, pool, estimand, unique(c(recipes, mse_arms)))
  if (!length(pool) || !setequal(unique(sub$cell), pool)) return(NULL)
  stat <- \(m) {
    per_cell <- \(arm, metric) {
      vapply(pool, \(c) m[[metric]][m$cell == c & m$arm == arm], numeric(1))
    }
    ratio <- \(a, b)
      exp(mean(log(per_cell(a, "sq_error") / per_cell(b, "sq_error"))))
    ce <- \(arm) {
      x <- m[m$arm == arm & m$cell %in% pool, ]
      calibration_error(x$coverage, x$coverage_se)
    }
    c(
      mse_ratio_block = ratio("ncvb_cl2", "ncv_cl2"),
      mse_ratio_mismatch = if (mismatch) ratio("ncva_cl2", "ncv_cl2"),
      # cost (§2.5): coverage lost by clustering CL2 on blocks instead of curves, within
      # the same fit (the larger of the NCV and block-NCV fits'), and between the recipes
      coverage_loss = max(
        mean(
          per_cell("ncv_cl2_bias", "coverage") -
            per_cell("ncv_cl2b_bias", "coverage")
        ),
        mean(
          per_cell("ncvb_cl2_bias", "coverage") -
            per_cell("ncvb_cl2b_bias", "coverage")
        )
      ),
      coverage_loss_recipe = mean(
        per_cell(recipes[["standard"]], "coverage") -
          per_cell(recipes[["block"]], "coverage")
      ),
      setNames(
        vapply(recipes, \(a) mean(per_cell(a, "coverage")), numeric(1)),
        paste0("coverage:", names(recipes))
      ),
      setNames(
        vapply(recipes, ce, numeric(1)),
        paste0("ce:", names(recipes))
      ),
      setNames(
        vapply(recipes, \(a) mean(per_cell(a, "score")), numeric(1)),
        paste0("score:", names(recipes))
      )
    )
  }
  boot_stat(sub, stat)
}

# §2.5 verdict from the pooled block-flip gain and the curve-flip cost: "helps" (gain at
# acceptable cost), "cost only" (no gain, acceptable cost), "hurts" (cost exceeded), or
# "inconclusive" when an interval straddles its threshold.
block_verdict <- function(stats) {
  get <- \(cells, name) stats[stats$cells == cells & stats$name == name, ]
  gain <- get("block_flips", "mse_ratio_block")
  cost_mse <- get("curve_flips", "mse_ratio_block")
  cost_cov <- get("curve_flips", "coverage_loss")
  if (!nrow(gain) || !nrow(cost_mse) || !nrow(cost_cov)) return("incomplete")
  gains <- if (gain$estimate <= BLOCK_GAIN && gain$upper < 1) {
    TRUE
  } else if (gain$lower > BLOCK_GAIN || gain$lower >= 1) {
    FALSE
  } else {
    NA
  }
  within <- \(x, bound) {
    if (straddles(x$lower, x$upper, bound)) return(NA)
    x$estimate <= bound
  }
  cost_ok <- within(cost_mse, BLOCK_COST_MSE) &&
    within(cost_cov, BLOCK_COST_COVERAGE)
  if (isFALSE(cost_ok)) return("hurts")
  if (is.na(gains) || is.na(cost_ok)) return("inconclusive")
  if (gains) "helps" else "cost only"
}
