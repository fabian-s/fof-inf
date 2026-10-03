# One-off: write every public plasmode dataset as a plain .rds (fitted grids, blocks,
# provenance, content hash) under PLASMODE_PUBLIC_DATA/<app>/data.rds, with a README per
# set. Runs locally (needs tf, roahd, fda.usc, fds, refund and the downloaded files); the
# written files are what plasmode/sources.R reads, locally and on LRZ. Public data only.
# Run from the repository root:  Rscript plasmode/export-public.R
source("scripts/setup.R")
source("plasmode/sources.R")

DOWNLOADED <- "2026-09-28" # date the files/package versions were obtained
# running.rda: the running data of the tidyfundata package (tidyfun project);
# ocean.RData: the Hawaii ocean profiles of the refund comment on penalized FoFR.
RUNNING_FILE <- Sys.getenv("PLASMODE_RUNNING_FILE", "running.rda")
OCEAN_FILE <- Sys.getenv("PLASMODE_OCEAN_FILE", "ocean.RData")
GAIT_DIR <- file.path(PLASMODE_PUBLIC_DATA, "vancriekinge2023")
OCEAN_BLOCK <- 3L # consecutive casts per block (a quarter)
ADELAIDE_BLOCK <- 2L # consecutive frame curves (10 weeks)
ADELAIDE_STEP <- 5L # every 5th Monday
AEMET_STEP <- 5L # every 5th day
AEMET_CUT_KM <- c(block = 250, block_alt = 150) # complete-linkage cuts
EARTH_RADIUS_KM <- 6371

equidistant <- \(n) seq(0, 1, length.out = n)

# Rows of a matrix interpolated from an equidistant grid on [0, 1] to `grid`.
regrid <- function(M, grid) {
  u <- equidistant(ncol(M))
  t(apply(M, 1, \(v) stats::approx(u, v, xout = grid)$y))
}

consecutive_blocks <- \(n, size) as.integer(ceiling(seq_len(n) / size))

great_circle_km <- function(lat, lon) {
  rad <- pi / 180
  a <- outer(lat, lat, \(p, q) sin((q - p) * rad / 2)^2) +
    outer(cos(lat * rad), cos(lat * rad)) *
      outer(lon, lon, \(p, q) sin((q - p) * rad / 2)^2)
  2 * EARTH_RADIUS_KM * asin(pmin(sqrt(a), 1))
}

# Loaders: raw public data on the fitted grids ---------------------------------------------

PUBLIC_LOADERS <- list(
  # Fukuchi et al. 2017 via tidyfundata; exactly the original plasmode_data("running").
  running = function() {
    e <- new.env()
    load(path.expand(RUNNING_FILE), envir = e)
    r <- e$running |>
      dplyr::filter(!grepl("_post_", id)) |>
      dplyr::group_by(id) |>
      dplyr::slice(1) |>
      dplyr::ungroup()
    list(
      Y = on_grid(r$knee_moment_ml, t_grid(61)),
      X = on_grid(r$knee_angle_ml, S_GRID),
      z = standardize(r$wt),
      yind = t_grid(61),
      ids = r$id,
      source = "tidyfundata::running (Fukuchi, Fukuchi & Duarte 2017, PeerJ 5:e3298); sagittal knee moment ~ knee angle + body mass, one row per person",
      licence = "tidyfundata: MIT; original data public (Fukuchi et al. 2017)"
    )
  },
  # refund::DTI baseline visit: rcst (55, native) ~ cca + case (MS = 1).
  dti = function() {
    utils::data("DTI", package = "refund", envir = environment())
    b <- DTI[DTI$visit == 1, ]
    ok <- stats::complete.cases(b$rcst) & stats::complete.cases(b$cca)
    b <- b[ok, ]
    list(
      Y = unclass(b$rcst),
      X = regrid(unclass(b$cca), S_GRID),
      z = as.numeric(b$case),
      yind = equidistant(ncol(b$rcst)),
      ids = as.character(b$ID),
      source = "refund::DTI, baseline visit (Goldsmith et al. 2011); rcst ~ cca + case",
      licence = "refund: GPL (>= 2)"
    )
  },
  # Van Criekinge et al. 2023 (figshare 24192489): knee moment (every 10th of 1001
  # points) ~ knee angle + age; subject order = sheets Sub01..Sub138 = sub_char rows.
  gait_vc = function() {
    curves <- readRDS(file.path(GAIT_DIR, "gait_vc-curves.rds"))
    chars <- utils::read.csv(file.path(GAIT_DIR, "sub_char.csv"))
    stopifnot(
      length(curves) == nrow(chars),
      identical(chars$sub, seq_len(nrow(chars)))
    )
    column <- \(name) t(vapply(curves, \(m) m[, name], numeric(1001)))
    list(
      Y = column("KneeMoment")[, seq(1, 1001, by = 10)],
      X = regrid(column("KneeAngles"), S_GRID),
      z = standardize(chars$Age),
      yind = equidistant(101),
      ids = sprintf("Sub%03d", chars$sub),
      source = "Van Criekinge et al. 2023, Sci Data 10:852, figshare 24192489 (Excel curves; age from the MAT file's sub_char); knee moment ~ knee angle (sagittal) + age",
      licence = "CC0"
    )
  },
  # roahd ECG: lead V5 (series 5) ~ lead I (series 7, presumed) + LBBB; 1024 -> 101.
  ecg_roahd = function() {
    groups <- list(healthy = roahd::mfD_healthy, lbbb = roahd::mfD_LBBB)
    lead <- \(k) do.call(rbind, lapply(groups, \(g) g$fDList[[k]]$values))
    n <- vapply(groups, \(g) g$N, integer(1))
    list(
      Y = regrid(lead(5), equidistant(101)),
      X = regrid(lead(7), S_GRID),
      z = rep(c(0, 1), n),
      yind = equidistant(101),
      ids = paste0(rep(names(groups), n), "_", sequence(n)),
      source = "roahd::mfD_healthy, roahd::mfD_LBBB (Ieva et al.); lead V5 ~ lead I (presumed) + LBBB",
      licence = "roahd: GPL-3"
    )
  },
  # Hawaii Ocean Time-series profiles (FRegSigCom::ocean): oxygen ~ temperature + cast order.
  ocean = function() {
    e <- new.env()
    load(path.expand(OCEAN_FILE), envir = e)
    o <- e$ocean
    G <- nrow(o$Oxygen)
    list(
      Y = unclass(o$Oxygen),
      X = regrid(unclass(o$Temperature), S_GRID),
      z = standardize(seq_len(G)),
      yind = equidistant(ncol(o$Oxygen)),
      ids = sprintf("cast%03d", seq_len(G)),
      block = consecutive_blocks(G, OCEAN_BLOCK),
      source = "FRegSigCom::ocean (Hawaii Ocean Time-series casts, monthly); oxygen ~ temperature + cast order",
      licence = "FRegSigCom: GPL (>= 2); HOT data public"
    )
  },
  # AEMET 1980-2009 station means (fda.usc::aemet): log precipitation (every 5th day) ~
  # temperature + altitude; blocks = complete-linkage clusters of great-circle distance.
  aemet = function() {
    utils::data("aemet", package = "fda.usc", envir = environment())
    days <- seq(1, 365, by = AEMET_STEP)
    dist <- stats::as.dist(great_circle_km(
      aemet$df$latitude,
      aemet$df$longitude
    ))
    tree <- stats::hclust(dist, method = "complete")
    list(
      Y = unclass(aemet$logprec$data)[, days],
      X = regrid(unclass(aemet$temp$data), S_GRID),
      z = standardize(aemet$df$altitude),
      yind = equidistant(length(days)),
      ids = make.unique(aemet$df$ind),
      block = stats::cutree(tree, h = AEMET_CUT_KM[["block"]]),
      block_alt = stats::cutree(tree, h = AEMET_CUT_KM[["block_alt"]]),
      source = "fda.usc::aemet (AEMET open data, station means 1980-2009); log precipitation ~ temperature + altitude",
      licence = "fda.usc: GPL-2"
    )
  },
  # Adelaide half-hourly demand, every 5th Monday (fds): demand ~ temperature + week index.
  adelaide = function() {
    utils::data(
      list = c("mondaydemand", "mondaytempairport"),
      package = "fds",
      envir = environment()
    )
    frame <- seq(1, ncol(mondaydemand$y), by = ADELAIDE_STEP)
    G <- length(frame)
    list(
      Y = t(unclass(mondaydemand$y)[, frame]),
      X = regrid(t(unclass(mondaytempairport$y)[, frame]), S_GRID),
      z = standardize(frame),
      yind = equidistant(nrow(mondaydemand$y)),
      ids = sprintf("monday%03d", frame),
      block = consecutive_blocks(G, ADELAIDE_BLOCK),
      source = "fds::mondaydemand, fds::mondaytempairport (Magnano, Boland & Hyndman 2008); every 5th Monday 1997-2007",
      licence = "fds: GPL-3"
    )
  }
)

# Export -------------------------------------------------------------------------------------

export_set <- function(app) {
  raw <- PUBLIC_LOADERS[[app]]()
  d <- list(
    Y = unname(raw$Y),
    X = unname(raw$X),
    z = raw$z,
    yind = raw$yind,
    xind = S_GRID,
    eval_t = native_eval_t(raw$yind),
    ids = as.character(raw$ids),
    block = raw$block,
    block_alt = raw$block_alt
  )
  check_plasmode_data(c(d[!vapply(d, is.null, logical(1))]))
  d$source <- raw$source
  d$licence <- raw$licence
  d$downloaded <- DOWNLOADED
  d$hash <- digest::digest(d[PUBLIC_HASHED_FIELDS])
  dir <- file.path(PLASMODE_PUBLIC_DATA, app)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(d, file.path(dir, "data.rds"))
  writeLines(
    c(
      paste("dataset:", app),
      paste("source:", d$source),
      paste("licence:", d$licence),
      paste("obtained:", d$downloaded),
      sprintf(
        "curves G = %d, fitted D = %d, covariate points %d",
        length(d$z),
        length(d$yind),
        length(d$xind)
      ),
      if (!is.null(d$block))
        sprintf(
          "blocks: %d (sizes %s)",
          length(unique(d$block)),
          toString(range(table(d$block)))
        ),
      if (!is.null(d$block_alt))
        sprintf("alternative blocks: %d", length(unique(d$block_alt))),
      paste("content hash (digest of Y, X, z, grids, ids, blocks):", d$hash),
      "written by plasmode/export-public.R"
    ),
    file.path(dir, "README.txt")
  )
  tibble::tibble(
    app = app,
    G = length(d$z),
    D = length(d$yind),
    n_eval_t = length(d$eval_t),
    n_block = length(unique(d$block %||% NA)),
    n_block_alt = length(unique(d$block_alt %||% NA)),
    hash = d$hash
  )
}

apps <- commandArgs(trailingOnly = TRUE)
if (!length(apps)) apps <- names(PUBLIC_LOADERS)
print(dplyr::bind_rows(lapply(apps, export_set)), width = Inf)
