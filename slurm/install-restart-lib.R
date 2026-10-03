# Install the study's R library on LRZ. Run on the login node with the
# micromamba R 4.6.1 env on PATH (its compilers are needed for source builds):
#   PATH=$HOME/envs/r461/bin:$PATH TMPDIR=$HOME/tmp Rscript slurm/install-restart-lib.R
# Expects mgcv_1.9-5-patched.tar.gz (software/) and refund-79a346fb.tar.gz (a source
# tarball of refund at commit 79a346fb) in RESTART_SRC (default software/).
# mgcv: the 1.9-5 pre-release plus software/mgcv-1.9-5-ncv-vc-and-cxwyd.patch.
lib <- path.expand(Sys.getenv("RESTART_LIB", "~/R/restart-lib"))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4L)
src <- path.expand(Sys.getenv("RESTART_SRC", "software"))

missing_only <- function(pkgs) {
  pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE, lib.loc = lib)
  ]
}

# mgcv's own dependencies first, then the pre-release mgcv tarball.
install.packages(missing_only(c("nlme", "Matrix")), lib = lib)
install.packages(
  file.path(src, "mgcv_1.9-5-patched.tar.gz"),
  lib = lib,
  repos = NULL,
  type = "source"
)
stopifnot(packageVersion("mgcv", lib.loc = lib) == "1.9.5")

# refund's CRAN dependencies (mgcv already satisfied from lib).
deps <- c(
  "fda",
  "lattice",
  "boot",
  "MASS",
  "magic",
  "gamm4",
  "lme4",
  "RLRsim",
  "grpreg",
  "ggplot2",
  "pbs"
)
install.packages(missing_only(deps), lib = lib)
install.packages(
  file.path(src, "refund-79a346fb.tar.gz"),
  lib = lib,
  repos = NULL,
  type = "source"
)
stopifnot(packageVersion("mgcv", lib.loc = lib) == "1.9.5")
stopifnot(requireNamespace("refund", lib.loc = lib, quietly = TRUE))
cat("done\n")
