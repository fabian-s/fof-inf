#!/bin/bash
# Aggregates of the 2026-10-03 runs (LRZ; plasmode task files never leave it):
# Satterthwaite on the extension/new synthetic cells and the plasmode cells, the G = 20 and
# null-effect cells, and their relative error. From $PROJECT_DIR/code-r4:
#   sbatch satterthwaite/r4-analyse.slurm
#SBATCH -J r4-analyse
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24000mb
#SBATCH --time=06:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
export PLASMODE_CARDIO=deid
export RESTART_CORES=8
set -e
# optional argument: the steps to run (default all): syn plasmode sgn relerr
steps="${1:-syn plasmode sgn relerr}"
for step in $steps; do
  case $step in
    syn) Rscript satterthwaite/analyse-satt-ext.R ../satt-r4-syn syn satterthwaite/summaries ;;
    plasmode) Rscript satterthwaite/analyse-satt-ext.R ../satt-r4-pm plasmode satterthwaite/summaries ;;
    sgn) Rscript small-g-null/analyse.R ../r4-results ../satt-r4-syn small-g-null/summaries ;;
    relerr) RELERR_DIRS=../r4-results Rscript analysis/relative-error.R synthetic small-g-null/summaries/relative-error.csv ;;
  esac
done
echo all done
