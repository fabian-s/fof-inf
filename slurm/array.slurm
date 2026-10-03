#!/bin/bash
# One run (R = 100) of the study, or of a subset of its cells, as a SLURM array. Each
# array task runs CHUNK tasks serially; chunk membership depends only on the task table
# (scripts/run-array.R). Arguments: chunk run [cells]. Submit from $PROJECT_DIR/code
# (copy of the repository). Optional 4th argument: output directory (default ../results).
#   sbatch --array=1-222 slurm/array.slurm 50 1           # all 111 cells: 11100 tasks / 50
#   sbatch --array=1-30 slurm/array.slurm 20 2 100:105    # cells 100-105, run 2: 600 / 20
#SBATCH -J restart
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=03:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export MGCV_EXPECT=1.9.5
export RESTART_CELLS="${3:-}" # optional subset, e.g. 106:111

Rscript scripts/run-array.R "$1" "$2" "${4:-../results}"
