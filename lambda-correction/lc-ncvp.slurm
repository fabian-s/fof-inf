#!/bin/bash
# Lambda-correction on the plasmode cells: each array task is one worker claiming tasks
# dynamically (lock dirs). Submit from $PROJECT_DIR/code-lc-ncvp (frozen rsync of
# the repository): sbatch --array=1-60 lambda-correction/lc-plasmode.slurm <out_dir> [reps] [apps] [max_tasks]
# Reads the plasmode study's truth caches and records; cardiology apps from the
# de-identified input only (PLASMODE_CARDIO=deid); ecg_gls4 is excluded by the worker.
#SBATCH -J lc-ncvp
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000mb
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

OUT="${1:-../lc-ncvp-results}"
REPS="${2:-1:100}"
APPS="${3:-}"
MAX="${4:-Inf}"
Rscript lambda-correction/run-lc-ncv-plasmode.R "$OUT" ../plasmode-results "$MAX" "$REPS" $APPS
