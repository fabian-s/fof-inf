#!/bin/bash
# Satterthwaite df workers for the extension, small-G, null-effect and plasmode cells
# (2026-10-03): each array task is one worker that claims tasks dynamically (lock dirs).
# Submit from $PROJECT_DIR/code-r4 (frozen copy of the repository). Arguments:
# <out_dir> <syn|plasmode> <cells (R expression) | apps> <rep_from> <rep_to>
#   sbatch --array=1-40 satterthwaite/satt-ext.slurm ../satt-ext-results syn 25:152 1 200
#   sbatch --array=1-40 satterthwaite/satt-ext.slurm ../satt-pm-results plasmode af,dti,... 1 100
#SBATCH -J satt-ext
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000mb
#SBATCH --time=24:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
export PLASMODE_CARDIO=deid

# stored study records for the per-replicate reproduction check
if [ "$2" = "plasmode" ]; then
  export SATT_REF_DIRS="$PROJECT_DIR/plasmode-results"
else
  export SATT_REF_DIRS="$PROJECT_DIR/results:$PROJECT_DIR/results/373f6fd65cac"
fi

Rscript satterthwaite/run-satt-ext.R "$1" "$2" "$3" "$4" "$5"
