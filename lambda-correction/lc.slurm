#!/bin/bash
# Lambda-correction workers: each array task is one worker claiming tasks dynamically
# (lock dirs). Submit from $PROJECT_DIR/code-lc (frozen copy of the repository):
#   sbatch --array=1-60 lambda-correction/lc.slurm            (all replicates)
#   sbatch --array=1-18 lambda-correction/lc.slurm ../lc-pilot 1:1   (pilot)
#SBATCH -J lambdacorr
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000mb
#SBATCH --time=12:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

OUT="${1:-../lc-results}"
REPS="${2:-1:200}"
Rscript lambda-correction/run-lc.R "$OUT" ../results/373f6fd65cac Inf "$REPS"
