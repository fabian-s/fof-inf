#!/bin/bash
# fastFGEE comparison (competitors/fgee/run.R). Each array task is one worker that
# claims tasks dynamically until none are left. Submit from $PROJECT_DIR/code-fgee:
#   sbatch --array=1-20 competitors/fgee/fgee.slurm A
#   sbatch --array=1-40 competitors/fgee/fgee.slurm B
#SBATCH -J fgee
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500mb
#SBATCH --time=12:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache" # experiment P reads the cached truths

Rscript competitors/fgee/run.R "$1" "$SLURM_ARRAY_TASK_ID" "$SLURM_ARRAY_TASK_COUNT" \
  "../fgee-results/$1" "${2:-200}"
