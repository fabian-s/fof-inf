#!/bin/bash
# fastFGEE plasmode comparison on the de-identified cardiology apps (ecg_gls4b, af), LRZ only.
# claims tasks dynamically until none are left. Submit from $PROJECT_DIR/code-fgee:
#   sbatch --array=1-20 competitors/fgee/fgee-cardio.slurm [max_rep]
#SBATCH -J fgee-cardio
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
export PLASMODE_CARDIO=deid FGEE_P_APPS=ecg_gls4b,af # de-identified input only; raw task files stay on LRZ

Rscript competitors/fgee/run.R P "$SLURM_ARRAY_TASK_ID" "$SLURM_ARRAY_TASK_COUNT" \
  "../fgee-results/Pc" "${1:-200}"
