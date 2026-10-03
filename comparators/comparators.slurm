#!/bin/bash
# Comparator study workers: each array task is one worker that claims tasks dynamically
# (lock dirs) until none are left, so the array size is just the number of workers (<= 64).
# Submit from $PROJECT_DIR/code-comparators (frozen copy of the repository):
#   sbatch --array=1-4  comparators/comparators.slurm ../comparators-pilot s012,p081 2
#   sbatch --array=1-64 comparators/comparators.slurm ../comparators-results all 200 <boot_syn> <boot_plas>
# Arguments: <out_dir> <keys|all> <max_rep> [boot reps synthetic] [boot reps plasmode].
# Cardiology apps run from the de-identified input only (PLASMODE_CARDIO=deid); task files
# stay on LRZ.
#SBATCH -J comparators
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000mb
#SBATCH --time=24:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache" # the production truths
export PLASMODE_CARDIO=deid

out_dir="$1"
[ "$2" != "all" ] && export COMP_KEYS="$2"
export COMP_MAX_REP="${3:-200}"
[ -n "$4" ] && export COMP_BOOT_REPS_SYN="$4"
[ -n "$5" ] && export COMP_BOOT_REPS_PLAS="$5"

Rscript comparators/run-comparators.R "$out_dir" ../results/373f6fd65cac \
  ../plasmode-results
