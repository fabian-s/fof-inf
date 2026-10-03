#!/bin/bash
# CL2-arms workers (GLS raw / fpc + CL2, pcre + CL2): each array task is one worker that
# claims tasks dynamically (lock dirs) until none are left.
# Submit from $PROJECT_DIR/code-cl2arms (frozen copy of the repository):
#   sbatch --array=1-2  comparators/cl2-arms.slurm ../cl2arms-pilot s009,p082 1
#   sbatch --array=1-30 comparators/cl2-arms.slurm ../cl2arms-results all 200 100
# Arguments: <out_dir> <keys|all> <max_rep> [pcre reps].
#SBATCH -J cl2arms
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000mb
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

out_dir="$1"
[ "$2" != "all" ] && export CL2_KEYS="$2"
export CL2_MAX_REP="${3:-200}"
[ -n "$4" ] && export CL2_PCRE_REPS="$4"

Rscript comparators/run-cl2-arms.R "$out_dir" ../results/373f6fd65cac \
  ../plasmode-results ../comparators-results
