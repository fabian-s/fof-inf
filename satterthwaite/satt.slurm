#!/bin/bash
# Satterthwaite df workers: each array task is one worker that claims tasks dynamically
# (lock dirs), so the array size is just the number of workers (<= 64).
# Submit from $PROJECT_DIR/code-satt (frozen copy of the repository):
#   sbatch --array=1-50 satterthwaite/satt.slurm                 # production
#   sbatch --array=1 satterthwaite/satt.slurm ../satt-pilot 18   # pilot (rep 1 of every cell)
#SBATCH -J satt
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=04:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

Rscript satterthwaite/run-satt.R "${1:-../satt-results}" ../results/373f6fd65cac ${2:-}
