#!/bin/bash
# Misspecified-truth workers: each array task is one worker that claims tasks dynamically
# (lock dirs), so the array size is just the number of workers (<= 64).
# Submit from $PROJECT_DIR/code-misspec (frozen rsync of R/, scripts/, data/, misspec/):
#   sbatch --array=1-18 --time=01:30:00 misspec/misspec.slurm <ncv_xlarge_max_rep> 1 1.2   # pilot
#   sbatch --array=1-64 --time=08:00:00 misspec/misspec.slurm 100 200 7.2                  # run (as submitted)
#SBATCH -J misspec
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=23:30:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

Rscript misspec/run-misspec.R ../misspec-results "$1" "$2" "$3"
