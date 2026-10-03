#!/bin/bash
# Covariance ablation workers: each array task is one worker that claims ablation tasks
# dynamically (lock dirs), so the array size is just the number of workers (<= 64).
# Submit from $PROJECT_DIR/code-ablation (frozen copy of the repository):
#   sbatch --array=1-60 ablation/ablation.slurm
#SBATCH -J ablation
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=05:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

Rscript ablation/run-ablation.R ../ablation-results ../results/373f6fd65cac
