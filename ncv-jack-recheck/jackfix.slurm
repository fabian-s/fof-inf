#!/bin/bash
# Corrected NCV covariance on the study's replicates (cells 12/24/36, reps 1-100).
# From $PROJECT_DIR/code-jackfix:
#   sbatch --array=1-30 ncv-jack-recheck/jackfix.slurm
#SBATCH -J jackfix
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000mb
#SBATCH --time=04:00:00
#SBATCH --export=NONE
module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
Rscript ncv-jack-recheck/run-jackfix.R ../jackfix-results syn 12,24,36 1 100
