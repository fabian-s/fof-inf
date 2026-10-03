#!/bin/bash
# Bias-corrected NCV workers (claim tasks via lock dirs). From $PROJECT_DIR/code-rbc:
#   sbatch --array=1-40%30 --nice=5000000 ncv-rbc/rbc.slurm
#SBATCH -J ncv-rbc2
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=08:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

Rscript ncv-rbc/worker-rbc.R ../rbc-fresh-results 7,8,9,10,11,12,19,20,21,22,23,24,31,32,33,34,35,36,64,65,66,67,68,69,70,71,72,73,74,75,76,77 201 400 with_reml
