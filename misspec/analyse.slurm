#!/bin/bash
# Summaries of the misspecified-truth test (one core). Submit from $PROJECT_DIR/code-misspec:
#   sbatch misspec/analyse.slurm <ncv_xlarge_max_rep>
#SBATCH -J misspec-analyse
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12000mb
#SBATCH --time=01:30:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"

Rscript misspec/analyse-misspec.R ../misspec-results ../misspec-summaries "$1"
