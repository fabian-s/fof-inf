#!/bin/bash
# Analysis of the lambda-correction plasmode run on LRZ (task files of the cardiology apps
# stay there; only the aggregate summaries are copied back). Submit from
# $PROJECT_DIR/code-lc-ncvp: sbatch lambda-correction/lc-plasmode-analyse.slurm
#SBATCH -J lc-ncvp-analyse
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000mb
#SBATCH --time=02:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
Rscript lambda-correction/analyse-lc-ncv-plasmode.R ../lc-ncvp-results ../lc-ncvp-summaries
