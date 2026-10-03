#!/bin/bash
# Summaries of the plasmode jackfix run (aggregates only leave LRZ). From $PROJECT_DIR/code-jackfix:
#   sbatch ncv-jack-recheck/jackfix-pm-an.slurm
#SBATCH -J jackfix-pm-an
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000mb
#SBATCH --time=04:00:00
#SBATCH --export=NONE
module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export R_LIBS="$HOME/R/restart-lib"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
export PLASMODE_CARDIO=deid
Rscript ncv-jack-recheck/analyse-jackfix.R ../jackfix-pm-results ../jackfix-summaries/jackfix-pm plasmode
