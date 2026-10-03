#!/bin/bash
# CL2-arms summaries (one core, on LRZ: the cardiology task files never leave it).
# Submit from $PROJECT_DIR/code-cl2arms:
#   sbatch comparators/cl2-arms-analyse.slurm ../cl2arms-results ../cl2arms-summaries
#SBATCH -J cl2arms-analyse
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

Rscript comparators/analyse-cl2-arms.R "$1" ../comparators-results \
  ../results/373f6fd65cac ../plasmode-results "$2"
