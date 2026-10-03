#!/bin/bash
# Plasmode application on the public running data (LRZ): 40 workers over its 1,600 tasks.
# Submit from $PROJECT_DIR/code:  sbatch --array=0-39 slurm/plasmode-running.slurm
# Only public data may be run on LRZ (the cardiology plasmode runs locally).
#SBATCH -J plasmode-running
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

Rscript plasmode/run-plasmode.R running "$SLURM_ARRAY_TASK_ID" 40
