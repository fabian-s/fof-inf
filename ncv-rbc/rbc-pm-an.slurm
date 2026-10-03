#!/bin/bash
#SBATCH -J rbc-pm-an
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12000mb
#SBATCH --time=02:00:00
#SBATCH --export=NONE
module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS="$HOME/R/restart-lib"
Rscript ncv-rbc/analyse-rbc.R ../rbc-pm2-results ../relerr-summaries/rbc-plasmode-cells.csv
