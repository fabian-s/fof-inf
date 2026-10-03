#!/bin/bash
# Interval scores / matched coverage for the bias-correction and combination runs (paired replicates).
#SBATCH -J combo2-an
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12000mb
#SBATCH --time=04:00:00
#SBATCH --export=NONE
module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS="$HOME/R/restart-lib"
S=../relerr-summaries
Rscript ncv-rbc/analyse-rbc.R ../combo2-syn-results $S/combo2-syn-cells.csv
Rscript ncv-rbc/analyse-rbc.R ../combo2-ms-results $S/combo2-ms-cells.csv
Rscript ncv-rbc/analyse-rbc.R ../combo2-pm-results $S/combo2-pm-cells.csv
