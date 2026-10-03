#!/bin/bash
# Interval scores / matched coverage for the bias-correction and combination runs (paired replicates).
#SBATCH -J score-an
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
Rscript ncv-rbc/analyse-rbc.R ../rbc-fresh-results $S/score-fresh-cells.csv 201 300
Rscript ncv-rbc/analyse-rbc.R ../combo-syn-results $S/score-combo-syn-cells.csv 201 300
Rscript ncv-rbc/analyse-rbc.R ../rbc-pm2-results $S/score-pm2-cells.csv 1 50
Rscript ncv-rbc/analyse-rbc.R ../combo-pm-results $S/score-combo-pm-cells.csv 1 50
