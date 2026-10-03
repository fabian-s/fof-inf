#!/bin/bash
# Corrected NCV covariance on the plasmode cells of the six independent-curve datasets
# (all subjects, curve flips, attached errors). From $PROJECT_DIR/code-jackfix:
#   sbatch --array=1-40 ncv-jack-recheck/jackfix-pm.slurm
#SBATCH -J jackfix-pm
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000mb
#SBATCH --time=08:00:00
#SBATCH --export=NONE
module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
export PLASMODE_CARDIO=deid
Rscript ncv-jack-recheck/run-jackfix.R ../jackfix-pm-results plasmode ecg_gls4b,af,running,dti,gait_vc,ecg_roahd 1 100
