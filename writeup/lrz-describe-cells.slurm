#!/bin/bash
# Realized per-cell design effects of the synthetic study (analysis/describe-cells.R), run on
# LRZ from a writable copy of the frozen synthetic code ($PROJECT_DIR/code-analysis is
# read-only and must not be modified). Environment mirrors slurm/array.slurm.
# Setup (login node, file ops only):
#   ssh lrz "mkdir -p $PROJECT_DIR/writeup-describe/logs && \
#            cp -r $PROJECT_DIR/code-analysis/. $PROJECT_DIR/writeup-describe/ && \
#            chmod -R u+w $PROJECT_DIR/writeup-describe"
#   scp writeup/lrz-describe-cells.slurm lrz:$PROJECT_DIR/writeup-describe/
#   ssh lrz "cd $PROJECT_DIR/writeup-describe && sbatch lrz-describe-cells.slurm"
# Fetch: scp lrz:$PROJECT_DIR/writeup-describe/synthetic-design-effects.csv writeup/tables/
#SBATCH -J describe-cells
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000mb
#SBATCH --time=04:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export MGCV_EXPECT=1.9.5

Rscript analysis/describe-cells.R synthetic-design-effects.csv 20
