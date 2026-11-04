#!/bin/bash
# Summaries of the fastFGEE comparison (one core, on LRZ: the cardiology task files never
# leave it). Submit from $PROJECT_DIR/code-fgee, with B/, P/ and, if present, Pc/ under
# <fgee results>:
#   sbatch competitors/fgee/analyse.slurm <fgee results> <out dir> <study results> <plasmode results>
# Writes <out dir>/{B,P}-*.csv and <out dir>/cardio/*.csv. Experiment A is summarised the
# same way from its own results directory (analyse.R and report.R handle all experiments).
#SBATCH -J fgee-analyse
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12000mb
#SBATCH --time=06:00:00
#SBATCH --export=NONE

set -e
module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS="$HOME/R/restart-lib"

results="$(realpath "$1")"
out="$2"
study="$3"
plasmode="$4"
mkdir -p "$out"

Rscript competitors/fgee/analyse.R "$results" "$out"
Rscript competitors/fgee/crit-summary.R "$results/B" "$out/B-crit.csv"
Rscript competitors/fgee/crit-summary.R "$results/P" "$out/P-crit.csv"
Rscript competitors/fgee/pair-ncv.R "$study" "$results/B" "$out/B-paired-ncv.csv"
Rscript competitors/fgee/paired-summary.R "$out/B-paired-ncv.csv" \
  "$out/B-paired-ncv-summary.csv" "$out/B-summary.csv"
Rscript competitors/fgee/pair-ncv.R "$plasmode" "$results/P" "$out/P-paired-ncv.csv" P
Rscript competitors/fgee/paired-summary.R "$out/P-paired-ncv.csv" \
  "$out/P-paired-ncv-summary.csv" "$out/P-summary.csv"

if [ -d "$results/Pc" ]; then
  # analyse.R expects the experiment's name as directory; Pc runs experiment P
  cardio_in="$out/cardio-in"
  mkdir -p "$cardio_in" "$out/cardio"
  ln -sfn "$results/Pc" "$cardio_in/P"
  Rscript competitors/fgee/analyse.R "$cardio_in" "$out/cardio"
  rm "$cardio_in/P" # no link to the task files inside the output directory
  rmdir "$cardio_in"
  Rscript competitors/fgee/crit-summary.R "$results/Pc" "$out/cardio/P-crit.csv"
  Rscript competitors/fgee/pair-ncv.R "$plasmode" "$results/Pc" "$out/cardio/Pc-paired-ncv.csv" P
  Rscript competitors/fgee/paired-summary.R "$out/cardio/Pc-paired-ncv.csv" \
    "$out/cardio/Pc-paired-ncv-summary.csv" "$out/cardio/P-summary.csv"
fi
