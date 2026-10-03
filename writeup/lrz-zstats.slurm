#!/bin/bash
# Studentised-error diagnostics (writeup/lrz-zstats.R) from the stored task records on
# LRZ. Runs from a writable copy of the frozen code with writeup/lrz-zstats.R added
# (the frozen copies $PROJECT_DIR/code-analysis and code-plasmode-v1 stay untouched):
#   ssh lrz "mkdir -p $PROJECT_DIR/writeup-zstats/logs && \
#            cp -r $PROJECT_DIR/code-plasmode-v1/. $PROJECT_DIR/writeup-zstats/ && \
#            chmod -R u+w $PROJECT_DIR/writeup-zstats && mkdir -p $PROJECT_DIR/writeup-zstats/writeup"
#   scp writeup/lrz-zstats.R writeup/lrz-zstats.slurm lrz:$PROJECT_DIR/writeup-zstats/writeup/
#   ssh lrz "cd $PROJECT_DIR/writeup-zstats && sbatch writeup/lrz-zstats.slurm"
# Fetch: scp lrz:$PROJECT_DIR/writeup-zstats/zstats-*.csv writeup/tables/
# (The plasmode-v1 code copy contains the synthetic engine too; if the synthetic design
# table differs from code-analysis, run the synthetic half from a copy of code-analysis.)
#SBATCH -J zstats
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000mb
#SBATCH --time=06:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export MGCV_EXPECT=1.9.5

# synthetic: both runs of 100 replicates
Rscript writeup/lrz-zstats.R synthetic "$PROJECT_DIR/results" zstats-synthetic.csv
# plasmode: the nine applications
Rscript writeup/lrz-zstats.R plasmode "$PROJECT_DIR/plasmode-results" zstats-plasmode.csv \
  dti gait_vc ecg_roahd ocean aemet adelaide running ecg_gls4b af
