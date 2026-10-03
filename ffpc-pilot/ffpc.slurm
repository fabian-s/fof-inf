#!/bin/bash
# ffpc() workers (2026-10-03): each array task is one worker that claims tasks dynamically
# (lock dirs). Submit from the frozen code copy (of the repository). Arguments:
# <out_dir> <cells (R expression)> <rules> <methods> <rep_from> <rep_to>
#   sbatch --array=1-40 ffpc-pilot/ffpc.slurm ../ffpc-results "5,11,12,35,107,126,129" p095,p090,k3 REML,NCV 1 200
#SBATCH -J ffpc
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000mb
#SBATCH --time=24:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

# stored study records (ff REML/NCV) of the same cells and replicates (and subdirectories)
export FFPC_REF_DIRS="$PROJECT_DIR/results:$PROJECT_DIR/r4-results"

Rscript ffpc-pilot/run-ffpc.R "$1" "$2" "$3" "$4" "$5" "$6"
