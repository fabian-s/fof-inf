#!/bin/bash
# Lambda-correction plasmode workers on ONE cm4_tiny node (used when the serial cluster's
# per-user CPU limit is taken by other arrays). N workers claim tasks from the same lock
# directory as the serial array (lc-plasmode.slurm), so both can run at once.
# Submit from $PROJECT_DIR/code-lc-ncvp:
#   sbatch lambda-correction/lc-plasmode-cm4.slurm [n_workers] [out_dir]
#SBATCH -J lc-ncvp-cm4
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=110000mb
#SBATCH --time=06:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
export PLASMODE_CARDIO=deid

N="${1:-56}"
OUT="${2:-../lc-ncvp-results}"
logs="../logs/lc-ncvp-cm4-$SLURM_JOB_ID"
mkdir -p "$logs" "$OUT"
for w in $(seq 1 "$N"); do
  Rscript lambda-correction/run-lc-ncv-plasmode.R "$OUT" ../plasmode-results Inf 1:100 \
    > "$logs/worker$w.log" 2>&1 &
done
wait
echo "cm4 workers done"
