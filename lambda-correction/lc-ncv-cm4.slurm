#!/bin/bash
# Lambda-correction plasmode workers on ONE cm4_tiny node (used when the serial cluster's
# per-user CPU limit is taken by other arrays). N workers claim tasks from the same lock
# directory as the serial array (lc-plasmode.slurm), so both can run at once.
# Submit from $PROJECT_DIR/code-lc-ncv:
#   sbatch lambda-correction/lc-plasmode-cm4.slurm [n_workers] [out_dir]
#SBATCH -J lc-ncv-cm4
#SBATCH -o ../logs/%x.%j.out
#SBATCH -e ../logs/%x.%j.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=60000mb
#SBATCH --time=06:00:00
#SBATCH --export=NONE

module load slurm_setup
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export REFUND_SHA="$(cat "$RESTART_LIB/REFUND_SHA")"

N="${1:-28}"
OUT="${2:-../lc-ncv-results}"
logs="../logs/lc-ncv-cm4-$SLURM_JOB_ID"
mkdir -p "$logs" "$OUT"
for w in $(seq 1 "$N"); do
  Rscript lambda-correction/run-lc-ncv.R "$OUT" ../results/373f6fd65cac Inf 1:200 \
    > "$logs/worker$w.log" 2>&1 &
done
wait
echo "cm4 workers done"
