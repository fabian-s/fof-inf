#!/bin/bash
# Plasmode study on one dataset (LRZ): the public sets, and the cardiology apps ecg_gls4b
# and af from their de-identified input only (PLASMODE_CARDIO=deid; the files under
# $PLASMODE_CACHE/deid/<app>/). The original ECG app ecg_gls4 is refused.
# Code: $PROJECT_DIR/code-plasmode (a copy of the repository). Caches and results:
# $PROJECT_DIR/plasmode-cache and $PROJECT_DIR/plasmode-results/<app>.
# Arguments: <app> <n_workers> [max_rep] [tag]. A tag (e.g. smoke) puts caches and results
# under $PROJECT_DIR/plasmode-<tag>/ instead. Without --array the job builds the app's
# truth caches ("prepare"); with --array=0-<n_workers - 1> each task runs one worker's share.
#   prep=$(sbatch --parsable slurm/plasmode-public.slurm gait_vc 40)
#   sbatch --dependency=afterok:${prep%%;*} --array=0-39 slurm/plasmode-public.slurm gait_vc 40
#SBATCH -J plasmode-public
#SBATCH -o ../logs/%x.%A_%a.out
#SBATCH -e ../logs/%x.%A_%a.err
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500mb
#SBATCH --time=08:00:00
#SBATCH --export=NONE

module load slurm_setup
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"

app="$1"
n_workers="$2"
max_rep="${3:-}"
tag="${4:-}"
if [ -n "$tag" ]; then
  export PLASMODE_CACHE="$PROJECT_DIR/plasmode-$tag/cache"
  export PLASMODE_RESULTS="$PROJECT_DIR/plasmode-$tag/results"
else
  export PLASMODE_CACHE="$PROJECT_DIR/plasmode-cache"
  export PLASMODE_RESULTS="$PROJECT_DIR/plasmode-results"
fi
case "$app" in
  ecg_gls4)
    echo "refusing the original patient-data app $app on the cluster" >&2
    exit 1
    ;;
  ecg_gls4b | af)
    export PLASMODE_CARDIO=deid
    ;;
esac
Rscript plasmode/run-plasmode.R "$app" "${SLURM_ARRAY_TASK_ID:-prepare}" "$n_workers" $max_rep
