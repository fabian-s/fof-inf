#!/bin/bash
# Plasmode extension smoke run on all public datasets in ONE cm4_tiny job (used when the
# serial cluster's 200-job submit limit is taken by a production array): per app, the
# truth caches are built ("prepare"), then 4 workers run R = 2 replicates of every cell.
# Caches and results go to $PROJECT_DIR/plasmode-smoke/. Public data only.
# Submit from $PROJECT_DIR/code-plasmode:  sbatch slurm/plasmode-smoke-cm4.slurm
# Timing caveat: cm4 cores are slower per thread than the serial nodes the production
# runs use; calibrate with the running (REML, NCV) cells against summaries/plasmode/fits.csv.
#SBATCH -J plasmode-smoke
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
PROJECT_DIR="${PROJECT_DIR:-$(dirname "$PWD")}" # jobs run from a copy of the repository inside PROJECT_DIR
export PATH="$HOME/envs/r461/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RESTART_LIB="$HOME/R/restart-lib"
export PLASMODE_CACHE="$PROJECT_DIR/plasmode-smoke/cache"
export PLASMODE_RESULTS="$PROJECT_DIR/plasmode-smoke/results"
logs="$PROJECT_DIR/plasmode-smoke/logs"
mkdir -p "$PLASMODE_CACHE" "$PLASMODE_RESULTS" "$logs"

run_app() {
  local app="$1"
  /usr/bin/time -v Rscript plasmode/run-plasmode.R "$app" prepare 4 2 \
    > "$logs/$app-prepare.log" 2>&1 || return 1
  for worker in 0 1 2 3; do
    /usr/bin/time -v Rscript plasmode/run-plasmode.R "$app" "$worker" 4 2 \
      > "$logs/$app-worker$worker.log" 2>&1 &
  done
  wait
}

for app in running dti gait_vc ecg_roahd ocean aemet adelaide; do
  run_app "$app" &
done
wait
echo "smoke done"
