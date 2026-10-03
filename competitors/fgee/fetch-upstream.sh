#!/bin/bash
# Downloads the authors' data generator (no license, so not redistributed here) at the commit used in the study.
cd "$(dirname "$0")"
curl -sfL https://raw.githubusercontent.com/gloewing/fgee_onestep/de6540003f09ac736a04e6139f2a8a6724ee13f7/simulations/functions/sim_functional_data.R \
  -o upstream-sim_functional_data.R
