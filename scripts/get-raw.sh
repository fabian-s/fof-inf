#!/bin/bash
# Download task records from the Zenodo deposit into raw/ and check them. Run from the repository root:
#   bash scripts/get-raw.sh results r4-results study.rds     # the named parts
#   bash scripts/get-raw.sh list                             # names of all parts
#   bash scripts/get-raw.sh all                              # everything; see README.md for the disk space
# Each archive unpacks to raw/<name>/ and is removed after its checksums pass (KEEP_ARCHIVES=1 keeps it).
set -euo pipefail

# https://doi.org/10.5281/zenodo.23916388
ZENODO_RECORD="${ZENODO_RECORD:-23916388}"
BASE_URL="${ZENODO_BASE_URL:-https://zenodo.org/records/$ZENODO_RECORD/files}"
PARTS=(
  results r4-results plasmode-results ablation-results satt-results satt-r4-syn satt-r4-pm
  misspec-results comparators-results cl2arms-results fgee-results ffpc-results lc-results
  lc-plasmode-results lc-ncv-results lc-ncvg2-results lc-ncvp-results jackfix-results
  jackfix-pm-results rbc-fresh-results rbc-pm2-results fixlam-pm-results combo-syn-results
  combo-pm-results combo2-syn-results combo2-pm-results combo2-ms-results
)

fetch() {
  curl -fL --retry 5 --retry-delay 10 -o "raw/$1" "$BASE_URL/$1?download=1"
}

get_part() {
  local name=$1
  if [ "$name" = "study.rds" ]; then
    fetch study.rds
    return
  fi
  case " ${PARTS[*]} " in
    *" $name "*) ;;
    *)
      echo "unknown part: $name (see: bash scripts/get-raw.sh list)" >&2
      exit 1
      ;;
  esac
  fetch "$name.tar.xz"
  fetch "$name.sha256"
  tar -xJf "raw/$name.tar.xz" -C raw
  (cd "raw/$name" && sha256sum --quiet -c "../$name.sha256")
  echo "$name: $(wc -l < "raw/$name.sha256") records checked"
  rm -f "raw/$name.sha256"
  [ -n "${KEEP_ARCHIVES:-}" ] || rm -f "raw/$name.tar.xz"
}

[ $# -gt 0 ] || { sed -n 2,6p "$0"; exit 1; }
if [ "$1" = "list" ]; then
  printf '%s\n' "${PARTS[@]}" study.rds
  exit 0
fi
[ -f DESIGN.md ] || { echo "run from the repository root" >&2; exit 1; }
mkdir -p raw
if [ "$1" = "all" ]; then
  set -- "${PARTS[@]}" study.rds
fi
for part in "$@"; do
  get_part "$part"
done
