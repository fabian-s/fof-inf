#!/bin/bash
# Pack the allowlisted records of one results folder for the deposit: strip each .rds's gzip layer
# (readRDS reads the result unchanged), then tar + xz, which compresses across records.
# Usage: zenodo-pack.sh <results_dir> <keep.txt> <out_dir> <threads> <tmp_dir>
set -euo pipefail
results_dir=${1%/}
keep=$2
out_dir=$3
threads=$4
tmp_dir=$5
name=$(basename "$results_dir")
work="$tmp_dir/$name"
rm -rf "$work"
mkdir -p "$work/$name" "$out_dir"
unpack_one() {
  mkdir -p "$(dirname "$3/$1")"
  if [ "$(head -c2 "$2/$1" | od -An -tx1 | tr -d ' ')" = "1f8b" ]; then
    gzip -dc "$2/$1" > "$3/$1"
  else
    cp "$2/$1" "$3/$1"
  fi
}
export -f unpack_one
xargs -a "$keep" -d '\n' -P "$threads" -I{} bash -c 'unpack_one "$1" "$2" "$3"' _ {} "$results_dir" "$work/$name"
(cd "$work/$name" && find . -type f -printf '%P\n' | sort | xargs -d '\n' sha256sum) > "$out_dir/$name.sha256"
tar -C "$work" --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-10-03 00:00:00 UTC' -cf - "$name" |
  xz -9 -T"$threads" > "$out_dir/$name.tar.xz.partial"
mv "$out_dir/$name.tar.xz.partial" "$out_dir/$name.tar.xz"
rm -rf "$work"
