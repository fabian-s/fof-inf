#!/bin/bash
# Upload a prepared deposit folder to an existing Zenodo draft. Skips files already there with the same size
# and md5, so it can be rerun after a failure. Never publishes: the record owner does that on the website.
#   bash scripts/zenodo-upload.sh <deposit_dir> <deposition_id> [token_file]
# token_file (default ~/.zenodo-token) holds a personal access token with deposit scope; keep it mode 600.
set -euo pipefail
dir=${1%/}
id=$2
token_file=${3:-$HOME/.zenodo-token}
api="https://zenodo.org/api/deposit/depositions/$id"

[ -r "$token_file" ] || { echo "no token file: $token_file" >&2; exit 1; }
[ "$(stat -c %a "$token_file")" = "600" ] || { echo "$token_file must be mode 600" >&2; exit 1; }
# The token goes through a private header file, not the command line, which other users may see.
header_file=$(mktemp)
trap 'rm -f "$header_file"' EXIT
chmod 600 "$header_file"
printf 'Authorization: Bearer %s\n' "$(tr -d '[:space:]' < "$token_file")" > "$header_file"
auth=(-H "@$header_file")

meta=$(curl -sS --fail-with-body --retry 5 --retry-delay 20 "${auth[@]}" "$api")
state=$(jq -r '.state' <<< "$meta")
bucket=$(jq -r '.links.bucket // empty' <<< "$meta")
[ -n "$bucket" ] || { echo "draft $id has no upload bucket (state: $state)" >&2; exit 1; }
echo "draft $id, state $state"

existing=$(curl -sS --fail-with-body --retry 5 --retry-delay 20 "${auth[@]}" "$bucket" | jq -r '.contents[]? | "\(.key) \(.size) \(.checksum)"')

for path in "$dir"/*; do
  name=$(basename "$path")
  size=$(stat -c %s "$path")
  md5="md5:$(md5sum "$path" | cut -d' ' -f1)"
  if grep -qxF "$name $size $md5" <<< "$existing"; then
    echo "already there: $name"
    continue
  fi
  echo "uploading: $name ($size bytes)"
  reply=$(curl -sS --fail-with-body --retry 3 --retry-delay 30 "${auth[@]}" \
    -H "Content-Type: application/octet-stream" --upload-file "$path" "$bucket/$name")
  got=$(jq -r '.checksum' <<< "$reply")
  [ "$got" = "$md5" ] || { echo "checksum mismatch for $name: $got vs $md5" >&2; exit 1; }
done

# Final check against what Zenodo reports.
final=$(curl -sS --fail-with-body --retry 5 --retry-delay 20 "${auth[@]}" "$bucket" | jq -r '.contents[]? | "\(.key) \(.size) \(.checksum)"')
missing=0
for path in "$dir"/*; do
  name=$(basename "$path")
  line="$name $(stat -c %s "$path") md5:$(md5sum "$path" | cut -d' ' -f1)"
  grep -qxF "$line" <<< "$final" || { echo "missing or different: $name" >&2; missing=$((missing + 1)); }
done
[ "$missing" -eq 0 ] || exit 1
echo "all $(ls "$dir" | wc -l) files on Zenodo match; publish the draft on the website"
