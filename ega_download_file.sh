#!/usr/bin/env bash
# Fetch one EGA file (EGAF*) with pyega3. Controlled-access; not a public URL.
#
# Usage:
#   ega_download_file.sh EGAF000... OUTPUT_DIR [EXPECTED.bam]
#   ega_download_file.sh -l IDS.txt OUTPUT_DIR
#
# Credentials: PYEGA3_CREDENTIALS (default: ~/programs/pyega3_credentials.json)
# Connections: PYEGA3_CONNECTIONS (default: 8)
#
# Prints the BAM path to stdout. Logs go to stderr.
set -euo pipefail

PYEGA3_CREDENTIALS="${PYEGA3_CREDENTIALS:-$HOME/programs/pyega3_credentials.json}"
PYEGA3_CONNECTIONS="${PYEGA3_CONNECTIONS:-8}"

usage() {
  echo "Usage: $0 EGAF000... OUTPUT_DIR [EXPECTED.bam]" >&2
  echo "       $0 -l IDS.txt OUTPUT_DIR" >&2
  exit 1
}

locate_bam() {
  local dest=$1
  local expected=${2:-}
  if [[ -n "$expected" && -s "$dest/$expected" ]]; then
    printf '%s\n' "$dest/$expected"
    return 0
  fi
  if [[ -n "$expected" ]]; then
    local found
    found=$(find "$dest" -type f -name "$expected" ! -name '*.cip' ! -name '*.partial' 2>/dev/null | awk 'NR==1{print; exit}')
    if [[ -n "${found:-}" && -s "$found" ]]; then
      printf '%s\n' "$found"
      return 0
    fi
  fi
  local any
  any=$(find "$dest" -type f -name '*.bam' ! -name '*.bai' ! -name '*.cip' 2>/dev/null | awk 'NR==1{print; exit}')
  if [[ -n "${any:-}" && -s "$any" ]]; then
    printf '%s\n' "$any"
    return 0
  fi
  return 1
}

fetch_one() {
  local egaf=$1
  local dest=$2
  local expected=${3:-}

  if ! command -v pyega3 >/dev/null 2>&1; then
    echo "ERROR: pyega3 not found on PATH" >&2
    return 1
  fi
  if [[ ! -f "$PYEGA3_CREDENTIALS" ]]; then
    echo "ERROR: pyega3 credentials file not found: $PYEGA3_CREDENTIALS" >&2
    echo "Set PYEGA3_CREDENTIALS to your EGA JSON credentials." >&2
    return 1
  fi

  mkdir -p "$dest"
  if bam=$(locate_bam "$dest" "$expected"); then
    echo "Already present: $bam" >&2
    printf '%s\n' "$bam"
    return 0
  fi

  echo "pyega3 fetch $egaf -> $dest" >&2
  pyega3 \
    -c "$PYEGA3_CONNECTIONS" \
    -cf "$PYEGA3_CREDENTIALS" \
    fetch "$egaf" \
    --output-dir "$dest" >&2

  if ! bam=$(locate_bam "$dest" "$expected"); then
    echo "ERROR: pyega3 finished but no BAM found under $dest" >&2
    return 1
  fi
  printf '%s\n' "$bam"
}

if [[ $# -lt 1 ]]; then
  usage
fi

if [[ "$1" == "-l" ]]; then
  [[ $# -eq 3 ]] || usage
  list=$2
  dest=$3
  while read -r egaf || [[ -n "${egaf:-}" ]]; do
    egaf=${egaf%$'\r'}
    egaf=$(printf '%s' "$egaf" | tr -d '[:space:],')
    [[ -z "$egaf" || "$egaf" =~ ^# ]] && continue
    fetch_one "$egaf" "$dest"
  done < "$list"
  exit 0
fi

[[ $# -ge 2 ]] || usage
fetch_one "$1" "$2" "${3:-}"
