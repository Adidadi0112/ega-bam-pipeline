#!/usr/bin/env bash
# Download 1000 Genomes JPT exome BAMs and run Wave 2 GVCF calling.
#
# Default: the 70 modelled JPT samples (paper/config/wave2_jpt_urls.txt),
# HTTPS, resume, companion BAI, skip finished GVCFs, keep BAMs.
#
# This does not download EGA UC alignments. Point variant_calling.sh -d at
# the local EGA BAM/CRAM directory after those files exist.
#
# Optional env:
#   LINKS_FILE ALL_JPT BAM_DIR DELETE_BAM_AFTER_GVCF DOWNLOAD_ONLY DRY_RUN
set -euo pipefail

if [[ "${TRACE:-0}" == "1" ]]; then
  set -x
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALL_JPT="${ALL_JPT:-0}"
DOWNLOAD_ONLY="${DOWNLOAD_ONLY:-0}"
DELETE_BAM_AFTER_GVCF="${DELETE_BAM_AFTER_GVCF:-0}"
DRY_RUN="${DRY_RUN:-0}"
BAM_DIR="${BAM_DIR:-$SCRIPT_DIR/bams/jpt}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/gvcf_results}"

if [[ "$ALL_JPT" == "1" ]]; then
  LINKS_FILE="${LINKS_FILE:-$SCRIPT_DIR/jpt_mapped_only.txt}"
else
  LINKS_FILE="${LINKS_FILE:-$SCRIPT_DIR/paper/config/wave2_jpt_urls.txt}"
fi

if [[ ! -f "$LINKS_FILE" ]]; then
  echo "ERROR: URL list not found: $LINKS_FILE" >&2
  exit 1
fi

mkdir -p "$BAM_DIR" "$RESULTS_DIR/logs"

to_https() {
  printf '%s' "$1" | sed 's|^ftp://ftp.1000genomes.ebi.ac.uk|https://ftp.1000genomes.ebi.ac.uk|'
}

download_file() {
  local url=$1
  local dest=$2
  if [[ -s "$dest" ]]; then
    echo "Already present: $dest"
    return 0
  fi
  echo "Downloading $url"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -c -x 16 -s 16 -k 1M --auto-file-renaming=false -d "$(dirname "$dest")" -o "$(basename "$dest")" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 10 -C - -o "$dest" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -c -O "$dest" "$url"
  else
    echo "ERROR: need aria2c, curl, or wget" >&2
    return 1
  fi
}

n_urls=$(grep -cve '^[[:space:]]*$' "$LINKS_FILE" || true)
echo "JPT URL list: $LINKS_FILE ($n_urls files)"
echo "BAM directory: $BAM_DIR"
echo "GVCF directory: $RESULTS_DIR"

while read -r url || [[ -n "${url:-}" ]]; do
  url=${url%$'\r'}
  [[ -z "${url:-}" || "$url" =~ ^# ]] && continue
  url=$(to_https "$url")
  filename=$(basename "$url")
  sample_id=${filename%.bam}
  bam_path="$BAM_DIR/$filename"
  bai_path="${bam_path}.bai"
  gvcf_path="$RESULTS_DIR/${sample_id}.g.vcf.gz"

  echo "################################################"
  echo "Sample: $sample_id"
  echo "################################################"

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "DRY_RUN $url -> $bam_path"
    continue
  fi

  if [[ "${SKIP_EXISTING:-1}" == "1" && -s "$gvcf_path" && ( -f "${gvcf_path}.tbi" || -f "${gvcf_path}.idx" ) ]]; then
    echo "GVCF exists, skipping download and call: $gvcf_path"
    continue
  fi

  download_file "$url" "$bam_path"
  if ! download_file "${url}.bai" "$bai_path"; then
    echo "WARNING: BAI download failed; indexing with samtools"
    samtools index "$bam_path"
  fi

  if [[ "$DOWNLOAD_ONLY" == "1" ]]; then
    continue
  fi

  RESULTS_DIR="$RESULTS_DIR" bash "$SCRIPT_DIR/variant_calling.sh" \
    -f "$bam_path" \
    -s "$sample_id" \
    -o "$gvcf_path"

  if [[ "$DELETE_BAM_AFTER_GVCF" == "1" ]]; then
    echo "DELETE_BAM_AFTER_GVCF=1: removing $bam_path"
    rm -f "$bam_path" "$bai_path"
  fi
done < "$LINKS_FILE"

echo "Done. GVCFs in $RESULTS_DIR"
echo "Next: call UC BAMs with variant_calling.sh -d, then joint_genotype.sh (not merge_vcf.sh)."
