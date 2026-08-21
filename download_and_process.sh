#!/usr/bin/env bash
# Wave 2: download one alignment at a time, call a GVCF, delete the BAM.
#
# COHORT=both (default) | uc | jpt
#   uc  — 74 modelled germline WES BAMs from EGA EGAD00001005237 via pyega3
#   jpt — 70 modelled 1000 Genomes JPT exome BAMs via HTTPS
#
# There are no public HTTP URLs for the UC files. pyega3 fetches EGAF accessions
# listed in paper/config/wave2_uc_ega_files.tsv. Do not fetch the whole dataset
# (827 files, including crypts, tumours, and RNA-seq).
#
# Optional env: COHORT ALL_JPT LINKS_FILE UC_MAP BAM_DIR UC_BAM_DIR JPT_BAM_DIR
#   RESULTS_DIR DELETE_BAM_AFTER_GVCF DOWNLOAD_ONLY DRY_RUN PYEGA3_CREDENTIALS
set -euo pipefail

if [[ "${TRACE:-0}" == "1" ]]; then
  set -x
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT="${COHORT:-both}"
ALL_JPT="${ALL_JPT:-0}"
DOWNLOAD_ONLY="${DOWNLOAD_ONLY:-0}"
DELETE_BAM_AFTER_GVCF="${DELETE_BAM_AFTER_GVCF:-1}"
DRY_RUN="${DRY_RUN:-0}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/gvcf_results}"
JPT_BAM_DIR="${JPT_BAM_DIR:-${BAM_DIR:-$SCRIPT_DIR/bams/jpt}}"
UC_BAM_DIR="${UC_BAM_DIR:-$SCRIPT_DIR/bams/uc}"
UC_MAP="${UC_MAP:-$SCRIPT_DIR/paper/config/wave2_uc_ega_files.tsv}"
EGA_FETCH="$SCRIPT_DIR/ega_download_file.sh"

if [[ "$ALL_JPT" == "1" ]]; then
  LINKS_FILE="${LINKS_FILE:-$SCRIPT_DIR/jpt_mapped_only.txt}"
else
  LINKS_FILE="${LINKS_FILE:-$SCRIPT_DIR/paper/config/wave2_jpt_urls.txt}"
fi

mkdir -p "$RESULTS_DIR/logs"

to_https() {
  printf '%s' "$1" | sed 's|^ftp://ftp.1000genomes.ebi.ac.uk|https://ftp.1000genomes.ebi.ac.uk|'
}

download_http() {
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

gvcf_ready() {
  local gvcf=$1
  local sample=$2
  if [[ ! -s "$gvcf" ]]; then
    return 1
  fi
  if [[ ! -f "${gvcf}.tbi" && ! -f "${gvcf}.idx" ]]; then
    return 1
  fi
  if command -v bcftools >/dev/null 2>&1; then
    bcftools query -l "$gvcf" 2>/dev/null | grep -Fxq "$sample"
  fi
}

remove_alignment() {
  local bam_path=$1
  local extra_dir=${2:-}
  rm -f "$bam_path" "${bam_path}.bai" "${bam_path}.crai" \
    "${bam_path%.bam}.bai" "${bam_path}.aria2"
  if [[ -n "$extra_dir" && -d "$extra_dir" ]]; then
    rm -rf "$extra_dir"
  fi
  echo "Removed alignment to free disk: $bam_path"
}

call_gvcf() {
  local bam_path=$1
  local sample_id=$2
  local gvcf_path=$3
  RESULTS_DIR="$RESULTS_DIR" bash "$SCRIPT_DIR/variant_calling.sh" \
    -f "$bam_path" \
    -s "$sample_id" \
    -o "$gvcf_path"
  if ! gvcf_ready "$gvcf_path" "$sample_id"; then
    echo "ERROR: GVCF missing or incomplete for $sample_id; keeping $bam_path for retry" >&2
    exit 1
  fi
}

process_jpt() {
  if [[ ! -f "$LINKS_FILE" ]]; then
    echo "ERROR: JPT URL list not found: $LINKS_FILE" >&2
    exit 1
  fi
  mkdir -p "$JPT_BAM_DIR"
  local n_urls
  n_urls=$(grep -cve '^[[:space:]]*$' "$LINKS_FILE" || true)
  echo "JPT URL list: $LINKS_FILE ($n_urls files)"
  echo "JPT BAM directory: $JPT_BAM_DIR"

  while read -r url || [[ -n "${url:-}" ]]; do
    url=${url%$'\r'}
    [[ -z "${url:-}" || "$url" =~ ^# ]] && continue
    url=$(to_https "$url")
    local filename sample_id bam_path bai_path gvcf_path
    filename=$(basename "$url")
    sample_id=${filename%.bam}
    bam_path="$JPT_BAM_DIR/$filename"
    bai_path="${bam_path}.bai"
    gvcf_path="$RESULTS_DIR/${sample_id}.g.vcf.gz"

    echo "################################################"
    echo "JPT sample: $sample_id"
    echo "################################################"

    if [[ "$DRY_RUN" == "1" ]]; then
      echo "DRY_RUN $url -> $bam_path"
      continue
    fi

    if [[ "${SKIP_EXISTING:-1}" == "1" ]] && gvcf_ready "$gvcf_path" "$sample_id"; then
      echo "GVCF exists, skipping download and call: $gvcf_path"
      if [[ "$DELETE_BAM_AFTER_GVCF" == "1" && -e "$bam_path" ]]; then
        remove_alignment "$bam_path"
      fi
      continue
    fi

    download_http "$url" "$bam_path"
    if ! download_http "${url}.bai" "$bai_path"; then
      echo "WARNING: BAI download failed; indexing with samtools"
      samtools index "$bam_path"
    fi

    if [[ "$DOWNLOAD_ONLY" == "1" ]]; then
      continue
    fi

    call_gvcf "$bam_path" "$sample_id" "$gvcf_path"
    if [[ "$DELETE_BAM_AFTER_GVCF" == "1" ]]; then
      remove_alignment "$bam_path"
    fi
  done < "$LINKS_FILE"
}

process_uc() {
  if [[ ! -f "$UC_MAP" ]]; then
    echo "ERROR: UC EGA map not found: $UC_MAP" >&2
    exit 1
  fi
  if [[ ! -x "$EGA_FETCH" ]]; then
    echo "ERROR: $EGA_FETCH is not executable" >&2
    exit 1
  fi
  mkdir -p "$UC_BAM_DIR"
  local n_uc
  n_uc=$(awk 'NR>1 && $1!="" {c++} END{print c+0}' "$UC_MAP")
  echo "UC EGA map: $UC_MAP ($n_uc germline WXS files)"
  echo "UC BAM directory: $UC_BAM_DIR"
  echo "pyega3 credentials: ${PYEGA3_CREDENTIALS:-$HOME/programs/pyega3_credentials.json}"

  local header=1
  while IFS=$'\t' read -r sample_id egan egaf file_name || [[ -n "${sample_id:-}" ]]; do
    sample_id=${sample_id%$'\r'}
    [[ -z "${sample_id:-}" || "$sample_id" =~ ^# ]] && continue
    if [[ "$header" == "1" && "$sample_id" == "sample_id" ]]; then
      header=0
      continue
    fi
    header=0
    file_name=${file_name%$'\r'}
    local staging gvcf_path bam_path
    staging="$UC_BAM_DIR/$sample_id"
    gvcf_path="$RESULTS_DIR/${sample_id}.g.vcf.gz"

    echo "################################################"
    echo "UC sample: $sample_id"
    echo "  EGAN $egan  EGAF $egaf  file $file_name"
    echo "################################################"

    if [[ "$DRY_RUN" == "1" ]]; then
      echo "DRY_RUN pyega3 fetch $egaf -> $staging/$file_name"
      continue
    fi

    if [[ "${SKIP_EXISTING:-1}" == "1" ]] && gvcf_ready "$gvcf_path" "$sample_id"; then
      echo "GVCF exists, skipping download and call: $gvcf_path"
      if [[ "$DELETE_BAM_AFTER_GVCF" == "1" ]]; then
        rm -rf "$staging"
      fi
      continue
    fi

    bam_path=$("$EGA_FETCH" "$egaf" "$staging" "$file_name")
    if [[ ! -s "$bam_path" ]]; then
      echo "ERROR: EGA download produced no BAM for $sample_id" >&2
      exit 1
    fi
    if [[ ! -f "${bam_path}.bai" && ! -f "${bam_path%.*}.bai" ]]; then
      echo "Indexing $bam_path"
      samtools index "$bam_path"
    fi

    if [[ "$DOWNLOAD_ONLY" == "1" ]]; then
      continue
    fi

    call_gvcf "$bam_path" "$sample_id" "$gvcf_path"
    if [[ "$DELETE_BAM_AFTER_GVCF" == "1" ]]; then
      remove_alignment "$bam_path" "$staging"
    fi
  done < "$UC_MAP"
}

echo "Cohort: $COHORT"
echo "GVCF directory: $RESULTS_DIR"
echo "Delete BAM after successful GVCF: $DELETE_BAM_AFTER_GVCF"

case "$COHORT" in
  uc) process_uc ;;
  jpt) process_jpt ;;
  both)
    process_uc
    process_jpt
    ;;
  *)
    echo "ERROR: COHORT must be uc, jpt, or both (got $COHORT)" >&2
    exit 1
    ;;
esac

echo "Done. GVCFs in $RESULTS_DIR"
echo "Next: joint_genotype.sh (not merge_vcf.sh). See paper/docs/bam_gvcf_checklist.md"
