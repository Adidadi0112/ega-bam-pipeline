#!/usr/bin/env bash
# Wave 2 per-sample GVCF calling: BQSR + HaplotypeCaller -ERC GVCF.
#
# The previous version emitted variant-only VCFs (no reference blocks). Those
# files cannot be jointly genotyped and were the source of the Wave 1
# missing-as-zero artifact. Do not use -ERC NONE for this paper.
#
# Required tools: gatk, samtools
# Optional env:
#   REF, KNOWN_SITES, RESULTS_DIR, INTERVALS_BED, GATK_JAVA_OPTS,
#   HC_THREADS, SKIP_EXISTING (default 1), KEEP_WORK (default 0),
#   DELETE_BAM_AFTER_GVCF (default 0 here; download_and_process.sh sets 1)
set -euo pipefail

if [[ "${TRACE:-0}" == "1" ]]; then
  set -x
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="${REF:-$SCRIPT_DIR/ref/human_g1k_v37.fasta}"
KNOWN_SITES="${KNOWN_SITES:-$SCRIPT_DIR/ref/dbsnp_138.b37.vcf}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/gvcf_results}"
INTERVALS_BED="${INTERVALS_BED:-}"
GATK_JAVA_OPTS="${GATK_JAVA_OPTS:--Xmx8g}"
HC_THREADS="${HC_THREADS:-4}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
KEEP_WORK="${KEEP_WORK:-0}"
DELETE_BAM_AFTER_GVCF="${DELETE_BAM_AFTER_GVCF:-0}"

gatk_cmd() {
  gatk --java-options "$GATK_JAVA_OPTS" "$@"
}

ensure_index() {
  local vcf=$1
  if [[ -f "${vcf}.idx" || -f "${vcf}.tbi" || -f "${vcf}.csi" ]]; then
    return 0
  fi
  gatk_cmd IndexFeatureFile -I "$vcf"
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

maybe_delete_bam() {
  local input_bam=$1
  local gvcf=$2
  local sample=$3
  if [[ "$DELETE_BAM_AFTER_GVCF" != "1" ]]; then
    return 0
  fi
  if ! gvcf_ready "$gvcf" "$sample"; then
    echo "WARNING: not deleting $input_bam because GVCF is incomplete" >&2
    return 1
  fi
  rm -f "$input_bam" "${input_bam}.bai" "${input_bam}.crai" \
    "${input_bam%.*}.bai" "${input_bam}.aria2"
  echo "Removed alignment to free disk: $input_bam"
}

strip_alignment_ext() {
  local name=$1
  name=${name%.bam}
  name=${name%.cram}
  name=${name%.BAM}
  name=${name%.CRAM}
  printf '%s' "$name"
}

process_bam() {
  local input_bam=$1
  local sample_id=${2:-}
  local output_gvcf=${3:-}

  if [[ ! -f "$input_bam" ]]; then
    echo "ERROR: alignment not found: $input_bam" >&2
    return 1
  fi
  if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference FASTA not found: $REF" >&2
    echo "Set REF to human_g1k_v37.fasta (see repository README Drive folder)." >&2
    return 1
  fi
  if [[ ! -f "$KNOWN_SITES" ]]; then
    echo "ERROR: known-sites VCF not found: $KNOWN_SITES" >&2
    return 1
  fi

  local base_name
  base_name=$(strip_alignment_ext "$(basename "$input_bam")")
  sample_id=${sample_id:-$base_name}
  mkdir -p "$RESULTS_DIR"
  output_gvcf=${output_gvcf:-"$RESULTS_DIR/${sample_id}.g.vcf.gz"}
  mkdir -p "$(dirname "$output_gvcf")"

  if [[ "$SKIP_EXISTING" == "1" && -s "$output_gvcf" && ( -f "${output_gvcf}.tbi" || -f "${output_gvcf}.idx" ) ]]; then
    echo "Skipping $sample_id (GVCF already exists): $output_gvcf"
    maybe_delete_bam "$input_bam" "$output_gvcf" "$sample_id" || true
    return 0
  fi

  local interval_args=()
  if [[ -n "$INTERVALS_BED" ]]; then
    if [[ ! -f "$INTERVALS_BED" ]]; then
      echo "ERROR: INTERVALS_BED not found: $INTERVALS_BED" >&2
      echo "Use the union of UC and JPT capture/exome targets, not the UC gene BED." >&2
      return 1
    fi
    interval_args=(-L "$INTERVALS_BED")
  else
    echo "WARNING: INTERVALS_BED unset; calling whole chromosomes. For Wave 2 set the capture-union BED." >&2
    interval_args=(-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y -L MT)
  fi

  ensure_index "$KNOWN_SITES"

  local work_dir="$RESULTS_DIR/work/${sample_id}"
  mkdir -p "$work_dir"
  local log="$RESULTS_DIR/logs/${sample_id}.log"
  mkdir -p "$(dirname "$log")"

  local fixed_bam="$work_dir/read_groups.bam"
  local recal_table="$work_dir/bqsr.table"
  local recalibrated_bam="$work_dir/recalibrated.bam"

  {
    echo "-------------------------------------------------------"
    echo "GVCF calling: $sample_id"
    echo "  input:  $input_bam"
    echo "  output: $output_gvcf"
    echo "  ref:    $REF"
    echo "-------------------------------------------------------"

    if [[ ! -f "${input_bam}.bai" && ! -f "${input_bam}.crai" && ! -f "${input_bam%.*}.bai" ]]; then
      echo "Indexing input alignment (no companion BAI/CRAI found)"
      samtools index "$input_bam"
    fi

    gatk_cmd AddOrReplaceReadGroups \
      -I "$input_bam" \
      -O "$fixed_bam" \
      -R "$REF" \
      -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM "$sample_id" \
      --VALIDATION_STRINGENCY LENIENT
    samtools index "$fixed_bam"

    gatk_cmd BaseRecalibrator \
      -R "$REF" \
      -I "$fixed_bam" \
      --known-sites "$KNOWN_SITES" \
      "${interval_args[@]}" \
      -O "$recal_table" \
      --read-validation-stringency LENIENT

    gatk_cmd ApplyBQSR \
      -R "$REF" \
      -I "$fixed_bam" \
      --bqsr-recal-file "$recal_table" \
      -O "$recalibrated_bam" \
      --read-validation-stringency LENIENT
    samtools index "$recalibrated_bam"

    gatk_cmd HaplotypeCaller \
      -R "$REF" \
      -I "$recalibrated_bam" \
      "${interval_args[@]}" \
      -ERC GVCF \
      --native-pair-hmm-threads "$HC_THREADS" \
      -O "$output_gvcf" \
      --read-validation-stringency LENIENT

    if [[ ! -f "${output_gvcf}.tbi" && ! -f "${output_gvcf}.idx" ]]; then
      ensure_index "$output_gvcf"
    fi

    if [[ "$KEEP_WORK" != "1" ]]; then
      rm -rf "$work_dir"
    fi

    echo "Finished $sample_id -> $output_gvcf"
  } 2>&1 | tee -a "$log"

  maybe_delete_bam "$input_bam" "$output_gvcf" "$sample_id"
}

usage() {
  cat <<EOF
Usage:
  $0 -f path/to/sample.bam [-s SAMPLE_ID] [-o OUT.g.vcf.gz] [-L INTERVALS.bed]
  $0 -d path/to/bam_directory

Wave 2 GVCF caller. Outputs SAMPLE.g.vcf.gz under RESULTS_DIR (default: ./gvcf_results).
Deletes the input BAM only if DELETE_BAM_AFTER_GVCF=1 and the GVCF is indexed.

Environment:
  REF INTERVALS_BED KNOWN_SITES RESULTS_DIR GATK_JAVA_OPTS HC_THREADS SKIP_EXISTING
  DELETE_BAM_AFTER_GVCF=1  # delete each BAM only after its GVCF is indexed
EOF
  exit 1
}

sample_id_arg=""
output_arg=""
mode=""
target=""

while getopts "f:d:s:o:L:h" opt; do
  case $opt in
    f) mode="file"; target=$OPTARG ;;
    d) mode="dir"; target=$OPTARG ;;
    s) sample_id_arg=$OPTARG ;;
    o) output_arg=$OPTARG ;;
    L) INTERVALS_BED=$OPTARG ;;
    *) usage ;;
  esac
done

if [[ -z "$mode" ]]; then
  usage
fi

if [[ "$mode" == "file" ]]; then
  process_bam "$target" "$sample_id_arg" "$output_arg"
else
  if [[ ! -d "$target" ]]; then
    echo "ERROR: directory not found: $target" >&2
    exit 1
  fi
  shopt -s nullglob
  local_files=("$target"/*.bam "$target"/*.cram)
  if [[ ${#local_files[@]} -eq 0 ]]; then
    echo "ERROR: no BAM/CRAM files in $target" >&2
    exit 1
  fi
  for file in "${local_files[@]}"; do
    if [[ "$file" == *_fixed.bam || "$file" == *_recalibrated.bam ]]; then
      continue
    fi
    process_bam "$file"
  done
fi
