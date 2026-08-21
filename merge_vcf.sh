#!/usr/bin/env bash
# Wave 1 leftover: bcftools-merge of variant-only VCFs.
# Wave 2 must jointly genotype GVCFs instead.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOINT="$SCRIPT_DIR/uc-genepy-ml/scripts/upstream/joint_genotype.sh"

cat >&2 <<EOF
ERROR: merge_vcf.sh is the Wave 1 variant-only merge and must not be used for Wave 2.

That path encodes absent records as missing genotypes, which GenePy then scores
as zero. Joint genotyping from GVCFs is required.

Use:
  $JOINT SAMPLE_MAP.tsv INTERVALS.bed REFERENCE.fasta WORKSPACE_DIR OUTPUT.vcf.gz

Build SAMPLE_MAP.tsv with:
  $SCRIPT_DIR/uc-genepy-ml/scripts/upstream/build_sample_map.sh GVCF_DIR sample_map.tsv

Checklist: paper/docs/bam_gvcf_checklist.md
EOF
exit 1
