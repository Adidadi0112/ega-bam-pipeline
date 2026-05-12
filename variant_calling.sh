#!/bin/bash
set -exuo pipefail

# --- KONFIGURACJA ŚCIEŻEK ---
REF="/home/adam/projects/ega-bam-pipeline/ref/human_g1k_v37.fasta"
KNOWN_SITES="/home/adam/projects/ega-bam-pipeline/ref/dbsnp_138.b37.vcf"
RESULTS_DIR="/home/adam/projects/ega-bam-pipeline/control_results"

INTERVALS="-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y -L MT"

mkdir -p "$RESULTS_DIR"

# --- FUNKCJA PRZETWARZAJĄCA POJEDYNCZY PLIK ---
process_bam() {
    local input_bam=$1
    local base_name=$(basename "$input_bam" .bam)
    local dir_name=$(dirname "$input_bam")
    
    echo "-------------------------------------------------------"
    echo "Przetwarzanie próbki: $base_name"
    echo "-------------------------------------------------------"

    # Definicja plików tymczasowych
    local fixed_bam="${dir_name}/${base_name}_fixed.bam"
    local recal_table="${RESULTS_DIR}/${base_name}_recal_data.table"
    local recalibrated_bam="${dir_name}/${base_name}_recalibrated.bam"
    local output_vcf="${RESULTS_DIR}/${base_name}.vcf"

    # 1. AddOrReplaceReadGroups
    gatk AddOrReplaceReadGroups \
       -I "$input_bam" -O "$fixed_bam" \
       -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM "$base_name" \
       --VALIDATION_STRINGENCY LENIENT
    
    samtools index "$fixed_bam"

    # 2. BaseRecalibrator
    gatk BaseRecalibrator \
       -R "$REF" -I "$fixed_bam" --known-sites "$KNOWN_SITES" \
       $INTERVALS -O "$recal_table" \
       --read-validation-stringency LENIENT

    # 3. ApplyBQSR
    gatk ApplyBQSR \
       -R "$REF" -I "$fixed_bam" --bqsr-recal-file "$recal_table" \
       -O "$recalibrated_bam" \
       --read-validation-stringency LENIENT
    
    # 4. Indexing recalibrated BAM
    samtools index "$recalibrated_bam"

    # 5. HaplotypeCaller
    gatk HaplotypeCaller \
       -R "$REF" -I "$recalibrated_bam" \
       $INTERVALS -O "$output_vcf" \
       --read-validation-stringency LENIENT
    
    rm "$fixed_bam" "${fixed_bam}.bai" "$recalibrated_bam" "${recalibrated_bam}.bai" "$recal_table"
    
    echo "Zakończono: $base_name. Wynik w $output_vcf"
}

# --- LOGIKA WYBORU TRYBU ---
usage() {
    echo "Użycie:"
    echo "  $0 -f ścieżka/do/pliku.bam    # Tryb pojedynczego pliku"
    echo "  $0 -d ścieżka/do/folderu      # Tryb batch (wszystkie .bam w folderze)"
    exit 1
}

if [ ! -f "${KNOWN_SITES}.idx" ]; then
    gatk IndexFeatureFile -I "$KNOWN_SITES"
fi

echo "Argumenty: $@"

while getopts "f:d:" opt; do
  case $opt in
    f)
      echo "Tryb pojedynczego pliku: $OPTARG"
      process_bam "$OPTARG"
      ;;
    d)
      echo "Tryb katalogu: $OPTARG"
      for file in "$OPTARG"/*.bam; do
          # Pomijamy pliki, które nie istnieją, oraz pliki tymczasowe wygenerowane przez ten skrypt
          if [ -e "$file" ] && [[ ! "$file" == *"_fixed.bam" ]] && [[ ! "$file" == *"_recalibrated.bam" ]]; then
              process_bam "$file"
          fi
      done
      ;;
    *)
      usage
      ;;
  esac
done

if [ $# -eq 0 ]; then usage; fi
