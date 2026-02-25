#!/bin/bash

# --- KONFIGURACJA ŚCIEŻEK ---
REF="/Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.fasta"
KNOWN_SITES="/Users/adamgruda/Projects/uc_genome_analysis/ref/dbsnp_138.b37.vcf"
RESULTS_DIR="/Users/adamgruda/Projects/uc_genome_analysis/results"
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
       -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM "$base_name"
    
    # Usuwamy oryginalny BAM po utworzeniu 'fixed' (zgodnie z prośbą)
    rm "$input_bam"
    samtools index "$fixed_bam"

    # 2. BaseRecalibrator
    gatk BaseRecalibrator \
       -R "$REF" -I "$fixed_bam" --known-sites "$KNOWN_SITES" \
       $INTERVALS -O "$recal_table"

    # 3. ApplyBQSR
    gatk ApplyBQSR \
       -R "$REF" -I "$fixed_bam" --bqsr-recal-file "$recal_table" \
       -O "$recalibrated_bam"
    
    # Usuwamy fixed_bam (i jego indeks) po BQSR
    rm "$fixed_bam" "${fixed_bam}.bai"

    # 4. Indexing recalibrated BAM
    samtools index "$recalibrated_bam"

    # 5. HaplotypeCaller
    gatk HaplotypeCaller \
       -R "$REF" -I "$recalibrated_bam" \
       $INTERVALS -O "$output_vcf"

    # 6. Czyszczenie końcowe próbki
    rm "$recalibrated_bam" "${recalibrated_bam}.bai" "$recal_table"
    
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

while getopts "f:d:" opt; do
  case $opt in
    f)
      process_bam "$OPTARG"
      ;;
    d)
      for file in "$OPTARG"/*.bam; do
          [ -e "$file" ] || continue
          process_bam "$file"
      done
      ;;
    *)
      usage
      ;;
  esac
done

if [ $# -eq 0 ]; then usage; fi