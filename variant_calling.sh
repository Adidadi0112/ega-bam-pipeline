# download reference files
# wget -P /Users/adamgruda/Projects/uc_genome_analysis/ref https://storage.googleapis.com/gatk-legacy-bundles/b37/human_g1k_v37.fasta.gz
# gunzip /Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.fasta.gz

# index ref - .fai file before running haplotype caller
# samtools faidx /Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.fasta


# ref dict - .dict file before running haplotype caller
# gatk CreateSequenceDictionary R=/Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.fasta O=/Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.dict

# wget -P /Users/adamgruda/Projects/uc_genome_analysis/ref https://storage.googleapis.com/gatk-legacy-bundles/b37/dbsnp_138.b37.vcf.gz
# gunzip /Users/adamgruda/Projects/uc_genome_analysis/ref/dbsnp_138.b37.vcf.gz

ref="/Users/adamgruda/Projects/uc_genome_analysis/ref/human_g1k_v37.fasta"
known_sites="/Users/adamgruda/Projects/uc_genome_analysis/ref/dbsnp_138.b37.vcf"
input_bam="/Users/adamgruda/Projects/uc_genome_analysis/data/EGAF00002549929/UC39-1_i.bam"
fixed_bam="/Users/adamgruda/Projects/uc_genome_analysis/data/EGAF00002549929/UC39-1_fixed.bam"
recal_table="/Users/adamgruda/Projects/uc_genome_analysis/results/recal_data.table"
recalibrated_bam="/Users/adamgruda/Projects/uc_genome_analysis/data/EGAF00002549929/UC39-1_recalibrated.bam"
intervals="-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y -L MT"

if [ ! -f "${known_sites}.idx" ]; then
    gatk IndexFeatureFile -I ${known_sites}
fi

# gatk AddOrReplaceReadGroups \
#    -I ${input_bam} \
#    -O ${fixed_bam} \
#    -RGID 1 \
#    -RGLB lib1 \
#    -RGPL illumina \
#    -RGPU unit1 \
#    -RGSM UC39-1

# samtools index ${fixed_bam}

# gatk BaseRecalibrator \
#    -R ${ref} \
#    -I ${fixed_bam} \
#    --known-sites ${known_sites} \
#    -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 \
#    -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 \
#    -L 21 -L 22 -L X -L Y -L MT \
#    -O ${recal_table}

# gatk ApplyBQSR \
#    -R ${ref} \
#    -I ${fixed_bam} \
#    --bqsr-recal-file ${recal_table} \
#    -O ${recalibrated_bam}


echo "Step 5: Indexing recalibrated BAM..."
samtools index /Users/adamgruda/Projects/uc_genome_analysis/data/EGAF00002549929/UC39-1_recalibrated.bam

echo "Step 6: Running HaplotypeCaller..."
gatk HaplotypeCaller \
   -R ${ref} \
   -I /Users/adamgruda/Projects/uc_genome_analysis/data/EGAF00002549929/UC39-1_recalibrated.bam \
   ${intervals} \
   -O /Users/adamgruda/Projects/uc_genome_analysis/results/UC39-1.vcf