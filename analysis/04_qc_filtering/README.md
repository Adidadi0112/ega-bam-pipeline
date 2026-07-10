Stage 04: Quality Control & Sample Cleaning
Overview

Removing technical artifacts and unreliable samples to ensure data integrity.
What we have (Input)

    uc_targeted_DP5.vcf.gz (Stage 03).

    samples_comparison.csv: Statistical report on heterozygosity and coverage.

What we do (Process)

    Outlier Removal: Exclude 36 FFPE samples (chemical artifacts) and failed samples (NA19004).  

    Hard Filtering: Apply thresholds for Variant Quality (QUAL > 30), Genotype Quality (GQ > 20), and Missingness (F_MISSING < 0.2).  

    Validation: Calculate the Ti/Tv ratio to ensure biological realism.

What we get (Output)

    uc_final_clean.vcf.gz: 864 high-quality SNPs across 107 samples.

    samples_to_keep.txt: Final list of validated samples.

Theoretical Background

FFPE (Formalin-Fixed Paraffin-Embedded) samples often introduce artificial C>T transitions. This stage acts as a "biological filter," ensuring that our Ti/Tv ratio (~2.62) matches expected human exome values, confirming that our remaining 864 variants are likely real.