Stage 03: Technical Harmonization
Overview

Ensuring that technical differences between sequencing runs do not bias the biological results.
What we have (Input)

    all_samples_merged_raw.vcf.gz (Stage 01).

    uc_genes.bed (Stage 02).

What we do (Process)

    Extraction: Extract variants only within the targeted BED regions.

    Depth Filtering: Recode genotypes with low coverage (DP < 5) as missing (./.).

What we get (Output)

    uc_targeted_DP5.vcf.gz: A 45 MB filtered subset of the data.

Theoretical Background

UC samples (~34x) and JPT controls (~24x) were sequenced at different depths. Without harmonization, an ML model might "classify" samples based on sequencing quality rather than actual mutations. Setting a minimum DP threshold levels the playing field.