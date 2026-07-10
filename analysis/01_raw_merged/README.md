Stage 01: Raw VCF Merging
Overview

The starting point of the project involves consolidating individual genomic files into a single dataset for joint analysis.
What we have (Input)

    Individual VCF files from two distinct cohorts: UC patients and JPT (Japanese) controls.

What we do (Process)

    Merge: Use bcftools merge to combine samples into a single multi-sample VCF.

    Indexing: Generate .tbi files for fast random access.

What we get (Output)

    all_samples_merged_raw.vcf.gz: A unified 1.9 GB file containing all 144 initial samples.

Theoretical Background

In population genetics, joint calling or merging is necessary to identify which variants are shared and which are unique to specific groups. By merging first, we ensure that every genomic position is accounted for across all individuals, preparing the data for frequency comparisons.