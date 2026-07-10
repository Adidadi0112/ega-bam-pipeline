Stage 05: Functional Annotation (VEP/CADD)
Overview

The transition from "genomic coordinates" to "biological impact scores."
What we have (Input)

    uc_for_vep.vcf: Unique variant positions.

    Ensembl VEP Online (GRCh37).

What we do (Process)

    Pathogenicity Scoring: Annotate variants using CADD (v1.6) and dbNSFP (SIFT, PolyPhen, REVEL).

    Population Frequency: Add gnomAD and ALFA East Asian allele frequencies.

What we get (Output)

    vep_results.txt: A comprehensive table containing functional weights for every variant.

Theoretical Background

To calculate the GenePy (Gene Pathogenicity) score, we need a weight (W) for each variant. CADD (Combined Annotation Dependent Depletion) provides a PHRED-scaled score that integrates multiple layers of information (conservation, epigenetics, protein structure) into a single metric of deleteriousness.
GenePy_Score=∑(CADD_PHRED×Genotype)