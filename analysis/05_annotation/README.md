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

## Exact final-VCF CADD input

The original `cadd_scores.tsv` covers only 61 of the 537 biallelic SNVs in
`uc_vqsr_rare_included.vcf.gz`. Generate a sample-free VCF containing the exact
final alleles with:

```bash
python analysis/scripts/prep/prepare_final_vcf_cadd.py
```

Score `final_biallelic_snvs_for_cadd.vcf` with CADD GRCh37-v1.6 and save the
tabular result as `cadd_scores_final_vcf.gz`. The CADD website may give the
plain tabular file a `.gz` suffix; the pipeline accepts that output directly. Running the preparation script
again validates that all 537 alleles have a CADD RawScore and PHRED value.

The expanded GenePy builder uses this exact CADD result and the frequency
priority `gnomADe_EAS_AF`, `gnomADe_AF`, `EAS_AF`, `AF`, then `1e-5` floor.
