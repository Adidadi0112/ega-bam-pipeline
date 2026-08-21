# Paper reconstruction (Wave 1)

This directory holds the source-aware reconstruction of EXP-001. It is **not** a conversion of thesis Chapters 1–3.

## Start here

1. [protocol/CLAIMS_AND_PROTOCOL.md](protocol/CLAIMS_AND_PROTOCOL.md) — locked titles, stop rules, Wave 1 settings  
2. [docs/wave1_decision.md](docs/wave1_decision.md) — provisional **Paper A**  
3. [manuscript/wave1_draft.md](manuscript/wave1_draft.md) — draft article  
4. [results/wave1/](results/wave1/) — tables and figures  

## Wave 1 result in one sentence

A QC/source Random Forest (AUROC 0.817) outperformed the 215-gene UC panel (0.765); making missing genotypes missing rather than zero dropped AUROC to 0.562.

## What this build did not do

No BAM download, no GVCF joint calling, no autoimmune/WES scoring (current VCF is UC-targeted). Wave 2 checklist: [docs/bam_gvcf_checklist.md](docs/bam_gvcf_checklist.md).

## Scripts

| Script | Role |
|--------|------|
| `analysis/scripts/paper/wave1_identifiability.py` | QC vs gene encodings vs random subsets |
| `analysis/scripts/paper/wave1_pathways.py` | Pathway-mean features + rank enrichment |
| `analysis/scripts/paper/plot_wave1.py` | Figures |
| `analysis/scripts/paper/wave2_nested_panels.py` | Stafford-style panels after joint GenePy exists |
| `download_and_process.sh` | JPT BAM download + Wave 2 GVCF calling (70 modelled samples) |
| `variant_calling.sh` | Per-sample BQSR + HaplotypeCaller `-ERC GVCF` |
| `uc-genepy-ml/scripts/upstream/call_gvcf.sh` | Same GVCF caller, portable CLI |
| `uc-genepy-ml/scripts/upstream/joint_genotype.sh` | GenomicsDB import + GenotypeGVCFs |
