# BAM / gVCF checklist (Wave 2)

This build does **not** download alignments by itself. Joint genotyping starts only after every box below is true.

## Blocking inputs

- [ ] EGA account + `pyega3` credentials for `EGAD00001005237` (74 modelled germline WES files in `paper/config/wave2_uc_ega_files.tsv`)
- [ ] 1000 Genomes JPT exome BAM files for the 70 modelled references (`paper/config/wave2_jpt_urls.txt`; do not use the 104-URL `jpt_mapped_only.txt` pool unless `ALL_JPT=1`)
- [ ] GRCh37 reference FASTA + index (repository README Drive folder)
- [ ] dbSNP VCF for BQSR
- [ ] Capture / target BED for the **union** of UC and JPT exome kits, not `config/targets.bed` (UC genes only)
- [ ] Disk and CPU for 144 GVCF calls + GenomicsDB import

## What was wrong in the original bash scripts

`download_and_process.sh` + `variant_calling.sh` + `merge_vcf.sh` downloaded JPT BAMs, called **variant-only** HaplotypeCaller (no `-ERC GVCF`), deleted the BAM, and `bcftools merge`d the VCFs. That is the Wave 1 artifact. Those scripts now emit GVCFs, restrict JPT to the modelled 70, fetch only the 74 modelled UC germline WES files via `pyega3`, refuse the merge, and still delete each BAM **after** the GVCF exists and is indexed (one BAM on disk at a time).

UC files are EGAF accessions, not HTTP URLs. The dataset also contains crypt, tumour, non-UC `CP*_germline`, and RNA-seq BAMs; do not `pyega3 fetch EGAD00001005237`.

## Commands

```bash
# 0. Environment (edit paths first)
source paper/config/wave2.env.example

# 1. UC (pyega3, 74 germline WES) then JPT (HTTPS, 70). One BAM at a time.
#    COHORT=uc or COHORT=jpt to run a single arm. DRY_RUN=1 lists IDs only.
source paper/config/wave2.env.example   # edit REF, INTERVALS_BED, PYEGA3_CREDENTIALS
bash download_and_process.sh

# 2. Sample map + joint genotyping
uc-genepy-ml/scripts/upstream/build_sample_map.sh "$RESULTS_DIR" sample_map.tsv
uc-genepy-ml/scripts/upstream/joint_genotype.sh sample_map.tsv "$INTERVALS_BED" \
  "$REF" data/genomicsdb data/joint_called.vcf.gz
# Use data/joint_called.filtered_nocall.vcf.gz for downstream masks (LowGQ -> ./.)

# 3. Mutual callable mask (keeps all samples; site intersection only)
uc-genepy-ml/scripts/upstream/callable_mask.sh \
  data/joint_called.filtered_nocall.vcf.gz \
  paper/config/wave2_uc_samples.txt paper/config/wave2_jpt_samples.txt \
  8 0.88 data/joint_callable.vcf.gz

# 4. Drop sites whose missingness predicts source
uc-genepy-ml/scripts/upstream/filter_source_missingness.sh data/joint_callable.vcf.gz \
  paper/config/wave2_uc_samples.txt paper/config/wave2_jpt_samples.txt \
  data/joint_callable_balanced.vcf.gz

# 5. Off-panel ancestry PCA (exclude UC target BED)
uc-genepy-ml/scripts/upstream/ancestry_pca.sh data/joint_callable_balanced.vcf.gz \
  uc-genepy-ml/config/targets.bed data/ancestry/offpanel 0.05
```

Then rebuild GenePy from the joint-called VCF (do not merge variant-only VCFs) and run `analysis/scripts/paper/wave2_nested_panels.py`.

## Do not

- Call only the 488-gene UC BED.
- Treat `variant_calling.sh` (old behaviour) or `scripts/upstream/call_variants.sh` as joint genotyping. Variant-only mode now refuses unless `ALLOW_VARIANT_ONLY=1`.
- Run `merge_vcf.sh` / `merge_variant_vcfs.sh` for Wave 2.
- Compute PCA from the 537 UC-panel SNVs.
- Expand the gene list on the current variant-only merge.
- Delete a BAM before its GVCF exists and is indexed. `DELETE_BAM_AFTER_GVCF` defaults to 1 and only runs after that check; set it to 0 if you have space to keep alignments.
- Fetch the whole EGA dataset or crypt/tumour/RNA-seq BAMs. Use the 74 EGAF IDs in `paper/config/wave2_uc_ega_ids.txt`.
