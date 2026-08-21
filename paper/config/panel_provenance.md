# Autoimmune / IBD panel provenance (frozen for Wave 2)

**Not scored in Wave 1.** The current VCF is restricted to the UC Open Targets BED.

## Primary (thesis continuity)

- Source: Open Targets `EFO_0000729` ulcerative colitis, export 30 April 2026, overall score > 0.15
- File: `uc-genepy-ml/config/genes_uc.txt`
- Count: 488 symbols (460 with GRCh37 intervals)

## Expansion (Stafford-style nested comparison)

Stafford et al. 2023 (*J Crohns Colitis* 17:1672–1680) found an HTG autoimmune panel (~1,540 genes) beat an IBD panel (489) and all genes (~15,000) for **CD vs UC in one sequenced cohort**. That result is not copied as a primary analysis here. It motivates a **nested comparison after joint calling**.

Because the HTG EdgeSeq gene list is commercial and may not be redistributable, this freeze uses:

1. The 488 UC Open Targets genes.
2. A documented supplement of innate/adaptive immunity, monogenic IBD, interferon, TNF, autophagy, barrier, and HLA class II genes that are not already in the UC list (`paper/config/autoimmune_supplement_only.txt`).

Combined unique symbols are in `paper/config/autoimmune_ibd_panel.txt`.

**Wave 2 evaluation rule:** the autoimmune panel must beat both the UC panel and a size/callability-matched random panel on the joint-called matrix. Whole-exome (all callable genes) is a control, not the primary representation.

If HTG or another licensed autoimmune list becomes available, replace the supplement and re-freeze hashes before scoring.
