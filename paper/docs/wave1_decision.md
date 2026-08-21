# Wave 1 decision (provisional)

**Date:** 2026-08-20  
**Rule set:** `paper/protocol/CLAIMS_AND_PROTOCOL.md` (locked before these AUROCs were inspected)

## Pre-specified rules applied

| Rule | Observation | Decision |
|------|-------------|----------|
| 1. QC-only AUROC ≥ UC-panel missing-as-zero | QC 0.817 ± 0.011 vs UC panel 0.765 ± 0.030 | **Provisional Paper A** |
| 2. Callability-aware AUROC ≤ 0.55 while missing-as-zero stays high | Callability-aware 0.562 ± 0.032 vs missing-as-zero 0.765 | **Provisional Paper A** (missing-as-zero encoding is a major driver) |
| 3. UC panel beats random 30-gene subsets by ≥ 0.05 | Yes (0.765 vs mean ~0.574) | Residual *on this callset* only. **Not Paper B.** Joint genotyping still required. |
| 4. Pathway AUROC exceeds gene-level and random gene-set features | Pathway 0.509; random gene-set features 0.624; gene-level 0.765 | No endotype-representation gain on this callset |
| 5. Rank-based enrichment | g:Profiler ordered: 0 FDR-significant terms. Prerank BH FDR ≈ 0.30 at the top | Exploratory null, as expected on a UC-preselected panel |

## Product selected for drafting

**Paper A** (methods / cautionary), Wave 1 evidence only.

Paper B remains closed until GVCF joint calling. Paper C remains closed until source-matched blood controls or an independent internally matched cohort.

## What Wave 1 does *not* establish

- That there is no UC biology in these exomes. It establishes that **on the current variant-only merge**, technical features and missing-as-zero encoding are sufficient to match or beat the gene-panel model.
- Incremental WES value over EHR (Q04) or IBD activity prediction (H01).

## Next required experiment

Wave 2: BAMs → GVCFs → joint genotype → callable mask → repeat QC vs gene-panel comparison. If QC still wins, Paper A is the submission product.
