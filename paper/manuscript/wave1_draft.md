# Source and callability artifacts can fully separate two Japanese exome cohorts labelled as ulcerative colitis versus population reference

**Working product:** Paper A (methods / cautionary). Wave 1 results only.  
**Not a diagnostic or biomarker paper.**  
**Intended use:** gene-level rare-variant burdens are evaluated as a leakage-safe representation for later digital-twin work. This draft tests whether that representation carries identifiable UC-related signal after source and callability controls.

## Abstract

**Background.** Gene-level aggregation of whole-exome variants (GenePy) is a candidate endotype-oriented input for later multimodal models. A master’s thesis reported that a Random Forest discriminated 74 Japanese ulcerative colitis (UC) exomes from 70 1000 Genomes JPT population-reference exomes (mean AUROC 0.765). Disease status is perfectly aligned with sequencing source, and genotypes were merged from variant-only VCFs.

**Methods.** We froze three mutually exclusive paper claims before inspecting new AUROCs. On the existing 144 × 215 GenePy matrix we compared, under identical repeated nested partitions (5 × 10, seed 42), (i) missing-as-zero gene burdens, (ii) a callability-aware encoding in which any missing locus sets the gene to missing, (iii) QC/source features only, (iv) random 30-gene subsets of the UC panel, and (v) pathway-mean GenePy features. Rank-based enrichment used a custom 215-gene background.

**Results.** A QC/source Random Forest outperformed the UC-panel model (AUROC 0.817 ± 0.011 vs 0.765 ± 0.030). Callability-aware encoding collapsed discrimination to 0.562 ± 0.032. Restoring *GJA3* slightly increased AUROC (0.787). Pathway-mean features (AUROC 0.509) did not beat size-matched random gene-set features (0.624) or gene-level scores. No GO Biological Process, Reactome, or KEGG term was significant in an ordered g:Profiler query after FDR correction.

**Conclusions.** On this variant-only merge, technical features and missing-as-zero encoding are sufficient to match or beat gene-panel discrimination. The result is internal cohort separation, not evidence of UC biology. Joint genotyping and source-matched controls remain required before any residual-signal or biological claim.

**Keywords:** ulcerative colitis; whole-exome sequencing; GenePy; callability; source confounding; nested cross-validation

## Introduction

Rare-variant whole-exome data are high-dimensional relative to typical clinical sample sizes. GenePy collapses per-gene zygosity, population frequency, and predicted deleteriousness into a sample-by-gene matrix (Mossotto et al. 2019). Stafford et al. (2023) used that representation to classify Crohn’s disease versus UC *within one sequenced IBD cohort*, and found that an autoimmune gene panel outperformed both an IBD panel and all exome genes.

The present cohort is different. All 74 UC samples come from EGA `EGAD00001005237` (Kyoto germline UC) and all 70 comparators from 1000 Genomes JPT exomes. Labels equal source. Independently called variant-only VCFs were merged, so an absent record, a no-call, and a homozygous-reference genotype can share a gene-level value of zero. Depth already differed by source (UC ~34×, JPT ~24×). Under those conditions, a gene-panel AUROC cannot be read as disease biology.

This Wave 1 study therefore asks a validity question: **how much of the published discrimination is recoverable from callability, missingness, and other QC features, and does gene-set aggregation add identifiable structure inside the UC panel?** Incremental value of WES over EHR, and prediction of IBD activity, are out of scope.

## Methods

### Claims freeze

Three paper products and stop/go rules were locked in `paper/protocol/CLAIMS_AND_PROTOCOL.md` before new AUROCs were computed. JPT samples are described as population-reference, not healthy controls.

### Cohort and GenePy matrix

The analysis uses the thesis modelling matrix: 74 UC and 70 JPT individuals; Open Targets UC panel (`EFO_0000729`, score > 0.15, 488 symbols); 537 biallelic SNVs; 215 gene-level predictors after dropping *GJA3*. GenePy used CADD RawScore scaled to [0, 1] and East Asian then global allele frequencies (`gnomADe_EAS_AF` → `gnomADe_AF` → `1000G_EAS_AF` → `1000G_AF` → 10⁻⁵ floor). Missing genotypes contributed zero burden in the primary matrix. A callability-aware matrix set a gene to missing if any mapped locus had a missing or unsupported genotype.

### Representations compared

All used Random Forest (`n_estimators=150`, `max_depth=5`, `min_samples_leaf=2`, `class_weight=balanced`) on `RepeatedStratifiedKFold` partitions identical to the thesis (5 splits, 10 repeats, `random_state=42`). Wave 1 used all columns of each representation (`k = all`) so that QC and gene panels were compared as whole inputs. Thesis TPE algorithm ranking was not repeated.

1. Missing-as-zero without *GJA3* (primary continuity matrix).
2. Missing-as-zero with *GJA3*.
3. Callability-aware matrix without *GJA3*.
4. QC/source features: mean callability, fully callable gene fraction, total GenePy burden, non-zero gene count, GenePy standard deviation, and targeted-VCF missing rate, mean depth, Ti/Tv, and called-site count.
5. Three random 30-gene subsets of the 215 genes (not off-panel genes; the current VCF is UC-targeted).
6. Pathway-mean GenePy for Enrichr Reactome 2022, GO Biological Process 2023, and KEGG 2021 sets with ≥3 modelled genes (712 features), versus 24 size-matched random gene-set means.

### Enrichment

Genes were ranked by Mann–Whitney *p* (UC vs JPT), by existing held-out SHAP and permutation importance, and by a callability-residualized mean difference. Tests: permutation prerank enrichment and g:Profiler ordered query with the 215 modelled genes as custom background. FDR < 0.05 and ≥2 overlapping genes were required for a reported hit.

## Results

### QC and source audit

Univariate AUROC for individual QC features ranged from 0.51 (total GenePy burden) to 0.70 (Ti/Tv). UC samples had higher mean depth (81 vs 70), higher mean gene callability (0.972 vs 0.944), and lower missing genotype rate (0.027 vs 0.059). Total gene burden did not separate classes (AUROC 0.51). A Random Forest using only these QC features achieved repetition-level AUROC **0.817 ± 0.011**, higher than the UC-panel gene model.

### Gene-panel encodings

| Representation | AUROC mean ± SD |
|----------------|-----------------|
| QC / source features | 0.817 ± 0.011 |
| UC panel + *GJA3* (missing-as-zero) | 0.787 ± 0.013 |
| UC panel without *GJA3* (thesis matrix) | 0.765 ± 0.030 |
| Random 30-gene subset (seed 0) | 0.616 ± 0.030 |
| Callability-aware without *GJA3* | 0.562 ± 0.032 |
| Random 30-gene subsets (seeds 1–2) | 0.560, 0.545 |

The locked-RF AUROC on the thesis matrix (0.765) matched the published TPE Random Forest (0.765 ± 0.027) to two decimals, indicating that algorithm ranking is not required to recover the discrimination. Treating missing loci as missing rather than zero removed most of that discrimination. The 215-gene panel still beat random 30-gene subsets drawn from the same panel; that comparison cannot separate disease genes from “more opportunities to encode source-specific missingness.”

*GJA3* was documented as a technical exclusion without an auditable pre-modelling rationale. Restoring it increased AUROC slightly. It is reported as a sensitivity analysis (`paper/docs/gja3_audit.md`).

### Pathway representation and enrichment

Mean GenePy across 712 Reactome/GO/KEGG sets produced AUROC 0.509 ± 0.030, below 24 size-matched random gene-set features (0.624 ± 0.028) and below gene-level scores. Ordered g:Profiler returned **no** FDR-significant term. Prerank permutation tests yielded minimum BH FDR ≈ 0.30. The thesis ORA null is therefore not rescued by rank-based methods on this panel.

## Discussion

Wave 1 supports Paper A. The principal uncertainty in the thesis was never ordinary ML overfitting: preprocessing was fold-contained. It was **non-identifiability of disease and source**. QC features alone outclassify the gene panel, and the gene panel’s advantage over missing-aware encoding is exactly the encoding that equates “not called” with “no rare burden.”

Pathway aggregation did not create an endotype-oriented advantage on this callset. That is consistent with aggregating a source-confounded gene matrix: if the predictive structure is missingness, averaging genes into Reactome sets can dilute rather than concentrate it. Stafford’s autoimmune-panel result was obtained in a same-study CD vs UC design with joint genotyping and an order-of-magnitude more samples; it is not a licence to treat whole-exome GenePy as the primary analysis here.

Limitations of Wave 1 itself: no joint genotyping; no off-panel ancestry PCA; random “negative panels” are subsets of UC genes, not matched non-UC genes; locked RF hyperparameters differ from thesis TPE. None of these limitations reverse the QC vs gene-panel ordering.

## Wave 2 / 3 (not done)

Scripts and checklists are in `uc-genepy-ml/scripts/upstream/` and `paper/docs/`. Joint-called GVCFs, a mutual callable mask, nested UC vs autoimmune vs all-callable panels, NBDC same-study blood controls, and All of Us/dbGaP access are required before Paper B or C language.

## Data and code

Controlled EGA and 1000 Genomes alignments are not redistributed. Aggregate Wave 1 tables and figures: `paper/results/wave1/`. Protocol: `paper/protocol/CLAIMS_AND_PROTOCOL.md`. Public modelling package: https://github.com/Adidadi0112/uc-genepy-ml

## TRIPOD+AI (Wave 1 self-audit, abbreviated)

| Item | Status |
|------|--------|
| Title identifies the study as methodological cohort discrimination | Yes |
| Abstract states source confounding and intended use | Yes |
| Intended use: representation validity, not diagnosis | Yes |
| Participants: 74 UC EGA + 70 JPT; labels = source | Yes |
| Predictors: GenePy and QC features; no EHR | Yes |
| Outcome: research-cohort membership, not clinical UC status of JPT | Yes |
| Sample size: n=144, reused across 10 CV repeats | Yes |
| Missing data: missing-as-zero vs callability-aware compared | Yes |
| Analytical methods: repeated nested partitions, locked RF | Yes |
| Overfitting: fold-contained processing; residual risk is source aliasing | Yes |
| Fairness / relatedness: not assessed (no kinship); wording avoids “unrelated” | Yes |
| Validation: internal only; external required for transport | Yes |
| Model output: AUROC with repetition SD; native threshold metrics secondary | Yes |
| Interpretability: QC univariate AUROC; gene ranks exploratory | Yes |
| Availability: code and aggregate outputs; individual genotypes controlled | Yes |

PROBAST+AI: high risk of bias for any *clinical prediction* use, because predictors and outcome are aliased to data source. That risk is the finding, not a defect to be hidden.
