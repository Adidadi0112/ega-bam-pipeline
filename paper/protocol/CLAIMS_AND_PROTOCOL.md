# Frozen claims and Wave 1 protocol

**Status:** locked before inspection of new AUROC values.  
**Date:** 2026-08-20  
**Experiment:** EXP-001 reconstruction, Wave 1  
**Intended use:** gene-level rare-variant burdens are evaluated as a leakage-safe representation for later digital-twin work. This study tests whether that representation carries identifiable ulcerative-colitis (UC)–related signal after source and callability controls. It is not a diagnostic model.

Claims below are mutually exclusive products. Wave 1 may only support a *provisional* lean toward Paper A. Papers B and C remain forbidden until later waves.

---

## Paper products (locked titles)

### Paper A — methods / cautionary

**Title:** Source and callability artifacts can fully separate two Japanese exome cohorts labelled as ulcerative colitis versus population reference

**Primary claim:** After common downstream processing, discrimination between 74 Kyoto UC exomes and 70 1000 Genomes JPT references is explained by source, callability, or other technical features.

**Venue if selected:** *BMC Research Notes* or a methods note.

### Paper B — source-aware residual signal

**Title:** Gene-level rare-variant burdens distinguish Japanese UC and 1000 Genomes JPT exomes after joint genotyping: a source-aware methodological study

**Primary claim:** After joint genotyping, a mutual callable mask, and technical baselines, some UC-panel signal remains on this 74 vs 70 split. The residual is not identified as disease biology.

**Venue if selected:** *BMC Medical Genomics*.

**Not licensed by Wave 1.** Requires Wave 2 gVCF joint calling.

### Paper C — cautious UC biology

**Title:** Gene-level exome burdens retain UC-associated signal after source-matched blood controls and callability masking

**Primary claim:** Discrimination persists after technical controls *and* after labels are no longer aliased to sequencing source.

**Venue if selected:** *BMC Medical Genomics* or an IBD/genetics journal if replication is solid.

**Not licensed by Wave 1.** Requires Wave 3 same-study blood controls and/or an independent internally matched cohort.

---

## Claims never allowed from EXP-001 alone

- Clinical diagnosis of UC.
- Incremental value of WES beyond EHR (Q04).
- Prediction of current IBD activity (H01).
- Validated biomarkers (including *NDUFAF2*, *S1PR5*, *ADCY3*).
- Transportability to another population.
- “Healthy controls” for JPT samples.

---

## Wave 1 stop / go rules (pre-specified)

Let AUROC denote Random Forest repetition-level mean AUROC on identical outer partitions (RepeatedStratifiedKFold, 5 splits × 10 repeats, `random_state=42`).

1. If QC-only AUROC ≥ UC-panel missing-as-zero AUROC, **provisional Paper A**.
2. If callability-aware AUROC falls to chance (≤ 0.55) while missing-as-zero remains high, **provisional Paper A** (missing-as-zero encoding is the driver).
3. If UC-panel AUROC exceeds both QC-only and the mean of five random 30-gene subsets by ≥ 0.05, record a *residual-on-this-callset* observation. Do **not** upgrade to Paper B. Joint genotyping is still required.
4. If pathway-level AUROC exceeds gene-level AUROC *and* exceeds size-matched random gene-set features, record a representation finding. It remains exploratory and source-confounded.
5. Rank-based enrichment is exploratory. FDR-significant terms are not biomarkers. Genome-wide background on this Open Targets panel is not discovery.

Final A/B/C selection occurs only after Wave 2 (and Wave 3 for C).

---

## Wave 1 analysis protocol

### Data

- Matrix A: `analysis/data/genepy_expanded/genepy_original_missing_as_zero_without_GJA3.csv` (144 × 215). Primary continuity matrix.
- Matrix B: same scores with *GJA3* restored (`genepy_original_missing_as_zero.csv`).
- Matrix C: `genepy_original_callability_aware.csv`, *GJA3* dropped to match Matrix A.
- Callability: `genepy_callability.csv`.
- Labels: 1 = UC (EGA `EGAD00001005237` germline), 0 = JPT identifier containing `.mapped.ILLUMINA.bwa.JPT.exome`.

### Classifier (representation comparison, not algorithm ranking)

- Random Forest only (representation comparison, not algorithm ranking).
- Outer partitions identical to the thesis: `RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)`.
- Locked classifier hyperparameters, chosen before any Wave 1 AUROC was inspected: `n_estimators=150`, `max_depth=5`, `min_samples_leaf=2`, `min_samples_split=4`, `max_features="sqrt"`, `class_weight="balanced"`, `random_state` = outer-split seed. TPE is retained as the thesis algorithm-comparison protocol and as the Wave 2 default; it is not re-run for every Wave 1 matrix.
- Fold-contained median imputation and standardisation. Wave 1 uses all columns of each representation (`k = all`) so that QC, pathway, and gene-panel matrices are compared as whole inputs. Mutual-information k-search remains the thesis algorithm-comparison protocol.
- Primary metric: AUROC. Secondary: balanced accuracy, sensitivity, specificity, precision, F1, average precision from native `predict()`.

### Experiments

1. Matrix A vs Matrix C (missing-as-zero vs callability-aware).
2. Matrix A vs Matrix B (*GJA3* exclusion audit).
3. QC/source baseline: per-sample mean callability, fully-callable gene fraction, total GenePy burden, non-zero gene count, GenePy standard deviation, and VCF-derived missing genotype rate / mean depth / Ti/Tv when the targeted VCF is readable.
4. Five random 30-gene subsets of the 215 genes (seeds 0–4), not off-panel genes.
5. Pathway-level mean GenePy for Reactome and GO Biological Process sets with ≥ 3 modelled genes; comparison with size-matched random gene-set features.
6. Rank-based enrichment: genes ranked by univariate Mann–Whitney U (UC vs JPT), by existing held-out permutation importance and SHAP where available, and by callability-residualized univariate ranks. Tests: g:Profiler ordered query (custom 215-gene background) plus a permutation prerank statistic. Exploratory.

### GJA3

Report with and without. If documentation of a pre-specified technical exclusion cannot be produced, treat the drop as a sensitivity analysis, not a prespecified filter.

---

## Wave 2 / 3 (not executed in this build)

- Joint-call GVCFs over the capture/exome union; intersection callable mask; ancestry PCA on off-panel SNPs.
- Nested UC vs autoimmune vs all-callable panels with size/callability-matched random panels.
- NBDC `hum0201` / `JGAS000199` peripheral-blood controls; All of Us / dbGaP access in parallel.
