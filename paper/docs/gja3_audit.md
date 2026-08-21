# GJA3 exclusion audit

**Wave 1 finding:** restoring *GJA3* raised locked-RF repetition-level AUROC from 0.765 to 0.787. Dropping it therefore did **not** inflate the published result; if anything it slightly reduced discrimination.

## What the repository says

- `TECHNICAL_EXCLUSIONS = {"GJA3"}` in `create_genepy_matrix.py`
- Public docs: “prespecified technical exclusion”
- Methods workflow SVG (thesis): “GJA3 removed after technical audit”
- Legacy `extract_top_variants.py` listed GJA3 among SHAP-focused genes (`GALC`, `SCN7A`, `IP6K3`, `GJA3`) on an older 62-gene matrix

## What is missing

No dated protocol, coverage metric, or HWE/call-rate threshold that uniquely identifies *GJA3* before modelling. Timing relative to seeing importance is not auditable.

## Wave 1 reporting rule

Treat *GJA3* as a **sensitivity analysis**, not a prespecified filter, until a contemporaneous technical note is located. Report both 215- and 216-gene matrices. Do not describe the gene as a biomarker.
