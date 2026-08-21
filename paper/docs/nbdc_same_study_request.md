# NBDC same-study blood controls (Wave 3 request)

**Status:** draft only. Not part of the Wave 1 local build.  
**Purpose:** repair the EGA-versus-1000 Genomes source split. This is **not** independent external replication.

## Dataset

- NBDC `hum0201` / `JGAS000199`
- Same Kyoto ulcerative-colitis WES study family as EGA `EGAD00001005237`
- Public description: peripheral-blood exomes for a small number of non-neoplastic UC, UC-with-neoplasia, and colonoscopy controls; other samples are organoids or tissue

## Request

Ask the existing EGA/Japanese data contact for:

1. Participant/sample manifest linking `EGAD00001005237` to `JGAS000199`
2. DNA source (peripheral blood vs organoid vs diseased tissue)
3. Permission to use the **peripheral-blood colonoscopy controls** and blood-derived UC samples only

## Analysis rules if access is granted

- Germline classifier: blood vs blood only. Do not mix organoids or lesional tissue.
- Check sample-ID overlap with the 74 modelled UC cases.
- Apply the Wave 2 joint-called, callable-mask pipeline.
- Treat the result as a source-control sensitivity analysis. Sample size is expected to be too small for a definitive performance claim.
- Do not label this as external validation.

## Contact log

- Date requested:
- Contact:
- Outcome:
