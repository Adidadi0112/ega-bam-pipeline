"""Create deterministic GenePy matrices using the original score definition."""

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import gzip
import json
import math

import pandas as pd


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
VEP_FILE = Path("analysis/05_annotation/vep_results_vqsr.txt")
CADD_FILE = Path("analysis/05_annotation/cadd_scores_final_vcf.gz")
VCF_FILE = Path("analysis/02_targeting/uc_vqsr_rare_included.vcf.gz")
PANEL_FILE = Path("analysis/02_targeting/genes_uc.txt")
OUTPUT_DIR = Path("analysis/data/genepy_expanded")

MAF_FLOOR = 1e-5
CADD_RAW_MIN = -7.53
CADD_RAW_MAX = 35.79

PTV_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
}
TECHNICAL_EXCLUSIONS = {"GJA3"}


def file_sha256(path):
    digest = sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text_auto(path):
    """Open text based on file contents rather than a possibly misleading suffix."""
    with open(path, "rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else open(path, "rt")


def find_header(path, prefix):
    with open_text_auto(path) as handle:
        for line_number, line in enumerate(handle):
            if line.startswith(prefix):
                return line_number, line.rstrip().lstrip("#").split("\t")
    raise ValueError(f"Could not find header {prefix!r} in {path}")


def variant_id(chromosome, position, reference, alternative):
    chromosome = str(chromosome).removeprefix("chr")
    return f"{chromosome}_{int(position)}_{reference}/{alternative}"


def load_cadd(path):
    header_line, _ = find_header(path, "#Chrom")
    with open_text_auto(path) as handle:
        cadd = pd.read_csv(handle, sep="\t", skiprows=header_line)
    cadd.columns = [column.lstrip("#").strip() for column in cadd.columns]

    required = {"Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"}
    missing = required.difference(cadd.columns)
    if missing:
        raise ValueError(f"Missing CADD columns: {sorted(missing)}")

    cadd["Variant"] = [
        variant_id(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in zip(
            cadd["Chrom"], cadd["Pos"], cadd["Ref"], cadd["Alt"]
        )
    ]

    # Extended CADD output contains several transcript annotations per variant,
    # but RawScore and PHRED are constant for a given allele in this dataset.
    score_variation = cadd.groupby("Variant")[["RawScore", "PHRED"]].nunique()
    if (score_variation > 1).any().any():
        raise ValueError("CADD scores are not constant across variant annotations")

    return cadd[["Variant", "RawScore", "PHRED"]].drop_duplicates("Variant")


def extract_extra_number(extra, field):
    values = extra.astype("string").str.extract(fr"(?:^|;){field}=([^;]+)")[0]
    return pd.to_numeric(values, errors="coerce")


def load_gene_panel(path):
    symbols = {
        symbol.strip()
        for symbol in path.read_text().replace("\n", ",").split(",")
        if symbol.strip()
    }
    if not symbols:
        raise ValueError(f"No gene symbols found in {path}")
    return symbols


def load_vep(path, wanted_variants=None, allowed_genes=None):
    header_line, header = find_header(path, "#Uploaded_variation")
    chunks = pd.read_csv(
        path,
        sep="\t",
        skiprows=header_line + 1,
        header=None,
        names=header,
        low_memory=False,
        chunksize=100_000,
    )

    retained = []
    for chunk in chunks:
        if wanted_variants is not None:
            normalized = (
                chunk["Uploaded_variation"]
                .astype(str)
                .str.replace(":", "_", regex=False)
                .str.replace("chr", "", regex=False)
            )
            chunk = chunk[normalized.isin(wanted_variants)]
        if not chunk.empty:
            retained.append(chunk)
    if not retained:
        raise ValueError("No VEP annotations matched the final VCF variants")
    vep = pd.concat(retained, ignore_index=True)

    if "Extra" not in vep.columns:
        raise ValueError("VEP output does not contain the Extra column")

    symbol = vep["Extra"].astype("string").str.extract(
        r"(?:^|;)SYMBOL=([^;]+)"
    )[0]
    vep["Gene_Name"] = symbol.fillna(vep["Gene"])
    if allowed_genes is not None:
        vep = vep[vep["Gene_Name"].isin(allowed_genes)].copy()
    frequency_columns = [
        ("gnomADe_EAS_AF", "gnomADe_EAS_AF"),
        ("gnomADe_AF", "gnomADe_AF"),
        ("EAS_AF", "1000G_EAS_AF"),
        ("AF", "1000G_AF"),
    ]
    selected_frequency = pd.Series(pd.NA, index=vep.index, dtype="Float64")
    frequency_source = pd.Series("floor", index=vep.index, dtype="string")
    for vep_field, source_name in frequency_columns:
        candidate = extract_extra_number(vep["Extra"], vep_field)
        use_candidate = selected_frequency.isna() & candidate.notna()
        selected_frequency.loc[use_candidate] = candidate.loc[use_candidate]
        frequency_source.loc[use_candidate] = source_name

    use_floor = selected_frequency.isna() | selected_frequency.le(0)
    selected_frequency.loc[use_floor] = MAF_FLOOR
    frequency_source.loc[use_floor] = "floor"
    vep["MAF"] = selected_frequency.astype(float)
    vep["MAF"] = vep["MAF"].clip(lower=MAF_FLOOR, upper=1.0 - MAF_FLOOR)
    vep["Frequency_Source"] = frequency_source
    vep["Variant"] = (
        vep["Uploaded_variation"]
        .astype(str)
        .str.replace(":", "_", regex=False)
        .str.replace("chr", "", regex=False)
    )
    vep["Is_PTV"] = vep["Consequence"].fillna("").apply(
        lambda value: bool(PTV_TERMS.intersection(str(value).split("&")))
    )

    # A variant may have several transcript consequences for the same gene.
    # Retain one frequency and mark the variant as truncating if any transcript
    # consequence for that gene is protein truncating.
    return (
        vep[["Variant", "Gene_Name", "MAF", "Frequency_Source", "Is_PTV"]]
        .dropna(subset=["Gene_Name"])
        .groupby(["Variant", "Gene_Name"], as_index=False)
        .agg(
            MAF=("MAF", "first"),
            Frequency_Source=("Frequency_Source", "first"),
            Is_PTV=("Is_PTV", "max"),
        )
    )


def scale_cadd_raw(raw_score):
    scaled = (float(raw_score) - CADD_RAW_MIN) / (CADD_RAW_MAX - CADD_RAW_MIN)
    return min(max(scaled, 0.0), 1.0)


def build_variant_mapping(vep, cadd):
    mapping = vep.merge(cadd, on="Variant", how="inner", validate="many_to_one")
    mapping["Deleteriousness"] = mapping["RawScore"].map(scale_cadd_raw)
    mapping.loc[mapping["Is_PTV"], "Deleteriousness"] = 1.0
    mapping = mapping.sort_values(["Variant", "Gene_Name"]).reset_index(drop=True)
    return mapping


def normalize_gt(gt):
    if gt is None:
        return "missing"
    gt = str(gt).replace("|", "/")
    if gt in {".", "./.", ""}:
        return "missing"
    if gt == "0/0":
        return "reference"
    if gt in {"0/1", "1/0"}:
        return "heterozygous"
    if gt in {"1/1", "1"}:
        return "homozygous_alt"
    return "unsupported"


def genepy_contribution(genotype, maf, deleteriousness):
    if genotype == "heterozygous":
        frequency_product = maf * (1.0 - maf)
    elif genotype == "homozygous_alt":
        frequency_product = maf ** 2
    else:
        return 0.0

    frequency_product = max(frequency_product, MAF_FLOOR ** 2)
    return -deleteriousness * math.log10(frequency_product)


def read_samples(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]
    raise ValueError(f"No #CHROM header found in {path}")


def read_biallelic_variant_ids(path):
    variants = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if "," in columns[4]:
                continue
            variants.add(variant_id(columns[0], columns[1], columns[3], columns[4]))
    return variants


def process_vcf(path, mapping):
    samples = read_samples(path)
    genes = sorted(mapping["Gene_Name"].unique())
    score = pd.DataFrame(0.0, index=samples, columns=genes)
    called = pd.DataFrame(0, index=samples, columns=genes, dtype=int)
    total = pd.Series(0, index=genes, dtype=int)
    missing_by_sample_gene = pd.DataFrame(False, index=samples, columns=genes)

    lookup = defaultdict(list)
    for row in mapping.itertuples(index=False):
        lookup[row.Variant].append(row)

    genotype_counts = Counter()
    audit_counts = Counter()
    matched_variants = set()

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            columns = line.rstrip().split("\t")
            chromosome, position, reference, alternatives = (
                columns[0], columns[1], columns[3], columns[4]
            )

            if "," in alternatives:
                audit_counts["multiallelic_vcf_records_excluded"] += 1
                continue

            current_variant = variant_id(
                chromosome, position, reference, alternatives
            )
            if current_variant not in lookup:
                audit_counts["vcf_records_without_annotation_mapping"] += 1
                continue

            format_fields = columns[8].split(":")
            if "GT" not in format_fields:
                audit_counts["mapped_records_without_GT"] += 1
                continue
            gt_index = format_fields.index("GT")
            matched_variants.add(current_variant)

            genotypes = []
            for sample_data in columns[9:]:
                fields = sample_data.split(":")
                gt = fields[gt_index] if gt_index < len(fields) else None
                genotype_counts[str(gt)] += 1
                genotype = normalize_gt(gt)
                if genotype == "unsupported":
                    audit_counts["unsupported_genotype_calls"] += 1
                genotypes.append(genotype)

            for annotation in lookup[current_variant]:
                gene = annotation.Gene_Name
                total[gene] += 1

                for sample, genotype in zip(samples, genotypes):
                    if genotype == "missing":
                        missing_by_sample_gene.at[sample, gene] = True
                        continue
                    if genotype == "unsupported":
                        missing_by_sample_gene.at[sample, gene] = True
                        continue

                    called.at[sample, gene] += 1
                    score.at[sample, gene] += genepy_contribution(
                        genotype,
                        annotation.MAF,
                        annotation.Deleteriousness,
                    )

    active_genes = total[total > 0].index.tolist()
    score = score[active_genes]
    called = called[active_genes]
    missing_by_sample_gene = missing_by_sample_gene[active_genes]
    total = total[active_genes]

    callability = called.div(total.astype(float), axis="columns")
    callability_aware = score.mask(missing_by_sample_gene)

    audit_counts["samples"] = len(samples)
    audit_counts["genes"] = len(active_genes)
    audit_counts["mapped_biallelic_variants_observed_in_vcf"] = len(matched_variants)
    audit_counts["mapping_variants_not_observed_in_vcf"] = (
        mapping["Variant"].nunique() - len(matched_variants)
    )
    audit_counts["sample_gene_cells_with_missing_calls"] = int(
        missing_by_sample_gene.to_numpy().sum()
    )

    return (
        score,
        callability_aware,
        callability,
        genotype_counts,
        audit_counts,
    )


def add_target(matrix):
    result = matrix.copy()
    is_control = result.index.str.contains(
        ".mapped.ILLUMINA.bwa.JPT.exome", regex=False
    )
    result["target"] = 1
    result.loc[is_control, "target"] = 0
    return result


def write_outputs(
    score,
    callability_aware,
    callability,
    mapping,
    genotype_counts,
    audit_counts,
):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    paths = {
        "missing_as_zero_matrix": OUTPUT_DIR / "genepy_original_missing_as_zero.csv",
        "missing_as_zero_without_GJA3": OUTPUT_DIR / "genepy_original_missing_as_zero_without_GJA3.csv",
        "callability_aware_matrix": OUTPUT_DIR / "genepy_original_callability_aware.csv",
        "callability_matrix": OUTPUT_DIR / "genepy_callability.csv",
        "variant_audit": OUTPUT_DIR / "genepy_variant_audit.csv",
        "metadata": OUTPUT_DIR / "genepy_build_metadata.json",
    }

    add_target(score).to_csv(paths["missing_as_zero_matrix"])
    final_score = score.drop(
        columns=[gene for gene in TECHNICAL_EXCLUSIONS if gene in score.columns]
    )
    add_target(final_score).to_csv(paths["missing_as_zero_without_GJA3"])
    add_target(callability_aware).to_csv(paths["callability_aware_matrix"])
    callability.to_csv(paths["callability_matrix"])
    mapping.to_csv(paths["variant_audit"], index=False)

    metadata = {
        "formula": "S_gh = -sum_i D_i * log10(f_i1 * f_i2)",
        "heterozygous_frequency_product": "MAF * (1 - MAF)",
        "homozygous_alt_frequency_product": "MAF^2",
        "deleteriousness": {
            "source": "CADD RawScore",
            "scaling": "clipped linear transformation to [0, 1]",
            "lower_bound": CADD_RAW_MIN,
            "upper_bound": CADD_RAW_MAX,
            "protein_truncating_value": 1.0,
        },
        "maf": {
            "priority": [
                "gnomADe_EAS_AF",
                "gnomADe_AF",
                "1000G_EAS_AF",
                "1000G_AF",
                "floor",
            ],
            "floor": MAF_FLOOR,
        },
        "missing_policies": {
            "missing_as_zero_matrix": "Missing GT contributes no observed burden",
            "callability_aware_matrix": "Gene score is NA if any mapped locus has missing/unsupported GT",
        },
        "gene_length_correction": "not applied",
        "genotype_counts": dict(sorted(genotype_counts.items())),
        "audit_counts": dict(sorted(audit_counts.items())),
        "inputs": {
            str(VEP_FILE): file_sha256(VEP_FILE),
            str(CADD_FILE): file_sha256(CADD_FILE),
            str(VCF_FILE): file_sha256(VCF_FILE),
            str(PANEL_FILE): file_sha256(PANEL_FILE),
        },
    }

    metadata["outputs"] = {
        str(path): file_sha256(path)
        for key, path in paths.items()
        if key != "metadata"
    }
    with open(paths["metadata"], "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)

    return paths


def main():
    if not CADD_FILE.exists():
        raise FileNotFoundError(
            f"Missing exact-final-VCF CADD scores: {CADD_FILE}. "
            "Run analysis/scripts/prep/prepare_final_vcf_cadd.py, score the "
            "generated VCF with CADD GRCh37-v1.6, and save the result at this path."
        )
    print("Loading CADD annotations...")
    cadd = load_cadd(CADD_FILE)
    print(f"Loaded {len(cadd):,} unique CADD alleles")

    print("Loading VEP annotations...")
    final_variants = read_biallelic_variant_ids(VCF_FILE)
    panel_genes = load_gene_panel(PANEL_FILE)
    vep = load_vep(
        VEP_FILE,
        wanted_variants=final_variants,
        allowed_genes=panel_genes,
    )
    print(f"Loaded {len(vep):,} unique VEP variant-gene pairs")

    mapping = build_variant_mapping(vep, cadd)
    print(
        f"Matched {mapping['Variant'].nunique():,} variants to "
        f"{mapping['Gene_Name'].nunique():,} genes"
    )

    print("Calculating original GenePy scores from observed GT fields...")
    outputs = process_vcf(VCF_FILE, mapping)
    paths = write_outputs(*outputs[:3], mapping, *outputs[3:])

    print("\nCreated:")
    for path in paths.values():
        print(f"  {path}")


if __name__ == "__main__":
    main()
