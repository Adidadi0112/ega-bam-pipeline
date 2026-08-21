"""Prepare and validate exact-final-VCF input for local/manual CADD scoring."""

from pathlib import Path
import gzip

import pandas as pd


VCF_FILE = Path("analysis/02_targeting/uc_vqsr_rare_included.vcf.gz")
CADD_INPUT = Path("analysis/05_annotation/final_biallelic_snvs_for_cadd.vcf")
CADD_OUTPUT = Path("analysis/05_annotation/cadd_scores_final_vcf.gz")


def variant_id(chromosome, position, reference, alternative):
    return f"{str(chromosome).removeprefix('chr')}_{int(position)}_{reference}/{alternative}"


def open_text_auto(path):
    with open(path, "rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else open(path, "rt")


def final_variants():
    rows = []
    with gzip.open(VCF_FILE, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if "," in fields[4]:
                continue
            rows.append((fields[0].removeprefix("chr"), int(fields[1]), fields[3], fields[4]))
    return sorted(set(rows), key=lambda row: (str(row[0]), row[1], row[2], row[3]))


def write_cadd_input(variants):
    CADD_INPUT.parent.mkdir(parents=True, exist_ok=True)
    with CADD_INPUT.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh37\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chromosome, position, reference, alternative in variants:
            handle.write(
                f"{chromosome}\t{position}\t.\t{reference}\t{alternative}\t.\tPASS\t.\n"
            )


def validate_cadd_output(variants):
    if not CADD_OUTPUT.exists():
        return False

    header_line = None
    with open_text_auto(CADD_OUTPUT) as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("#Chrom") or line.startswith("Chrom\t"):
                header_line = line_number
                break
    if header_line is None:
        raise ValueError(f"No CADD header found in {CADD_OUTPUT}")

    with open_text_auto(CADD_OUTPUT) as handle:
        scores = pd.read_csv(handle, sep="\t", skiprows=header_line)
    scores.columns = [column.lstrip("#") for column in scores.columns]
    required = {"Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"}
    if not required.issubset(scores.columns):
        raise ValueError(f"CADD output is missing columns: {sorted(required - set(scores.columns))}")

    expected = {variant_id(*row) for row in variants}
    observed = {
        variant_id(row.Chrom, row.Pos, row.Ref, row.Alt)
        for row in scores.itertuples(index=False)
    }
    missing = expected - observed
    extra = observed - expected
    print(f"CADD output coverage: {len(expected - missing)}/{len(expected)}")
    if missing:
        print(f"Missing CADD variants: {len(missing)}")
    if extra:
        print(f"Extra CADD variants: {len(extra)}")
    return not missing


def main():
    variants = final_variants()
    write_cadd_input(variants)
    print(f"Created {CADD_INPUT} with {len(variants)} biallelic SNVs and no sample columns")
    if validate_cadd_output(variants):
        print("CADD output is complete; create_genepy_matrix.py can now be run.")
    else:
        print(f"Score this file with CADD GRCh37-v1.6 and save the result as {CADD_OUTPUT}")


if __name__ == "__main__":
    main()
