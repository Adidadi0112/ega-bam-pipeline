"""Functional interpretation and pathway audit for the final no-GJA3 models.

The analysis deliberately uses the 75 tested matrix features as the custom
enrichment background.  It also audits callability and contributing variants,
because model importance alone is not evidence of biological association.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import gzip
import json
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy.stats import fisher_exact, hypergeom, mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[3]
MATRIX_FILE = ROOT / "analysis/data/genepy_original/genepy_original_missing_as_zero_without_GJA3.csv"
CALLABILITY_FILE = ROOT / "analysis/data/genepy_original/genepy_callability.csv"
VARIANT_AUDIT_FILE = ROOT / "analysis/data/genepy_original/genepy_variant_audit.csv"
VCF_FILE = ROOT / "analysis/02_targeting/uc_vqsr_rare_included.vcf.gz"
VEP_FILE = ROOT / "analysis/05_annotation/vep_results_vqsr.txt"
IMPORTANCE_FILE = ROOT / "analysis/results/final_model_explainability_without_GJA3/heldout_permutation_importance_summary.csv"
SHAP_FILE = ROOT / "analysis/results/final_model_explainability_without_GJA3/shap_summary.csv"
SELECTION_FILE = ROOT / "analysis/results/tpe_model_evaluation_missing_as_zero_without_GJA3/gene_selection_frequency.csv"
OPEN_TARGETS_FILE = ROOT / "OT-EFO_0000729-associated-targets-4_30_2026-v26_03.tsv"
OUTPUT_DIR = ROOT / "analysis/results/functional_interpretation_without_GJA3"

GPROFILER_BASE = "https://biit.cs.ut.ee/gprofiler/api"
PRIMARY_GENES = ["SP140", "GALC", "MYO3B", "ZNF831", "COL28A1"]
# CACNA1C-AS1 is identical to CACNA1C in the matrix and is therefore represented
# once, by the protein-coding CACNA1C symbol, in the sensitivity query.
SENSITIVITY_GENES = PRIMARY_GENES + ["CACNA1C"]
GENE_SETS = {
    "primary_cross_model": PRIMARY_GENES,
    "primary_plus_CACNA1C_locus": SENSITIVITY_GENES,
}
SOURCES = {"GO:BP": "Gene Ontology Biological Process", "REAC": "Reactome"}
MIN_TERM_SIZE = 3
FDR_THRESHOLD = 0.05


FUNCTIONAL_CONTEXT = [
    {
        "Gene": "SP140",
        "Set": "Primary",
        "Molecular_function": "Immune-enriched nuclear chromatin reader and transcriptional regulator.",
        "Mechanistic_theme": "Immune-cell chromatin and macrophage transcription",
        "IBD_evidence_level": "Strong for IBD/Crohn's disease; not UC-specific",
        "Evidence_summary": "Human risk variants alter SP140 expression or splicing, and experimental loss perturbs macrophage programs and worsens colitis; this does not validate a UC biomarker.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/11262",
        "Disease_source": "https://pubmed.ncbi.nlm.nih.gov/28783698/; https://pmc.ncbi.nlm.nih.gov/articles/PMC9442451/",
    },
    {
        "Gene": "GALC",
        "Set": "Primary",
        "Molecular_function": "Lysosomal galactosylceramidase involved in glycosphingolipid catabolism.",
        "Mechanistic_theme": "Lysosomal lipid metabolism",
        "IBD_evidence_level": "Low-to-moderate locus evidence; GALC attribution uncertain",
        "Evidence_summary": "GALC lies in an established IBD locus, but fine-mapping and functional evidence favour neighbouring GPR65 as the effector gene.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/2581",
        "Disease_source": "https://www.nature.com/articles/nature11582; https://pmc.ncbi.nlm.nih.gov/articles/PMC9536022/",
    },
    {
        "Gene": "MYO3B",
        "Set": "Primary",
        "Molecular_function": "Class III actin-activated myosin ATPase with an N-terminal kinase domain.",
        "Mechanistic_theme": "Actin-based motor and protrusion dynamics",
        "IBD_evidence_level": "No direct UC/IBD evidence identified",
        "Evidence_summary": "The cytoskeletal function is biologically plausible, but a direct IBD genetic or mechanistic link was not identified in the focused evidence search.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/140469",
        "Disease_source": "https://www.ebi.ac.uk/gwas/search?query=MYO3B",
    },
    {
        "Gene": "ZNF831",
        "Set": "Primary",
        "Molecular_function": "Poorly characterized zinc-finger protein with predicted zinc-ion binding.",
        "Mechanistic_theme": "Putative transcriptional regulation",
        "IBD_evidence_level": "Moderate IBD-locus/eQTL evidence; mechanism unresolved",
        "Evidence_summary": "An IBD-associated locus has been linked to blood ZNF831 expression, but the causal gene, intestinal mechanism, and UC specificity remain unresolved.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/128611",
        "Disease_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4931595/",
    },
    {
        "Gene": "COL28A1",
        "Set": "Primary",
        "Molecular_function": "Type XXVIII extracellular-matrix collagen containing a von Willebrand factor A domain.",
        "Mechanistic_theme": "Extracellular matrix",
        "IBD_evidence_level": "No direct UC/IBD evidence identified",
        "Evidence_summary": "Its matrix function offers only general plausibility; a barrier or fibrosis role in UC would be speculative without direct evidence.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/340267",
        "Disease_source": "https://www.ebi.ac.uk/gwas/genes/COL28A1",
    },
    {
        "Gene": "CACNA1C",
        "Set": "Sensitivity",
        "Molecular_function": "Pore-forming CaV1.2 L-type voltage-gated calcium-channel subunit.",
        "Mechanistic_theme": "Calcium entry and smooth-muscle function",
        "IBD_evidence_level": "Preclinical colitis-response evidence; no direct human susceptibility evidence",
        "Evidence_summary": "Experimental colitis reduces colonic CACNA1C expression and contractility, supporting a downstream response rather than inherited UC susceptibility.",
        "Function_source": "https://www.ncbi.nlm.nih.gov/gene/775",
        "Disease_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3186840/",
    },
]


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment, retaining NaNs."""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p_values), np.nan)
    valid = np.flatnonzero(np.isfinite(p_values))
    if not len(valid):
        return adjusted
    order = valid[np.argsort(p_values[valid])]
    ranked = p_values[order] * len(valid) / np.arange(1, len(valid) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def normalize_gt(raw: str) -> str:
    gt = str(raw).replace("|", "/")
    if gt in {".", "./.", ""}:
        return "missing"
    if gt == "0/0":
        return "reference"
    if gt in {"0/1", "1/0"}:
        return "heterozygous"
    if gt in {"1/1", "1"}:
        return "homozygous_alt"
    return "unsupported"


def variant_id(chromosome: str, position: str, reference: str, alternative: str) -> str:
    return f"{chromosome.removeprefix('chr')}_{int(position)}_{reference}/{alternative}"


def extra_value(extra: str, key: str) -> str | None:
    for field in str(extra).split(";"):
        if field.startswith(f"{key}="):
            return field.split("=", 1)[1]
    return None


def candidate_vep_context(wanted_variants: set[str]) -> pd.DataFrame:
    """Retrieve alternative frequency fields that were and were not used."""
    rows = []
    header = None
    with VEP_FILE.open("rt") as handle:
        for line in handle:
            if line.startswith("#Uploaded_variation"):
                header = line.rstrip().lstrip("#").split("\t")
                continue
            if header is None or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            record = dict(zip(header, fields))
            current = record["Uploaded_variation"].replace(":", "_").replace("chr", "")
            if current not in wanted_variants:
                continue
            extra = record.get("Extra", "")
            gene = extra_value(extra, "SYMBOL") or record.get("Gene")
            if gene not in SENSITIVITY_GENES:
                continue
            rows.append(
                {
                    "Variant_GRCh37": current,
                    "Gene": gene,
                    "VEP_existing_variation": record.get("Existing_variation"),
                    "VEP_consequence": record.get("Consequence"),
                    "VEP_AF": extra_value(extra, "AF"),
                    "VEP_EAS_AF": extra_value(extra, "EAS_AF"),
                    "VEP_gnomADe_AF": extra_value(extra, "gnomADe_AF"),
                    "VEP_gnomADe_EAS_AF": extra_value(extra, "gnomADe_EAS_AF"),
                    "VEP_clinical_significance": extra_value(extra, "CLIN_SIG"),
                }
            )

    context = pd.DataFrame(rows)
    if context.empty:
        return context

    def first_present(values):
        values = [value for value in values if pd.notna(value) and str(value) not in {"", "-"}]
        return values[0] if values else np.nan

    return (
        context.groupby(["Variant_GRCh37", "Gene"], as_index=False)
        .agg(
            VEP_existing_variation=("VEP_existing_variation", first_present),
            VEP_consequences=("VEP_consequence", lambda values: ";".join(sorted(set(values)))),
            VEP_AF=("VEP_AF", first_present),
            VEP_EAS_AF=("VEP_EAS_AF", first_present),
            VEP_gnomADe_AF=("VEP_gnomADe_AF", first_present),
            VEP_gnomADe_EAS_AF=("VEP_gnomADe_EAS_AF", first_present),
            VEP_clinical_significance=("VEP_clinical_significance", first_present),
        )
    )


def exact_duplicate_audit(matrix: pd.DataFrame) -> pd.DataFrame:
    features = matrix.drop(columns="target")
    remaining = list(features.columns)
    rows = []
    group_number = 0
    while remaining:
        first = remaining.pop(0)
        same = [first] + [gene for gene in remaining if features[first].equals(features[gene])]
        remaining = [gene for gene in remaining if gene not in same]
        if len(same) > 1:
            group_number += 1
            for gene in same:
                rows.append(
                    {
                        "Duplicate_group": group_number,
                        "Gene": gene,
                        "Group_members": ";".join(same),
                        "Group_size": len(same),
                        "Includes_interpretation_gene": gene in SENSITIVITY_GENES,
                    }
                )
    return pd.DataFrame(rows)


def gene_level_audit(matrix: pd.DataFrame, callability: pd.DataFrame) -> pd.DataFrame:
    y = matrix["target"]
    rows = []
    for gene in matrix.columns.drop("target"):
        uc = matrix.loc[y == 1, gene]
        jpt = matrix.loc[y == 0, gene]
        uc_call = callability.loc[y == 1, gene]
        jpt_call = callability.loc[y == 0, gene]
        test = mannwhitneyu(uc, jpt, alternative="two-sided")
        rows.append(
            {
                "Gene": gene,
                "UC_n": len(uc),
                "JPT_n": len(jpt),
                "UC_mean_score": uc.mean(),
                "JPT_mean_score": jpt.mean(),
                "UC_median_score": uc.median(),
                "JPT_median_score": jpt.median(),
                "UC_nonzero_fraction": (uc != 0).mean(),
                "JPT_nonzero_fraction": (jpt != 0).mean(),
                "UC_mean_callability": uc_call.mean(),
                "JPT_mean_callability": jpt_call.mean(),
                "Callability_difference_UC_minus_JPT": uc_call.mean() - jpt_call.mean(),
                "UC_fully_callable_fraction": (uc_call == 1).mean(),
                "JPT_fully_callable_fraction": (jpt_call == 1).mean(),
                "Mann_Whitney_U": test.statistic,
                "Rank_biserial_UC_higher_positive": 2 * test.statistic / (len(uc) * len(jpt)) - 1,
                "Raw_p": test.pvalue,
            }
        )
    result = pd.DataFrame(rows)
    result["BH_q_across_75_genes"] = bh_adjust(result["Raw_p"].to_numpy())
    return result.sort_values(["BH_q_across_75_genes", "Raw_p", "Gene"])


def model_evidence_table(gene_audit: pd.DataFrame) -> pd.DataFrame:
    importance = pd.read_csv(IMPORTANCE_FILE)
    shap = pd.read_csv(SHAP_FILE)
    selection = pd.read_csv(SELECTION_FILE)

    result = gene_audit.set_index("Gene")
    for model in ["Naive Bayes", "Random Forest"]:
        label = "NB" if model == "Naive Bayes" else "RF"
        imp = importance[importance["Model"] == model].set_index("Gene")
        shp = shap[shap["Model"] == model].set_index("Gene")
        sel = selection[selection["Model"] == model].set_index("Gene")
        result[f"{label}_heldout_permutation_importance"] = imp["Mean_heldout_importance"]
        result[f"{label}_positive_fold_fraction"] = imp["Positive_fold_fraction"]
        result[f"{label}_mean_absolute_SHAP"] = shp["Mean_absolute_SHAP"]
        result[f"{label}_selection_frequency"] = sel["Outer_fold_selection_frequency"]

    result["Mean_cross_model_heldout_importance"] = result[
        ["NB_heldout_permutation_importance", "RF_heldout_permutation_importance"]
    ].mean(axis=1)
    result["Positive_importance_in_both_models"] = (
        (result["NB_heldout_permutation_importance"] > 0)
        & (result["RF_heldout_permutation_importance"] > 0)
    )
    return result.reset_index().sort_values(
        "Mean_cross_model_heldout_importance", ascending=False
    )


def observed_variant_context(matrix: pd.DataFrame) -> pd.DataFrame:
    mapping = pd.read_csv(VARIANT_AUDIT_FILE)
    mapping = mapping[mapping["Gene_Name"].isin(SENSITIVITY_GENES)]
    by_variant = {key: group for key, group in mapping.groupby("Variant")}
    target_by_sample = matrix["target"].to_dict()
    rows = []

    with gzip.open(VCF_FILE, "rt") as handle:
        samples = []
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                missing_samples = set(samples).difference(target_by_sample)
                if missing_samples:
                    raise ValueError(f"VCF samples absent from matrix: {sorted(missing_samples)}")
                continue

            fields = line.rstrip().split("\t")
            if "," in fields[4]:
                continue
            current = variant_id(fields[0], fields[1], fields[3], fields[4])
            if current not in by_variant:
                continue
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")
            calls = []
            for sample, sample_field in zip(samples, fields[9:]):
                sample_fields = sample_field.split(":")
                raw_gt = sample_fields[gt_index] if gt_index < len(sample_fields) else "./."
                calls.append((sample, normalize_gt(raw_gt)))

            for annotation in by_variant[current].itertuples(index=False):
                row = {
                    "Gene": annotation.Gene_Name,
                    "Variant_GRCh37": current,
                    "External_ALT_AF_used": annotation.MAF,
                    "Frequency_floor_used": math.isclose(annotation.MAF, 1e-5),
                    "CADD_RawScore": annotation.RawScore,
                    "CADD_PHRED": annotation.PHRED,
                    "Is_PTV": annotation.Is_PTV,
                    "Deleteriousness": annotation.Deleteriousness,
                }
                for target, cohort in [(0, "JPT"), (1, "UC")]:
                    cohort_calls = [gt for sample, gt in calls if target_by_sample[sample] == target]
                    counts = pd.Series(cohort_calls).value_counts()
                    for category in ["reference", "heterozygous", "homozygous_alt", "missing", "unsupported"]:
                        row[f"{cohort}_{category}_n"] = int(counts.get(category, 0))
                    called = (
                        row[f"{cohort}_reference_n"]
                        + row[f"{cohort}_heterozygous_n"]
                        + row[f"{cohort}_homozygous_alt_n"]
                    )
                    alt_count = row[f"{cohort}_heterozygous_n"] + 2 * row[f"{cohort}_homozygous_alt_n"]
                    row[f"{cohort}_observed_ALT_AF"] = alt_count / (2 * called) if called else np.nan
                    row[f"{cohort}_ALT_carrier_fraction_among_called"] = (
                        (row[f"{cohort}_heterozygous_n"] + row[f"{cohort}_homozygous_alt_n"]) / called
                        if called
                        else np.nan
                    )
                rows.append(row)
    result = pd.DataFrame(rows).sort_values(["Gene", "Variant_GRCh37"])
    vep_context = candidate_vep_context(set(result["Variant_GRCh37"]))
    result = result.merge(
        vep_context,
        on=["Variant_GRCh37", "Gene"],
        how="left",
        validate="one_to_one",
    )
    result["Selected_frequency_source"] = np.select(
        [
            result["VEP_gnomADe_EAS_AF"].notna(),
            result["VEP_gnomADe_AF"].notna(),
        ],
        ["gnomADe_EAS_AF", "gnomADe_AF"],
        default="1e-5 floor",
    )
    vep_eas = pd.to_numeric(result["VEP_EAS_AF"], errors="coerce")
    vep_global = pd.to_numeric(result["VEP_AF"], errors="coerce")
    observed_max = result[["UC_observed_ALT_AF", "JPT_observed_ALT_AF"]].max(axis=1)
    result["Floor_conflicts_with_common_observed_or_VEP_AF"] = (
        result["Frequency_floor_used"]
        & (vep_eas.fillna(vep_global).fillna(0).ge(0.01) | observed_max.ge(0.01))
    )
    return result


def request_json(method: str, url: str, **kwargs) -> dict:
    headers = {"User-Agent": "GenePy-master-thesis-functional-audit/1.0"}
    response = requests.request(method, url, headers=headers, timeout=120, **kwargs)
    response.raise_for_status()
    return response.json()


def convert_background(symbols: list[str]) -> tuple[pd.DataFrame, dict]:
    raw = request_json(
        "POST",
        f"{GPROFILER_BASE}/convert/convert/",
        json={"organism": "hsapiens", "query": symbols, "target": "ENSG"},
    )
    converted = pd.DataFrame(raw["result"])
    converted["Mapped"] = converted["converted"].notna() & converted["converted"].ne("None")
    return converted, raw


def run_enrichment(
    symbol_to_ensg: dict[str, str], background_ensg: list[str]
) -> tuple[pd.DataFrame, dict]:
    rows = []
    raw_responses = {}
    ensg_to_symbol = {ensg: symbol for symbol, ensg in symbol_to_ensg.items()}

    for set_name, symbols in GENE_SETS.items():
        query_ensg = [symbol_to_ensg[symbol] for symbol in symbols if symbol in symbol_to_ensg]
        for source, source_name in SOURCES.items():
            payload = {
                "organism": "hsapiens",
                "query": query_ensg,
                "sources": [source],
                "user_threshold": FDR_THRESHOLD,
                "all_results": True,
                "ordered": False,
                "domain_scope": "custom",
                "background": background_ensg,
                "significance_threshold_method": "fdr",
                "no_evidences": False,
            }
            raw = request_json(
                "POST", f"{GPROFILER_BASE}/gost/profile/", json=payload
            )
            raw_responses[f"{set_name}__{source}"] = {"payload": payload, "response": raw}

            for result in raw.get("result", []):
                overlap_ids = [
                    gene_id
                    for gene_id, evidence in zip(query_ensg, result.get("intersections", []))
                    if evidence
                ]
                overlap_symbols = [ensg_to_symbol.get(gene_id, gene_id) for gene_id in overlap_ids]
                n = int(result["effective_domain_size"])
                k_term = int(result["term_size"])
                n_query = int(result["query_size"])
                overlap = int(result["intersection_size"])
                table = [[overlap, n_query - overlap], [k_term - overlap, n - k_term - n_query + overlap]]
                odds_ratio = fisher_exact(table, alternative="greater").statistic
                rows.append(
                    {
                        "Gene_set": set_name,
                        "Source": source,
                        "Source_name": source_name,
                        "Term_ID": result["native"],
                        "Term_name": result["name"],
                        "Effective_background_size": n,
                        "Mapped_query_size": n_query,
                        "Background_term_size": k_term,
                        "Overlap_size": overlap,
                        "Overlap_genes": ";".join(overlap_symbols),
                        "Odds_ratio": odds_ratio,
                        "Calculated_raw_hypergeometric_p": hypergeom.sf(overlap - 1, n, k_term, n_query),
                        "FDR_adjusted_p_from_gProfiler": float(result["p_value"]),
                        "Meets_minimum_term_size": k_term >= MIN_TERM_SIZE,
                        "Significant_FDR_0.05": (
                            k_term >= MIN_TERM_SIZE and float(result["p_value"]) < FDR_THRESHOLD
                        ),
                    }
                )
    return pd.DataFrame(rows), raw_responses


def add_open_targets(functional: pd.DataFrame) -> pd.DataFrame:
    targets = pd.read_csv(OPEN_TARGETS_FILE, sep="\t", dtype=str)
    wanted = ["symbol", "globalScore", "gwasCredibleSets", "europepmc", "expressionAtlas"]
    targets = targets[wanted].rename(
        columns={
            "symbol": "Gene",
            "globalScore": "Open_Targets_global_score",
            "gwasCredibleSets": "Open_Targets_GWAS_credible_sets_score",
            "europepmc": "Open_Targets_Europe_PMC_score",
            "expressionAtlas": "Open_Targets_expression_score",
        }
    )
    result = functional.merge(targets, on="Gene", how="left", validate="one_to_one")
    result["Open_Targets_interpretation"] = (
        "Panel-source evidence is descriptive, not independent validation, because the model background was selected from Open Targets."
    )
    return result


def plot_candidate_distributions(matrix: pd.DataFrame) -> None:
    plot_data = (
        matrix[SENSITIVITY_GENES + ["target"]]
        .rename(columns={"target": "Cohort"})
        .assign(Cohort=lambda frame: frame["Cohort"].map({0: "JPT reference", 1: "UC"}))
        .melt(id_vars="Cohort", var_name="Gene", value_name="GenePy score")
    )
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
    for gene, axis in zip(SENSITIVITY_GENES, axes.flat):
        subset = plot_data[plot_data["Gene"] == gene]
        sns.boxplot(data=subset, x="Cohort", y="GenePy score", hue="Cohort", ax=axis, showfliers=False, legend=False)
        sns.stripplot(data=subset, x="Cohort", y="GenePy score", ax=axis, color="black", alpha=0.32, size=2.5, jitter=0.22)
        axis.set_title(gene, fontstyle="italic")
        axis.set_xlabel("")
        axis.tick_params(axis="x", rotation=12)
    fig.suptitle("Candidate GenePy distributions (descriptive; post-selection)", y=1.01)
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(OUTPUT_DIR / f"candidate_genepy_distributions.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_callability(candidate_audit: pd.DataFrame) -> None:
    plot_data = candidate_audit[
        ["Gene", "UC_fully_callable_fraction", "JPT_fully_callable_fraction"]
    ].melt(id_vars="Gene", var_name="Cohort", value_name="Fully callable fraction")
    plot_data["Cohort"] = plot_data["Cohort"].map(
        {"UC_fully_callable_fraction": "UC", "JPT_fully_callable_fraction": "JPT reference"}
    )
    fig, axis = plt.subplots(figsize=(8, 4.8))
    sns.barplot(data=plot_data, y="Gene", x="Fully callable fraction", hue="Cohort", ax=axis)
    axis.set_xlim(0, 1.04)
    axis.set_xlabel("Fraction of samples with all contributing loci called")
    axis.set_ylabel("")
    axis.set_title("Candidate-gene callability differs by source cohort")
    axis.legend(title="")
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(OUTPUT_DIR / f"candidate_callability_audit.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_importance_vs_callability(evidence: pd.DataFrame) -> float:
    rho, _ = spearmanr(
        evidence["Callability_difference_UC_minus_JPT"],
        evidence["Mean_cross_model_heldout_importance"],
    )
    fig, axis = plt.subplots(figsize=(8, 5.5))
    highlighted = evidence["Gene"].isin(SENSITIVITY_GENES)
    axis.scatter(
        evidence.loc[~highlighted, "Callability_difference_UC_minus_JPT"],
        evidence.loc[~highlighted, "Mean_cross_model_heldout_importance"],
        color="#8d99ae",
        alpha=0.65,
        label="Other final features",
    )
    axis.scatter(
        evidence.loc[highlighted, "Callability_difference_UC_minus_JPT"],
        evidence.loc[highlighted, "Mean_cross_model_heldout_importance"],
        color="#c1121f",
        s=52,
        label="Interpreted genes",
    )
    for row in evidence.loc[highlighted].itertuples(index=False):
        axis.annotate(
            row.Gene,
            (row.Callability_difference_UC_minus_JPT, row.Mean_cross_model_heldout_importance),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    axis.axvline(0, color="black", linewidth=0.8, linestyle="--")
    axis.axhline(0, color="black", linewidth=0.8, linestyle=":")
    axis.set_xlabel("Mean callability difference (UC - JPT reference)")
    axis.set_ylabel("Mean held-out importance across NB and RF")
    axis.set_title(f"Feature importance is associated with callability imbalance (Spearman rho={rho:.2f})")
    axis.legend(frameon=False)
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(OUTPUT_DIR / f"importance_vs_callability.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return float(rho)


def plot_enrichment_null_or_hits(enrichment: pd.DataFrame) -> None:
    significant = enrichment[enrichment["Significant_FDR_0.05"]].copy()
    fig, axis = plt.subplots(figsize=(9, 4.8))
    if significant.empty:
        tested = enrichment.groupby(["Gene_set", "Source"]).size()
        lines = ["No GO Biological Process or Reactome term reached FDR q < 0.05.", ""]
        lines.extend(f"{gene_set}, {source}: {count} overlapping terms returned" for (gene_set, source), count in tested.items())
        axis.text(0.5, 0.5, "\n".join(lines), ha="center", va="center", fontsize=12)
        axis.set_axis_off()
        axis.set_title("Custom-background over-representation analysis: null result")
    else:
        significant = significant.nsmallest(12, "FDR_adjusted_p_from_gProfiler").sort_values(
            "FDR_adjusted_p_from_gProfiler", ascending=False
        )
        significant["minus_log10_q"] = -np.log10(significant["FDR_adjusted_p_from_gProfiler"])
        sns.barplot(data=significant, y="Term_name", x="minus_log10_q", hue="Source", ax=axis)
        axis.axvline(-np.log10(FDR_THRESHOLD), color="black", linestyle="--", linewidth=0.9)
        axis.set_xlabel("-log10(FDR-adjusted p)")
        axis.set_ylabel("")
        axis.set_title("Significant custom-background pathway terms")
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(OUTPUT_DIR / f"pathway_enrichment_summary.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    matrix = pd.read_csv(MATRIX_FILE, index_col=0)
    callability = pd.read_csv(CALLABILITY_FILE, index_col=0).loc[matrix.index]
    background_symbols = matrix.columns.drop("target").tolist()

    duplicates = exact_duplicate_audit(matrix)
    duplicates.to_csv(OUTPUT_DIR / "exact_duplicate_feature_groups.csv", index=False)

    gene_audit = gene_level_audit(matrix, callability)
    gene_audit.to_csv(OUTPUT_DIR / "gene_level_group_and_callability_audit.csv", index=False)

    evidence = model_evidence_table(gene_audit)
    evidence.to_csv(OUTPUT_DIR / "model_evidence_and_callability_all_75.csv", index=False)
    candidate_evidence = evidence[evidence["Gene"].isin(SENSITIVITY_GENES)].copy()
    candidate_evidence["Interpretation_set"] = candidate_evidence["Gene"].map(
        {gene: "Primary" for gene in PRIMARY_GENES} | {"CACNA1C": "Sensitivity"}
    )
    candidate_evidence.to_csv(OUTPUT_DIR / "candidate_gene_model_and_technical_evidence.csv", index=False)

    variants = observed_variant_context(matrix)
    variants.to_csv(OUTPUT_DIR / "candidate_observed_variant_context.csv", index=False)

    functional = add_open_targets(pd.DataFrame(FUNCTIONAL_CONTEXT))
    functional.to_csv(OUTPUT_DIR / "candidate_function_and_disease_evidence.csv", index=False)

    conversion, raw_conversion = convert_background(background_symbols)
    conversion.to_csv(OUTPUT_DIR / "gprofiler_identifier_mapping_audit.csv", index=False)
    with (OUTPUT_DIR / "gprofiler_identifier_mapping_raw.json").open("w") as handle:
        json.dump(raw_conversion, handle, indent=2)

    mapped = conversion[conversion["Mapped"]].drop_duplicates("incoming")
    symbol_to_ensg = mapped.set_index("incoming")["converted"].to_dict()
    missing_candidates = set(SENSITIVITY_GENES).difference(symbol_to_ensg)
    if missing_candidates:
        raise ValueError(f"Candidate symbols not mapped by g:Profiler: {sorted(missing_candidates)}")
    background_ensg = sorted(set(mapped["converted"]))

    enrichment, raw_enrichment = run_enrichment(symbol_to_ensg, background_ensg)
    enrichment = enrichment.sort_values(
        ["Gene_set", "Source", "FDR_adjusted_p_from_gProfiler", "Term_ID"]
    )
    enrichment.to_csv(OUTPUT_DIR / "pathway_enrichment_all_terms.csv", index=False)
    enrichment[enrichment["Significant_FDR_0.05"]].to_csv(
        OUTPUT_DIR / "pathway_enrichment_significant_terms.csv", index=False
    )
    with (OUTPUT_DIR / "gprofiler_enrichment_raw.json").open("w") as handle:
        json.dump(raw_enrichment, handle, indent=2)

    data_versions = request_json(
        "GET", f"{GPROFILER_BASE}/util/data_versions?organism=hsapiens"
    )
    with (OUTPUT_DIR / "gprofiler_data_versions.json").open("w") as handle:
        json.dump(data_versions, handle, indent=2)

    plot_candidate_distributions(matrix)
    plot_callability(candidate_evidence)
    importance_callability_rho = plot_importance_vs_callability(evidence)
    plot_enrichment_null_or_hits(enrichment)

    metadata = {
        "analysis_scope": "exploratory functional interpretation after final model fitting",
        "primary_query": PRIMARY_GENES,
        "sensitivity_query": SENSITIVITY_GENES,
        "custom_background_matrix_features": len(background_symbols),
        "custom_background_successfully_mapped_unique_ensembl_genes": len(background_ensg),
        "exact_duplicate_groups": int(duplicates["Duplicate_group"].nunique()),
        "unique_feature_vectors": len(background_symbols) - len(duplicates) + int(duplicates["Duplicate_group"].nunique()),
        "enrichment_sources": SOURCES,
        "minimum_background_term_size": MIN_TERM_SIZE,
        "fdr_threshold": FDR_THRESHOLD,
        "multiple_testing": "g:Profiler Benjamini-Hochberg FDR, queried separately for each source and gene set",
        "spearman_importance_vs_callability_difference": importance_callability_rho,
        "significant_enrichment_terms": int(enrichment["Significant_FDR_0.05"].sum()),
        "important_interpretation_note": (
            "The candidate genes are post-selection explanations, and several are strongly coupled to differential callability. "
            "Functional annotation and enrichment do not validate association, mechanism, or biomarker status."
        ),
        "input_sha256": {
            str(path.relative_to(ROOT)): file_sha256(path)
            for path in [MATRIX_FILE, CALLABILITY_FILE, VARIANT_AUDIT_FILE, VCF_FILE, VEP_FILE, IMPORTANCE_FILE, SHAP_FILE]
        },
    }
    with (OUTPUT_DIR / "analysis_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)

    print(f"Created functional audit in {OUTPUT_DIR.relative_to(ROOT)}")
    print(f"Mapped background: {len(background_ensg)}/{len(background_symbols)} unique genes")
    print(f"Exact duplicate groups: {duplicates['Duplicate_group'].nunique()} (56 unique feature vectors)")
    print(f"Significant GO/Reactome terms at FDR < 0.05: {metadata['significant_enrichment_terms']}")
    print(f"Importance-callability Spearman rho: {importance_callability_rho:.3f}")


if __name__ == "__main__":
    main()
