"""Pathway-level GenePy features and rank-based enrichment for Wave 1."""

from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from analysis.scripts.paper.wave1_identifiability import (  # noqa: E402
    OUTPUT_DIR,
    RANDOM_STATE,
    evaluate_matrix,
    load_labelled_matrix,
    summarise,
)
from sklearn.model_selection import RepeatedStratifiedKFold

DATA_DIR = ROOT / "analysis/data/genepy_expanded"
CONFIG_DIR = ROOT / "paper/config"
MIN_SET_SIZE = 3
ENRICHR_LIBRARIES = {
    "Reactome_2022": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022",
    "GO_Biological_Process_2023": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023",
    "KEGG_2021_Human": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
}


def download_gmt(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    with urllib.request.urlopen(url, timeout=60) as response:
        destination.write_bytes(response.read())


def parse_gmt(path, allowed_genes):
    sets = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            members = [gene for gene in parts[2:] if gene in allowed_genes]
            if len(members) >= MIN_SET_SIZE:
                sets[name] = sorted(set(members))
    return sets


def pathway_matrix(genepy, gene_sets, prefix):
    rows = {}
    for name, members in gene_sets.items():
        column = f"{prefix}__{name}"
        rows[column] = genepy[members].mean(axis=1)
    return pd.DataFrame(rows, index=genepy.index)


def random_gene_set_matrix(genepy, sizes, n_sets, seed):
    rng = np.random.default_rng(seed)
    genes = list(genepy.columns)
    columns = {}
    manifest = []
    for set_index, size in enumerate(sizes):
        for draw in range(n_sets):
            members = sorted(rng.choice(genes, size=size, replace=False).tolist())
            name = f"random_size{size}_draw{draw}"
            columns[name] = genepy[members].mean(axis=1)
            manifest.append({"Set": name, "Size": size, "Genes": ",".join(members)})
    return pd.DataFrame(columns, index=genepy.index), pd.DataFrame(manifest)


def prerank_enrichment(ranked_genes, gene_sets, n_permutations=1000, seed=42):
    """Simple weighted KS-style prerank NES with label-permutation FDR."""
    rank_map = {gene: i for i, gene in enumerate(ranked_genes)}
    n_genes = len(ranked_genes)
    scores = np.linspace(1, -1, n_genes)

    def enrichment_score(members):
        hits = np.array([rank_map[gene] for gene in members if gene in rank_map])
        if hits.size < MIN_SET_SIZE:
            return np.nan
        running = np.zeros(n_genes)
        hit_weight = np.zeros(n_genes)
        hit_weight[hits] = np.abs(scores[hits])
        hit_total = hit_weight.sum()
        miss_total = n_genes - hits.size
        if hit_total == 0 or miss_total == 0:
            return np.nan
        increment = hit_weight / hit_total
        decrement = np.full(n_genes, 1 / miss_total)
        decrement[hits] = 0
        running = np.cumsum(increment - decrement)
        return float(running[np.argmax(np.abs(running))])

    observed = []
    for name, members in gene_sets.items():
        observed.append({
            "Set": name,
            "Size": len([gene for gene in members if gene in rank_map]),
            "ES": enrichment_score(members),
        })
    observed = pd.DataFrame(observed).dropna(subset=["ES"])

    rng = np.random.default_rng(seed)
    null_abs = []
    genes = np.array(ranked_genes)
    for _ in range(n_permutations):
        shuffled = genes.copy()
        rng.shuffle(shuffled)
        rank_map = {gene: i for i, gene in enumerate(shuffled)}
        scores = np.linspace(1, -1, n_genes)
        for members in gene_sets.values():
            value = enrichment_score(members)
            if np.isfinite(value):
                null_abs.append(abs(value))
    null_abs = np.array(null_abs) if null_abs else np.array([1.0])
    observed["NES"] = observed["ES"] / (null_abs.mean() if null_abs.mean() else 1)
    observed["Permutation_p"] = [
        (1 + (null_abs >= abs(es)).sum()) / (1 + len(null_abs))
        for es in observed["ES"]
    ]
    observed = observed.sort_values("Permutation_p")
    n_tests = len(observed)
    observed["BH_FDR"] = np.minimum(
        1.0,
        observed["Permutation_p"].to_numpy() * n_tests / np.arange(1, n_tests + 1),
    )
    observed["BH_FDR"] = np.minimum.accumulate(observed["BH_FDR"].iloc[::-1])[::-1]
    return observed


def gprofiler_ordered(ranked_genes, background):
    try:
        from gprofiler import GProfiler
    except ImportError:
        return pd.DataFrame()
    profiler = GProfiler(return_dataframe=True)
    result = profiler.profile(
        organism="hsapiens",
        query=ranked_genes,
        background=background,
        domain_scope="custom",
        ordered=True,
        sources=["REAC", "GO:BP", "KEGG"],
        significance_threshold_method="fdr",
        user_threshold=0.05,
        all_results=True,
        no_evidences=False,
    )
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)

    X, y = load_labelled_matrix(
        DATA_DIR / "genepy_original_missing_as_zero_without_GJA3.csv"
    )
    allowed = set(X.columns)

    all_sets = {}
    set_sizes = []
    for library, url in ENRICHR_LIBRARIES.items():
        gmt_path = CONFIG_DIR / f"{library}.gmt"
        print(f"Loading {library}...")
        try:
            download_gmt(url, gmt_path)
        except Exception as error:
            print(f"  skipped download ({error})")
            if not gmt_path.exists():
                continue
        parsed = parse_gmt(gmt_path, allowed)
        print(f"  {len(parsed)} sets with ≥{MIN_SET_SIZE} modelled genes")
        prefix = library.split("_")[0]
        for name, members in parsed.items():
            all_sets[f"{prefix}__{name}"] = members
            set_sizes.append(len(members))
        membership_rows = [
            {"Library": library, "Set": name, "Size": len(members), "Genes": ",".join(members)}
            for name, members in parsed.items()
        ]
        pd.DataFrame(membership_rows).to_csv(
            OUTPUT_DIR / f"gene_sets_{library}.csv", index=False
        )

    if not all_sets:
        raise SystemExit("No gene sets could be constructed.")

    X_path = pathway_matrix(X, all_sets, prefix="pathway")
    X_path.to_csv(OUTPUT_DIR / "pathway_feature_matrix.csv")

    size_values = [int(value) for value in np.quantile(set_sizes, [0.25, 0.5, 0.75])]
    X_random, random_manifest = random_gene_set_matrix(
        X, size_values, n_sets=8, seed=RANDOM_STATE
    )
    random_manifest.to_csv(OUTPUT_DIR / "random_gene_set_manifest.csv", index=False)
    X_random.to_csv(OUTPUT_DIR / "random_gene_set_feature_matrix.csv")

    outer_cv = RepeatedStratifiedKFold(
        n_splits=5, n_repeats=10, random_state=RANDOM_STATE
    )
    outer_splits = list(outer_cv.split(X, y))
    performances = []
    predictions = []
    for name, matrix in [
        ("pathway_mean_sets", X_path),
        ("random_matched_gene_sets", X_random),
    ]:
        performance, prediction = evaluate_matrix(name, matrix, y, outer_splits)
        performances.append(performance)
        predictions.append(prediction)
    performance = pd.concat(performances, ignore_index=True)
    prediction = pd.concat(predictions, ignore_index=True)
    repetition_means, summary = summarise(performance)
    performance.to_csv(OUTPUT_DIR / "pathway_outer_fold_performance.csv", index=False)
    prediction.to_csv(OUTPUT_DIR / "pathway_out_of_fold_predictions.csv", index=False)
    repetition_means.to_csv(
        OUTPUT_DIR / "pathway_repetition_mean_performance.csv", index=False
    )
    summary.to_csv(OUTPUT_DIR / "pathway_representation_summary.csv", index=False)

    rank_file = OUTPUT_DIR / "univariate_gene_ranks.csv"
    if rank_file.exists():
        ranks = pd.read_csv(rank_file)
    else:
        rows = []
        for gene in X.columns:
            statistic, p_value = mannwhitneyu(
                X.loc[y == 1, gene], X.loc[y == 0, gene], alternative="two-sided"
            )
            rows.append({"Gene": gene, "MannWhitney_p": p_value, "Signed_AUROC": roc_auc_score(y, X[gene])})
        ranks = pd.DataFrame(rows)
    ranked = ranks.sort_values("MannWhitney_p")["Gene"].tolist()
    residual_col = "Callability_residualized_mean_difference"
    if residual_col in ranks.columns:
        residual_ranked = (
            ranks.dropna(subset=[residual_col])
            .assign(abs_resid=lambda frame: frame[residual_col].abs())
            .sort_values("abs_resid", ascending=False)["Gene"]
            .tolist()
        )
    else:
        residual_ranked = ranked

    prerank = prerank_enrichment(ranked, all_sets)
    prerank.to_csv(OUTPUT_DIR / "prerank_enrichment_univariate.csv", index=False)
    residual_prerank = prerank_enrichment(residual_ranked, all_sets)
    residual_prerank.to_csv(
        OUTPUT_DIR / "prerank_enrichment_callability_residualized.csv", index=False
    )

    shap_path = ROOT / (
        "analysis/results/final_model_explainability_expanded_without_GJA3/"
        "shap_summary.csv"
    )
    perm_path = ROOT / (
        "analysis/results/final_model_explainability_expanded_without_GJA3/"
        "heldout_permutation_importance_summary.csv"
    )
    if shap_path.exists():
        shap = pd.read_csv(shap_path)
        if "Model" in shap.columns:
            shap = shap[shap["Model"] == "Random Forest"]
        gene_col = "Gene" if "Gene" in shap.columns else shap.columns[0]
        score_col = (
            "Mean_absolute_SHAP"
            if "Mean_absolute_SHAP" in shap.columns
            else shap.columns[1]
        )
        shap_ranked = shap.sort_values(score_col, ascending=False)[gene_col].tolist()
        shap_prerank = prerank_enrichment(shap_ranked, all_sets)
        shap_prerank.to_csv(OUTPUT_DIR / "prerank_enrichment_shap.csv", index=False)
    if perm_path.exists():
        perm = pd.read_csv(perm_path)
        perm = perm[perm["Model"] == "Random Forest"] if "Model" in perm.columns else perm
        perm_ranked = perm.sort_values(
            "Mean_heldout_importance", ascending=False
        )["Gene"].tolist()
        perm_prerank = prerank_enrichment(perm_ranked, all_sets)
        perm_prerank.to_csv(
            OUTPUT_DIR / "prerank_enrichment_permutation_importance.csv", index=False
        )

    ordered = gprofiler_ordered(ranked, list(X.columns))
    if not ordered.empty:
        ordered.to_csv(OUTPUT_DIR / "gprofiler_ordered_univariate.csv", index=False)
        significant = ordered[
            ordered.get("significant", False)
            & (ordered.get("intersection_size", 0) >= 2)
        ] if "significant" in ordered.columns else pd.DataFrame()
        significant.to_csv(
            OUTPUT_DIR / "gprofiler_ordered_univariate_significant.csv", index=False
        )

    metadata = {
        "min_set_size": MIN_SET_SIZE,
        "n_pathway_features": int(X_path.shape[1]),
        "n_random_gene_set_features": int(X_random.shape[1]),
        "libraries": list(ENRICHR_LIBRARIES),
        "interpretation": (
            "Exploratory. Custom 215-gene UC panel background. Not biomarker discovery."
        ),
    }
    with (OUTPUT_DIR / "pathway_analysis_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)
    print("\nPathway representation summary")
    print(summary.round(4).to_string(index=False))
    print(f"Wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
