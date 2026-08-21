"""ORA of stable Random Forest genes from the expanded no-GJA3 matrix."""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from gprofiler import GProfiler


ROOT = Path(__file__).resolve().parents[3]
IMPORTANCE_FILE = ROOT / (
    "analysis/results/final_model_explainability_expanded_without_GJA3/"
    "heldout_permutation_importance_summary.csv"
)
OUTPUT_DIR = ROOT / "analysis/results/rf_pathway_enrichment_expanded"
MAX_GENES = 20
MIN_POSITIVE_FOLD_FRACTION = 0.60
MIN_INTERSECTION_SIZE = 2
FDR_THRESHOLD = 0.05
SOURCES = ["REAC", "GO:BP"]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    importance = pd.read_csv(IMPORTANCE_FILE)
    rf = (
        importance[importance["Model"] == "Random Forest"]
        .sort_values("Mean_heldout_importance", ascending=False)
    )
    background = rf["Gene"].dropna().drop_duplicates().tolist()
    selected_frame = rf[
        (rf["Mean_heldout_importance"] > 0)
        & (rf["Positive_fold_fraction"] >= MIN_POSITIVE_FOLD_FRACTION)
    ].head(MAX_GENES)
    selected = selected_frame["Gene"].tolist()

    selected_frame.to_csv(OUTPUT_DIR / "rf_pathway_selected_genes.csv", index=False)
    pd.DataFrame({"Gene": background}).to_csv(
        OUTPUT_DIR / "rf_pathway_background_genes.csv", index=False
    )

    profiler = GProfiler(return_dataframe=True)
    results = profiler.profile(
        organism="hsapiens",
        query=selected,
        background=background,
        domain_scope="custom",
        sources=SOURCES,
        significance_threshold_method="fdr",
        user_threshold=FDR_THRESHOLD,
        all_results=True,
        no_evidences=False,
    )
    results.to_csv(OUTPUT_DIR / "rf_pathway_enrichment_all_terms.csv", index=False)

    if results.empty:
        significant = results.copy()
        summary = results.copy()
    else:
        significant = results[
            results["significant"]
            & (results["intersection_size"] >= MIN_INTERSECTION_SIZE)
        ].sort_values("p_value")
        summary_columns = [
            "source",
            "native",
            "name",
            "p_value",
            "intersection_size",
            "term_size",
            "intersections",
        ]
        summary = significant[summary_columns].head(10)

    significant.to_csv(
        OUTPUT_DIR / "rf_pathway_enrichment_significant_terms.csv", index=False
    )
    summary.to_csv(
        OUTPUT_DIR / "rf_pathway_enrichment_summary.csv", index=False
    )

    figure, axis = plt.subplots(figsize=(9, 5))
    if summary.empty:
        axis.text(
            0.5,
            0.5,
            "No GO Biological Process or Reactome pathway\n"
            "remained significant after FDR correction.",
            ha="center",
            va="center",
            fontsize=12,
        )
        axis.set_axis_off()
        axis.set_title("Random Forest pathway over-representation analysis")
    else:
        plot_data = summary.copy().sort_values("p_value", ascending=False)
        plot_data["minus_log10_fdr"] = -np.log10(plot_data["p_value"])
        sns.barplot(
            data=plot_data,
            y="name",
            x="minus_log10_fdr",
            hue="source",
            ax=axis,
        )
        axis.axvline(-np.log10(FDR_THRESHOLD), color="black", linestyle="--")
        axis.set_xlabel("-log10(FDR-adjusted p-value)")
        axis.set_ylabel("")
        axis.set_title("Random Forest pathway over-representation analysis")
    figure.tight_layout()
    for extension in ["svg", "png"]:
        figure.savefig(
            OUTPUT_DIR / f"rf_pathway_enrichment.{extension}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)

    metadata = {
        "analysis": "over-representation analysis",
        "model": "Random Forest",
        "selection": (
            "Mean held-out permutation importance > 0 and positive "
            f"fold fraction >= {MIN_POSITIVE_FOLD_FRACTION}; top {MAX_GENES}"
        ),
        "selected_gene_count": len(selected),
        "selected_genes": selected,
        "custom_background_gene_count": len(background),
        "sources": SOURCES,
        "multiple_testing": "g:Profiler FDR",
        "fdr_threshold": FDR_THRESHOLD,
        "minimum_reported_intersection_size": MIN_INTERSECTION_SIZE,
        "significant_reported_term_count": len(significant),
        "interpretation": (
            "Exploratory post-model functional interpretation conditional on "
            "the 215-gene Open Targets-derived modelling background."
        ),
    }
    with (OUTPUT_DIR / "analysis_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Selected genes ({len(selected)}): {', '.join(selected)}")
    print(f"Custom background: {len(background)} genes")
    print(f"Reported significant terms: {len(significant)}")
    if summary.empty:
        print("No term passed FDR < 0.05 with at least two selected genes.")
    else:
        print(summary.to_string(index=False))
    print(f"Results written to {OUTPUT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
