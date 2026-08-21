"""Publication figures for Wave 1 representation comparisons."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = ROOT / "paper/results/wave1"
FIGURE_DIR = OUTPUT_DIR / "figures"


def save(figure, stem):
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ["svg", "png"]:
        figure.savefig(
            FIGURE_DIR / f"{stem}.{extension}", dpi=300, bbox_inches="tight"
        )
    plt.close(figure)


def main():
    summary = pd.read_csv(OUTPUT_DIR / "representation_summary.csv")
    repetitions = pd.read_csv(OUTPUT_DIR / "repetition_mean_performance.csv")
    labels = {
        "missing_as_zero_without_GJA3": "UC panel (missing-as-zero)",
        "missing_as_zero_with_GJA3": "UC panel + GJA3",
        "callability_aware_without_GJA3": "Callability-aware",
        "qc_source_baseline": "QC / source features",
    }
    primary = [name for name in labels if name in set(summary.Representation)]
    random_reps = [
        name for name in summary.Representation
        if name.startswith("random_30")
    ]
    plot = repetitions[repetitions.Representation.isin(primary)].copy()
    plot["Label"] = plot["Representation"].map(labels)
    order = (
        summary[summary.Representation.isin(primary)]
        .sort_values("AUROC_mean", ascending=False)["Representation"]
        .map(labels)
        .tolist()
    )
    plt.figure(figsize=(8.5, 5))
    sns.boxplot(data=plot, x="Label", y="Outer_AUROC", order=order, color="white")
    sns.stripplot(
        data=plot, x="Label", y="Outer_AUROC", order=order, color="#1f4e79", size=5
    )
    plt.axhline(0.5, color="grey", linestyle="--", linewidth=1)
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("Repetition-level AUROC")
    plt.xlabel("")
    plt.title("Wave 1: technical vs gene-panel discrimination")
    plt.tight_layout()
    save(plt.gcf(), "wave1_auroc_primary")

    if random_reps:
        random_summary = summary[summary.Representation.isin(random_reps)]
        uc = summary.set_index("Representation").loc["missing_as_zero_without_GJA3"]
        figure, axis = plt.subplots(figsize=(7, 4.5))
        axis.errorbar(
            ["UC panel"],
            [uc.AUROC_mean],
            yerr=[uc.AUROC_std],
            fmt="o",
            capsize=4,
            label="215-gene UC panel",
        )
        axis.errorbar(
            ["Random 30-gene subsets"],
            [random_summary.AUROC_mean.mean()],
            yerr=[random_summary.AUROC_mean.std(ddof=1) if len(random_summary) > 1 else 0],
            fmt="o",
            capsize=4,
            label="Mean of random subsets",
        )
        axis.axhline(0.5, color="grey", linestyle="--")
        axis.set_ylabel("Mean repetition-level AUROC")
        axis.set_ylim(0.35, 0.9)
        axis.legend()
        save(figure, "wave1_uc_vs_random_subsets")

    pathway_path = OUTPUT_DIR / "pathway_representation_summary.csv"
    if pathway_path.exists():
        pathway = pd.read_csv(pathway_path)
        figure, axis = plt.subplots(figsize=(7.5, 4.8))
        axis.bar(
            pathway.Representation,
            pathway.AUROC_mean,
            yerr=pathway.AUROC_std,
            capsize=4,
            color="#4c6a92",
        )
        axis.axhline(0.5, color="grey", linestyle="--")
        axis.set_ylabel("Mean repetition-level AUROC")
        axis.set_ylim(0.35, 0.9)
        axis.tick_params(axis="x", rotation=18)
        save(figure, "wave1_pathway_vs_gene")

    qc = pd.read_csv(OUTPUT_DIR / "qc_feature_matrix.csv", index_col=0)
    genepy = pd.read_csv(
        ROOT / "analysis/data/genepy_expanded/genepy_original_missing_as_zero_without_GJA3.csv",
        index_col=0,
    )
    qc["Cohort"] = genepy["target"].map({0: "JPT reference", 1: "UC"})
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for axis, column in zip(axes, ["mean_callability", "total_genepy_burden"]):
        sns.boxplot(data=qc, x="Cohort", y=column, ax=axis, color="white")
        sns.stripplot(data=qc, x="Cohort", y=column, ax=axis, alpha=0.55, size=3)
        axis.set_xlabel("")
    save(figure, "wave1_qc_by_cohort")


if __name__ == "__main__":
    main()
