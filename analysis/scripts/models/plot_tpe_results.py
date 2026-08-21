"""Create publication figures from the finalized TPE model comparison."""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


RESULTS_DIR = Path(
    "analysis/results/tpe_model_evaluation_missing_as_zero_without_GJA3"
)
FIGURE_DIR = RESULTS_DIR / "figures"
EXPLAIN_DIR = Path(
    "analysis/results/final_model_explainability_expanded_without_GJA3"
)
MATRIX_PATH = Path(
    "analysis/data/genepy_expanded/"
    "genepy_original_missing_as_zero_without_GJA3.csv"
)
CALLABILITY_PATH = Path(
    "analysis/data/genepy_expanded/genepy_callability.csv"
)


def save(figure, stem, directory=FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ["svg", "png"]:
        figure.savefig(
            directory / f"{stem}.{extension}", dpi=300, bbox_inches="tight"
        )
    plt.close(figure)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(RESULTS_DIR / "model_summary.csv")
    repetitions = pd.read_csv(RESULTS_DIR / "repetition_mean_performance.csv")
    predictions = pd.read_csv(RESULTS_DIR / "out_of_fold_predictions.csv")
    order = summary.sort_values("AUROC_mean", ascending=False)["Model"].tolist()
    palette = dict(zip(order, sns.color_palette("colorblind", len(order))))

    plt.figure(figsize=(9, 5.5))
    sns.boxplot(
        data=repetitions,
        x="Model",
        y="Outer_AUROC",
        order=order,
        color="white",
        fliersize=0,
    )
    sns.stripplot(
        data=repetitions,
        x="Model",
        y="Outer_AUROC",
        order=order,
        hue="Model",
        palette=palette,
        size=5,
        jitter=0.16,
        legend=False,
    )
    plt.axhline(0.5, color="grey", linestyle="--", linewidth=1)
    plt.xlabel("")
    plt.ylabel("Repetition-level AUROC")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    save(plt.gcf(), "model_auroc_stability")

    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    for model in order:
        model_predictions = predictions[predictions["Model"] == model]
        pooled_auroc = roc_auc_score(
            model_predictions["True_label"],
            model_predictions["Probability_UC"],
        )
        pooled_ap = average_precision_score(
            model_predictions["True_label"],
            model_predictions["Probability_UC"],
        )
        fpr, tpr, _ = roc_curve(
            model_predictions["True_label"],
            model_predictions["Probability_UC"],
        )
        precision, recall, _ = precision_recall_curve(
            model_predictions["True_label"],
            model_predictions["Probability_UC"],
        )
        axes[0].plot(
            fpr,
            tpr,
            color=palette[model],
            label=f"{model} ({pooled_auroc:.3f})",
        )
        axes[1].plot(
            recall,
            precision,
            color=palette[model],
            label=f"{model} ({pooled_ap:.3f})",
        )

    axes[0].plot([0, 1], [0, 1], color="grey", linestyle="--")
    axes[0].set_xlabel("False-positive rate")
    axes[0].set_ylabel("True-positive rate")
    axes[0].set_title("Receiver operating characteristic")
    prevalence = predictions["True_label"].mean()
    axes[1].axhline(prevalence, color="grey", linestyle="--")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("Precision--recall")
    for axis in axes:
        axis.legend(fontsize=8)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1.02)
    plt.tight_layout()
    save(figure, "model_roc_pr_curves")

    # Direct comparison of the three principal performance metrics.
    long_rows = []
    for row in summary.itertuples(index=False):
        for metric in ["AUROC", "Average_precision", "Balanced_accuracy"]:
            long_rows.append({
                "Model": row.Model,
                "Metric": metric.replace("_", " "),
                "Mean": getattr(row, f"{metric}_mean"),
                "SD": getattr(row, f"{metric}_std"),
            })
    comparison = pd.DataFrame(long_rows)
    figure, axis = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(order))
    width = 0.24
    for offset, metric in enumerate(comparison["Metric"].unique()):
        values = comparison[comparison["Metric"] == metric].set_index("Model").loc[order]
        axis.bar(
            x + (offset - 1) * width,
            values["Mean"],
            width,
            yerr=values["SD"],
            capsize=3,
            label=metric,
        )
    axis.axhline(0.5, color="grey", linestyle="--", linewidth=1)
    axis.set_xticks(x, order, rotation=25, ha="right")
    axis.set_ylabel("Mean repetition-level score")
    axis.set_ylim(0.35, 0.85)
    axis.legend()
    save(figure, "model_performance_comparison")

    # Threshold-dependent behaviour from all repeated out-of-fold predictions.
    figure, axes = plt.subplots(2, 3, figsize=(11, 7))
    for axis, model in zip(axes.flat, order):
        model_predictions = predictions[predictions["Model"] == model]
        matrix = confusion_matrix(
            model_predictions["True_label"],
            model_predictions["Predicted_label"],
            normalize="true",
        )
        sns.heatmap(
            matrix,
            annot=True,
            fmt=".2f",
            cmap="Blues",
            vmin=0,
            vmax=1,
            cbar=False,
            xticklabels=["JPT reference", "UC"],
            yticklabels=["JPT reference", "UC"],
            ax=axis,
        )
        axis.set_title(model)
        axis.set_xlabel("Predicted")
        axis.set_ylabel("Observed")
    figure.tight_layout()
    save(figure, "model_confusion_matrices")

    selections = pd.read_csv(RESULTS_DIR / "gene_selection_frequency.csv")
    selected_models = [model for model in ["Random Forest", "Naive Bayes"] if model in selections.Model.unique()]
    figure, axes = plt.subplots(1, len(selected_models), figsize=(11, 6), squeeze=False)
    for axis, model in zip(axes.flat, selected_models):
        top = (
            selections[selections.Model == model]
            .nlargest(15, "Outer_fold_selection_frequency")
            .sort_values("Outer_fold_selection_frequency")
        )
        axis.barh(top.Gene, top.Outer_fold_selection_frequency, color=palette[model])
        axis.set_xlim(0, 1.02)
        axis.set_title(model)
        axis.set_xlabel("Selection frequency across outer folds")
    figure.tight_layout()
    save(figure, "gene_selection_frequency")

    matrix = pd.read_csv(MATRIX_PATH, index_col=0)
    callability = pd.read_csv(CALLABILITY_PATH, index_col=0).drop(columns="GJA3", errors="ignore")
    callability = callability.loc[matrix.index, matrix.drop(columns="target").columns]
    cohort = matrix["target"].map({0: "JPT control", 1: "UC"})
    callability_summary = pd.DataFrame({
        "Cohort": cohort,
        "Mean gene callability": callability.mean(axis=1),
        "Fully callable gene fraction": callability.eq(1).mean(axis=1),
    })
    callability_long = callability_summary.melt(id_vars="Cohort", var_name="Measure", value_name="Value")
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for axis, measure in zip(axes, callability_long.Measure.unique()):
        subset = callability_long[callability_long.Measure == measure]
        sns.boxplot(data=subset, x="Cohort", y="Value", color="white", ax=axis)
        sns.stripplot(data=subset, x="Cohort", y="Value", hue="Cohort", legend=False, alpha=.55, size=3, ax=axis)
        axis.set_title(measure)
        axis.set_xlabel("")
    figure.tight_layout()
    save(figure, "callability_cohort_diagnostic")

    audit = pd.read_csv(EXPLAIN_DIR / "gene_importance_with_callability_audit.csv")
    rf_audit = audit[audit.Model == "Random Forest"].copy()
    figure, axis = plt.subplots(figsize=(7, 5.5))
    axis.scatter(
        rf_audit.Absolute_callability_difference,
        rf_audit.Mean_absolute_SHAP,
        alpha=.55,
    )
    for row in rf_audit.nlargest(10, "Mean_absolute_SHAP").itertuples():
        axis.annotate(row.Gene, (row.Absolute_callability_difference, row.Mean_absolute_SHAP), fontsize=8)
    rho = rf_audit[["Absolute_callability_difference", "Mean_absolute_SHAP"]].corr(method="spearman").iloc[0, 1]
    axis.set_xlabel("Absolute UC–control callability difference")
    axis.set_ylabel("Mean absolute held-out SHAP")
    axis.set_title(f"RF importance versus callability imbalance (Spearman rho={rho:.2f})")
    save(figure, "importance_vs_callability")

    candidates = (
        audit[audit.Model == "Random Forest"]
        .nlargest(5, "Mean_heldout_importance")["Gene"]
        .tolist()
    )
    gene_long = (
        matrix[candidates]
        .assign(Cohort=cohort)
        .melt(id_vars="Cohort", var_name="Gene", value_name="GenePy score")
    )
    figure, axes = plt.subplots(1, len(candidates), figsize=(15, 4), sharey=False)
    for axis, gene in zip(axes, candidates):
        subset = gene_long[gene_long.Gene == gene]
        sns.boxplot(data=subset, x="Cohort", y="GenePy score", color="white", showfliers=False, ax=axis)
        sns.stripplot(data=subset, x="Cohort", y="GenePy score", hue="Cohort", legend=False, alpha=.45, size=2.5, ax=axis)
        axis.set_title(gene)
        axis.set_xlabel("")
    figure.tight_layout()
    save(figure, "candidate_genepy_distributions")

    candidate_callability = (
        callability[candidates]
        .assign(Cohort=cohort)
        .melt(id_vars="Cohort", var_name="Gene", value_name="Callability")
        .groupby(["Gene", "Cohort"], as_index=False)["Callability"]
        .mean()
    )
    figure, axis = plt.subplots(figsize=(8, 5))
    sns.barplot(
        data=candidate_callability,
        x="Gene",
        y="Callability",
        hue="Cohort",
        order=candidates,
        ax=axis,
    )
    axis.set_ylim(0, 1.03)
    axis.set_ylabel("Mean contributing-locus callability")
    axis.set_title("Callability of the leading RF features by source cohort")
    figure.tight_layout()
    save(figure, "candidate_callability_audit")

    print(f"Figures written to {FIGURE_DIR}")


if __name__ == "__main__":
    main()
