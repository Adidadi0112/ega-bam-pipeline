"""Held-out explanations for the expanded-matrix leading models."""

import argparse
import json
import os
from pathlib import Path
import sys
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
from analysis.scripts.models import rf_auc_experiment as protocol

warnings.filterwarnings(
    "ignore", category=FutureWarning, module="sklearn.linear_model._logistic"
)
warnings.filterwarnings(
    "ignore", category=UserWarning, module="sklearn.linear_model._logistic"
)


RESULTS_DIR = Path(
    "analysis/results/tpe_model_evaluation_missing_as_zero_without_GJA3"
)
DEFAULT_OUTPUT_DIR = Path(
    "analysis/results/final_model_explainability_expanded_without_GJA3"
)
MODELS = ["Random Forest", "Naive Bayes"]
PERMUTATION_REPEATS = 20
SHAP_BACKGROUND_SAMPLES = 20


def load_parameters(path):
    parameters = pd.read_csv(path)
    return {
        (row.Model, int(row.Repetition), int(row.Fold)): json.loads(row.Parameters)
        for row in parameters.itertuples(index=False)
        if row.Model in MODELS
    }


def aggregate_permutation(raw):
    return (
        raw.groupby(["Model", "Gene"], as_index=False)
        .agg(
            Mean_heldout_importance=("Importance_mean", "mean"),
            SD_across_outer_folds=("Importance_mean", "std"),
            Positive_fold_fraction=("Importance_mean", lambda x: (x > 0).mean()),
            Outer_folds=("Importance_mean", "size"),
        )
        .sort_values(
            ["Model", "Mean_heldout_importance"],
            ascending=[True, False],
        )
    )


def aggregate_shap(raw):
    return (
        raw.groupby(["Model", "Gene"], as_index=False)
        .agg(
            Mean_absolute_SHAP=("SHAP_value", lambda x: np.abs(x).mean()),
            Mean_signed_SHAP=("SHAP_value", "mean"),
            SD_SHAP=("SHAP_value", "std"),
            Observations=("SHAP_value", "size"),
        )
        .sort_values(
            ["Model", "Mean_absolute_SHAP"], ascending=[True, False]
        )
    )


def selected_feature_permutation(pipeline, X_test, y_test, seed):
    """Permute selected raw genes; unselected genes have exact zero importance."""
    rng = np.random.default_rng(seed)
    baseline = roc_auc_score(
        y_test, pipeline.predict_proba(X_test)[:, 1]
    )
    selected = set(
        X_test.columns[pipeline["selector"].get_support()]
    )
    rows = []
    for gene in X_test.columns:
        if gene not in selected:
            rows.append((gene, 0.0, 0.0))
            continue
        decreases = []
        for _ in range(PERMUTATION_REPEATS):
            permuted = X_test.copy()
            permuted[gene] = rng.permutation(permuted[gene].to_numpy())
            score = roc_auc_score(
                y_test, pipeline.predict_proba(permuted)[:, 1]
            )
            decreases.append(baseline - score)
        rows.append((gene, np.mean(decreases), np.std(decreases, ddof=1)))
    return rows


def plot_permutation(summary, output_dir):
    top_genes = (
        summary.groupby("Gene")["Mean_heldout_importance"]
        .max()
        .nlargest(15)
        .index
    )
    plot_data = summary[summary["Gene"].isin(top_genes)].copy()
    order = (
        plot_data.groupby("Gene")["Mean_heldout_importance"]
        .max()
        .sort_values()
        .index
    )

    plt.figure(figsize=(9, 7))
    sns.barplot(
        data=plot_data,
        x="Mean_heldout_importance",
        y="Gene",
        hue="Model",
        order=order,
    )
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Mean decrease in held-out AUROC")
    plt.ylabel("Gene")
    plt.title("Held-out permutation importance")
    plt.tight_layout()
    plt.savefig(output_dir / "heldout_permutation_importance.svg")
    plt.savefig(output_dir / "heldout_permutation_importance.png", dpi=300)
    plt.close()


def plot_shap(shap_frame, model_name, feature_names, output_dir):
    shap_frame = shap_frame[shap_frame["Model"] == model_name]
    values = shap_frame.pivot(
        index="Explanation_ID", columns="Gene", values="SHAP_value"
    ).reindex(columns=feature_names)
    data = shap_frame.pivot(
        index="Explanation_ID", columns="Gene", values="GenePy_value"
    ).reindex(columns=feature_names)
    explanation = shap.Explanation(
        values=values.to_numpy(),
        data=data.to_numpy(),
        feature_names=feature_names,
    )
    shap.plots.beeswarm(explanation, max_display=15, show=False)
    plt.title(f"{model_name} held-out SHAP values")
    plt.tight_layout()
    stem = model_name.lower().replace(" ", "_") + "_shap_beeswarm"
    plt.savefig(output_dir / f"{stem}.svg")
    plt.savefig(output_dir / f"{stem}.png", dpi=300)
    plt.close()


def main(max_splits=None, skip_shap=False, output_dir=DEFAULT_OUTPUT_DIR):
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(protocol.DATA_PATH, index_col=0)
    X = data.drop(columns="target")
    y = data["target"].astype(int)
    parameters = load_parameters(RESULTS_DIR / "best_parameters.csv")

    outer_cv = RepeatedStratifiedKFold(
        n_splits=protocol.N_OUTER_FOLDS,
        n_repeats=protocol.N_OUTER_REPEATS,
        random_state=protocol.RANDOM_STATE,
    )
    splits = list(outer_cv.split(X, y))
    if max_splits is not None:
        splits = splits[:max_splits]

    permutation_rows = []
    shap_rows = []
    explanation_id = 0

    for model_name in MODELS:
        print(f"\nExplaining {model_name}...")
        for split, (train_idx, test_idx) in enumerate(splits):
            repetition = split // protocol.N_OUTER_FOLDS + 1
            fold = split % protocol.N_OUTER_FOLDS + 1
            seed = protocol.RANDOM_STATE + split
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            class_ratio = float(
                (y_train == 0).sum() / (y_train == 1).sum()
            )
            params = parameters[(model_name, repetition, fold)]
            pipeline = protocol.make_pipeline(
                model_name, params, seed, class_ratio
            )
            if model_name == "Random Forest":
                pipeline["classifier"].set_params(n_jobs=1)
            pipeline.fit(X_train, y_train)

            importance = selected_feature_permutation(
                pipeline, X_test, y_test, seed
            )
            for gene, mean, std in importance:
                permutation_rows.append({
                    "Model": model_name,
                    "Repetition": repetition,
                    "Fold": fold,
                    "Gene": gene,
                    "Importance_mean": mean,
                    "Importance_std": std,
                })

            if not skip_shap:
                if model_name == "Random Forest":
                    selected_mask = pipeline["selector"].get_support()
                    transformed_test = pipeline[:-1].transform(X_test)
                    explanation = shap.TreeExplainer(
                        pipeline["classifier"]
                    )(transformed_test)
                    selected_values = np.asarray(explanation.values)
                    if selected_values.ndim == 3:
                        selected_values = selected_values[:, :, 1]
                    shap_values = np.zeros((len(X_test), X.shape[1]))
                    shap_values[:, selected_mask] = selected_values
                else:
                    background = X_train.sample(
                        n=min(SHAP_BACKGROUND_SAMPLES, len(X_train)),
                        random_state=seed,
                    )

                    def predict_uc(values):
                        frame = pd.DataFrame(values, columns=X.columns)
                        return pipeline.predict_proba(frame)[:, 1]

                    masker = shap.maskers.Independent(
                        background, max_samples=len(background)
                    )
                    explainer = shap.Explainer(
                        predict_uc,
                        masker,
                        algorithm="permutation",
                        feature_names=X.columns.tolist(),
                        seed=seed,
                    )
                    explanation = explainer(
                        X_test,
                        max_evals=2 * X.shape[1] + 1,
                        silent=True,
                    )
                    shap_values = explanation.values

                for sample_position, sample in enumerate(X_test.index):
                    for gene_position, gene in enumerate(X.columns):
                        shap_rows.append({
                            "Model": model_name,
                            "Explanation_ID": explanation_id,
                            "Repetition": repetition,
                            "Fold": fold,
                            "Sample": sample,
                            "True_label": int(y_test.iloc[sample_position]),
                            "Gene": gene,
                            "GenePy_value": X_test.iloc[
                                sample_position, gene_position
                            ],
                            "SHAP_value": shap_values[
                                sample_position, gene_position
                            ],
                        })
                    explanation_id += 1

            print(
                f"  repetition {repetition}, fold {fold}: "
                f"permutation complete"
                + ("; SHAP complete" if not skip_shap else "")
            )

    permutation_raw = pd.DataFrame(permutation_rows)
    permutation_summary = aggregate_permutation(permutation_raw)
    permutation_raw.to_csv(
        output_dir / "heldout_permutation_importance_by_fold.csv", index=False
    )
    permutation_summary.to_csv(
        output_dir / "heldout_permutation_importance_summary.csv", index=False
    )
    plot_permutation(permutation_summary, output_dir)

    if shap_rows:
        shap_raw = pd.DataFrame(shap_rows)
        shap_summary = aggregate_shap(shap_raw)
        shap_raw.to_csv(output_dir / "shap_values.csv", index=False)
        shap_summary.to_csv(
            output_dir / "shap_summary.csv", index=False
        )
        for model_name in MODELS:
            plot_shap(shap_raw, model_name, X.columns.tolist(), output_dir)

    print(f"\nExplainability results written to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-splits", type=int)
    parser.add_argument("--skip-shap", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    main(args.max_splits, args.skip_shap, args.output_dir)
