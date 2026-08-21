"""TPE-tuned comparison of supervised models using repeated nested CV."""

from functools import partial
import json
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier


DATA_PATH = Path(
    "analysis/data/genepy_expanded/genepy_original_missing_as_zero_without_GJA3.csv"
)
OUTPUT_DIR = Path(
    "analysis/results/tpe_model_evaluation_continuous_mi_without_GJA3"
)
N_OUTER_FOLDS = 5
N_OUTER_REPEATS = 10
N_INNER_FOLDS = 3
N_TRIALS = 20
N_ESTIMATORS = 300
RANDOM_STATE = 42
MODEL_NAMES = [
    "Random Forest",
    "Logistic Regression",
    "SVM",
    "KNN",
    "Naive Bayes",
    "XGBoost",
]


def continuous_mutual_information(X, y, random_state):
    """Estimate mutual information from continuous GenePy scores."""
    return mutual_info_classif(
        X,
        y,
        discrete_features=False,
        random_state=random_state,
    )


def feature_count_choices(n_features):
    choices = [value for value in [5, 10, 15, 20, 30, 45] if value < n_features]
    return choices + ["all"]


def suggest_params(trial, model_name, n_features):
    params = {
        "k": trial.suggest_categorical("k", feature_count_choices(n_features))
    }

    if model_name == "Random Forest":
        params.update({
            "criterion": trial.suggest_categorical(
                "criterion", ["gini", "entropy", "log_loss"]
            ),
            "max_depth": trial.suggest_categorical(
                "max_depth", [2, 3, 5, 8, None]
            ),
            "min_samples_split": trial.suggest_int("min_samples_split", 4, 16),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 2, 8),
            "max_features": trial.suggest_categorical(
                "max_features", ["sqrt", "log2", 0.25, 0.5, 1.0]
            ),
            "max_samples": trial.suggest_categorical(
                "max_samples", [0.6, 0.8, 1.0]
            ),
            "class_weight": trial.suggest_categorical(
                "class_weight", [None, "balanced", "balanced_subsample"]
            ),
        })
    elif model_name == "Logistic Regression":
        params.update({
            "C": trial.suggest_float("C", 1e-4, 1e4, log=True),
            "penalty": trial.suggest_categorical("penalty", ["l1", "l2"]),
            "class_weight": trial.suggest_categorical(
                "class_weight", [None, "balanced"]
            ),
        })
    elif model_name == "SVM":
        params.update({
            "C": trial.suggest_float("C", 1e-3, 1e3, log=True),
            "kernel": trial.suggest_categorical("kernel", ["linear", "rbf"]),
            "gamma": trial.suggest_float("gamma", 1e-4, 1.0, log=True),
            "class_weight": trial.suggest_categorical(
                "class_weight", [None, "balanced"]
            ),
        })
    elif model_name == "KNN":
        params.update({
            "n_neighbors": trial.suggest_categorical(
                "n_neighbors", [3, 5, 7, 9, 11, 15]
            ),
            "weights": trial.suggest_categorical(
                "weights", ["uniform", "distance"]
            ),
            "p": trial.suggest_categorical("p", [1, 2]),
        })
    elif model_name == "Naive Bayes":
        params["var_smoothing"] = trial.suggest_float(
            "var_smoothing", 1e-12, 1e-6, log=True
        )
    elif model_name == "XGBoost":
        params.update({
            "n_estimators": trial.suggest_int("n_estimators", 100, 500, step=50),
            "max_depth": trial.suggest_int("max_depth", 2, 6),
            "learning_rate": trial.suggest_float(
                "learning_rate", 0.01, 0.3, log=True
            ),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float(
                "colsample_bytree", 0.5, 1.0
            ),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
        })
    else:
        raise ValueError(f"Unknown model: {model_name}")

    return params


def make_classifier(
    model_name, params, seed, class_ratio, svm_probability=True
):
    if model_name == "Random Forest":
        return RandomForestClassifier(
            n_estimators=N_ESTIMATORS,
            criterion=params["criterion"],
            max_depth=params["max_depth"],
            min_samples_split=params["min_samples_split"],
            min_samples_leaf=params["min_samples_leaf"],
            max_features=params["max_features"],
            max_samples=params["max_samples"],
            class_weight=params["class_weight"],
            random_state=seed,
            n_jobs=-1,
        )
    if model_name == "Logistic Regression":
        return LogisticRegression(
            C=params["C"],
            penalty=params["penalty"],
            solver="liblinear",
            class_weight=params["class_weight"],
            max_iter=2000,
            random_state=seed,
        )
    if model_name == "SVM":
        return SVC(
            C=params["C"],
            kernel=params["kernel"],
            gamma=params["gamma"],
            class_weight=params["class_weight"],
            probability=svm_probability,
            max_iter=100000,
            random_state=seed,
        )
    if model_name == "KNN":
        return KNeighborsClassifier(
            n_neighbors=params["n_neighbors"],
            weights=params["weights"],
            p=params["p"],
        )
    if model_name == "Naive Bayes":
        return GaussianNB(var_smoothing=params["var_smoothing"])
    if model_name == "XGBoost":
        return XGBClassifier(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            learning_rate=params["learning_rate"],
            subsample=params["subsample"],
            colsample_bytree=params["colsample_bytree"],
            min_child_weight=params["min_child_weight"],
            reg_lambda=params["reg_lambda"],
            scale_pos_weight=class_ratio,
            eval_metric="logloss",
            random_state=seed,
            n_jobs=-1,
        )
    raise ValueError(f"Unknown model: {model_name}")


def make_pipeline(
    model_name, params, seed, class_ratio, svm_probability=True
):
    mi_score = partial(continuous_mutual_information, random_state=seed)
    return Pipeline([
        ("imputer", SimpleImputer(
            strategy="median",
            add_indicator=False,
            keep_empty_features=True,
        )),
        ("selector", SelectKBest(score_func=mi_score, k=params["k"])),
        ("scaler", StandardScaler()),
        ("classifier", make_classifier(
            model_name,
            params,
            seed,
            class_ratio,
            svm_probability=svm_probability,
        )),
    ])


def tune(X, y, model_name, seed):
    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS,
        shuffle=True,
        random_state=seed,
    )

    def objective(trial):
        params = suggest_params(trial, model_name, X.shape[1])
        scores = []
        for train_idx, valid_idx in inner_cv.split(X, y):
            inner_y_train = y.iloc[train_idx]
            inner_class_ratio = float(
                (inner_y_train == 0).sum() / (inner_y_train == 1).sum()
            )
            pipeline = make_pipeline(
                model_name,
                params,
                seed,
                inner_class_ratio,
                svm_probability=False,
            )
            pipeline.fit(X.iloc[train_idx], y.iloc[train_idx])
            if model_name == "SVM":
                validation_score = pipeline.decision_function(
                    X.iloc[valid_idx]
                )
            else:
                validation_score = pipeline.predict_proba(
                    X.iloc[valid_idx]
                )[:, 1]
            scores.append(roc_auc_score(y.iloc[valid_idx], validation_score))
        return float(np.mean(scores))

    study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=seed, n_startup_trials=5),
    )
    study.optimize(objective, n_trials=N_TRIALS, show_progress_bar=False)
    return study.best_params, study.best_value


def calculate_metrics(y_true, prediction, probability):
    return {
        "Outer_AUROC": roc_auc_score(y_true, probability),
        "Outer_average_precision": average_precision_score(y_true, probability),
        "Outer_balanced_accuracy": balanced_accuracy_score(y_true, prediction),
        "Outer_accuracy": accuracy_score(y_true, prediction),
        "Outer_precision": precision_score(y_true, prediction, zero_division=0),
        "Outer_sensitivity": recall_score(y_true, prediction, zero_division=0),
        "Outer_specificity": recall_score(
            y_true, prediction, pos_label=0, zero_division=0
        ),
        "Outer_F1": f1_score(y_true, prediction, zero_division=0),
    }


def paired_wilcoxon(repetition_means, summary):
    ordered = summary.sort_values("AUROC_mean", ascending=False)
    best_model = ordered.iloc[0]["Model"]
    runner_up = ordered.iloc[1]["Model"]
    paired = repetition_means.pivot(
        index="Repetition", columns="Model", values="Outer_AUROC"
    )[[best_model, runner_up]].dropna()
    differences = paired[best_model] - paired[runner_up]
    if np.allclose(differences, 0):
        statistic, p_value = 0.0, 1.0
    else:
        statistic, p_value = wilcoxon(
            differences, alternative="two-sided"
        )

    return pd.DataFrame([{
        "Best_model": best_model,
        "Runner_up": runner_up,
        "Repetitions": len(paired),
        "Mean_paired_AUROC_difference": differences.mean(),
        "Median_paired_AUROC_difference": differences.median(),
        "Best_model_wins": int((differences > 0).sum()),
        "Ties": int((differences == 0).sum()),
        "Wilcoxon_statistic": statistic,
        "Wilcoxon_p_value": p_value,
        "Interpretation": "Exploratory; winner and runner-up selected from the same data",
    }])


def rf_vs_all_wilcoxon(repetition_means):
    """Paired exploratory RF comparisons with Holm correction."""
    pivot = repetition_means.pivot(
        index="Repetition", columns="Model", values="Outer_AUROC"
    )
    rows = []
    for comparator in [
        "Naive Bayes",
        "XGBoost",
        "Logistic Regression",
        "KNN",
        "SVM",
    ]:
        paired = pivot[["Random Forest", comparator]].dropna()
        differences = paired["Random Forest"] - paired[comparator]
        if np.allclose(differences, 0):
            statistic, p_value = 0.0, 1.0
        else:
            statistic, p_value = wilcoxon(
                differences, alternative="two-sided"
            )
        rows.append({
            "Comparison": f"Random Forest vs {comparator}",
            "Repetitions": len(paired),
            "Mean_AUROC_difference": differences.mean(),
            "Median_AUROC_difference": differences.median(),
            "Wilcoxon_statistic": statistic,
            "Wilcoxon_p_value": p_value,
        })

    result = pd.DataFrame(rows)
    # Holm step-down family-wise error correction.
    order = np.argsort(result["Wilcoxon_p_value"].to_numpy())
    sorted_p = result.loc[order, "Wilcoxon_p_value"].to_numpy()
    adjusted_sorted = np.maximum.accumulate(
        (len(sorted_p) - np.arange(len(sorted_p))) * sorted_p
    ).clip(max=1.0)
    adjusted = np.empty(len(result))
    adjusted[order] = adjusted_sorted
    result["Holm_adjusted_p_value"] = adjusted
    result["Interpretation"] = (
        "Exploratory; RF selected as best model from the same comparison"
    )
    return result


def main():
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(DATA_PATH, index_col=0)
    X = data.drop(columns="target")
    y = data["target"].astype(int)
    outer_cv = RepeatedStratifiedKFold(
        n_splits=N_OUTER_FOLDS,
        n_repeats=N_OUTER_REPEATS,
        random_state=RANDOM_STATE,
    )
    outer_splits = list(outer_cv.split(X, y))

    performance_rows = []
    prediction_rows = []
    parameter_rows = []
    selection_rows = []

    for model_name in MODEL_NAMES:
        print(f"\nEvaluating {model_name}...")
        for split, (train_idx, test_idx) in enumerate(outer_splits):
            repetition = split // N_OUTER_FOLDS + 1
            fold = split % N_OUTER_FOLDS + 1
            seed = RANDOM_STATE + split

            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            class_ratio = float(
                (y_train == 0).sum() / (y_train == 1).sum()
            )

            best_params, inner_auc = tune(X_train, y_train, model_name, seed)
            pipeline = make_pipeline(
                model_name, best_params, seed, class_ratio
            )
            pipeline.fit(X_train, y_train)

            probability = pipeline.predict_proba(X_test)[:, 1]
            prediction = pipeline.predict(X_test)
            selected = X.columns[
                pipeline["selector"].get_support()
            ].tolist()
            metrics = calculate_metrics(y_test, prediction, probability)

            performance_rows.append({
                "Model": model_name,
                "Repetition": repetition,
                "Fold": fold,
                "Inner_AUROC": inner_auc,
                "Number_of_genes": len(selected),
                **metrics,
            })
            parameter_rows.append({
                "Model": model_name,
                "Repetition": repetition,
                "Fold": fold,
                "Parameters": json.dumps(best_params, sort_keys=True),
            })
            for gene in selected:
                selection_rows.append({
                    "Model": model_name,
                    "Repetition": repetition,
                    "Fold": fold,
                    "Gene": gene,
                })
            for sample, truth, prob, pred in zip(
                X_test.index, y_test, probability, prediction
            ):
                prediction_rows.append({
                    "Model": model_name,
                    "Repetition": repetition,
                    "Fold": fold,
                    "Sample": sample,
                    "True_label": truth,
                    "Probability_UC": prob,
                    "Predicted_label": pred,
                })

            print(
                f"  repetition {repetition}/{N_OUTER_REPEATS}, "
                f"fold {fold}/{N_OUTER_FOLDS}: "
                f"AUROC={metrics['Outer_AUROC']:.3f}, "
                f"genes={len(selected)}"
            )

    performance = pd.DataFrame(performance_rows)
    predictions = pd.DataFrame(prediction_rows)
    parameters = pd.DataFrame(parameter_rows)
    selections = pd.DataFrame(selection_rows)

    repetition_means = (
        performance.groupby(["Model", "Repetition"], as_index=False)
        .mean(numeric_only=True)
    )
    metric_columns = {
        "Outer_AUROC": "AUROC",
        "Outer_average_precision": "Average_precision",
        "Outer_balanced_accuracy": "Balanced_accuracy",
        "Outer_accuracy": "Accuracy",
        "Outer_precision": "Precision",
        "Outer_sensitivity": "Sensitivity",
        "Outer_specificity": "Specificity",
        "Outer_F1": "F1",
    }
    summary_rows = []
    for model_name, group in repetition_means.groupby("Model"):
        row = {"Model": model_name}
        for column, label in metric_columns.items():
            row[f"{label}_mean"] = group[column].mean()
            row[f"{label}_std"] = group[column].std(ddof=1)
        row["Mean_number_of_genes"] = group["Number_of_genes"].mean()
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values(
        "AUROC_mean", ascending=False
    )
    comparison = paired_wilcoxon(repetition_means, summary)
    rf_comparisons = rf_vs_all_wilcoxon(repetition_means)

    selection_frequency = (
        selections.groupby(["Model", "Gene"]).size()
        .div(N_OUTER_FOLDS * N_OUTER_REPEATS)
        .rename("Outer_fold_selection_frequency")
        .reset_index()
        .sort_values(
            ["Model", "Outer_fold_selection_frequency", "Gene"],
            ascending=[True, False, True],
        )
    )

    performance.to_csv(OUTPUT_DIR / "outer_fold_performance.csv", index=False)
    repetition_means.to_csv(
        OUTPUT_DIR / "repetition_mean_performance.csv", index=False
    )
    predictions.to_csv(OUTPUT_DIR / "out_of_fold_predictions.csv", index=False)
    parameters.to_csv(OUTPUT_DIR / "best_parameters.csv", index=False)
    selections.to_csv(OUTPUT_DIR / "selected_genes_by_fold.csv", index=False)
    selection_frequency.to_csv(
        OUTPUT_DIR / "gene_selection_frequency.csv", index=False
    )
    summary.to_csv(OUTPUT_DIR / "model_summary.csv", index=False)
    comparison.to_csv(
        OUTPUT_DIR / "best_vs_runner_up_wilcoxon.csv", index=False
    )
    rf_comparisons.to_csv(
        OUTPUT_DIR / "rf_vs_all_wilcoxon_holm.csv", index=False
    )

    print("\nTPE nested model comparison")
    print(summary.round(4).to_string(index=False))
    print("\nExploratory paired comparison")
    print(comparison.round(4).to_string(index=False))
    print("\nExploratory RF comparisons with Holm correction")
    print(rf_comparisons.round(4).to_string(index=False))
    print(f"\nResults written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
