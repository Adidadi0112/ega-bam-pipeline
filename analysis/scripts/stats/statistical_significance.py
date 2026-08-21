import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_PATH = "analysis/data/genepy_matrix_v2.csv"
TARGET_COLUMN = "target"
RESULTS_DIR = "results/stats"

N_SPLITS = 5
N_REPEATS = 10
RANDOM_STATE = 42
PRIMARY_METRIC = "AUC"

os.makedirs(RESULTS_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
data = pd.read_csv(DATA_PATH, index_col=0)
X = data.drop(columns=TARGET_COLUMN)
y = data[TARGET_COLUMN]

classes = np.sort(y.unique())
if len(classes) != 2:
    raise ValueError(
        f"This script expects binary classification, but found classes: {classes}"
    )

# Prefer label 1 as the positive class; otherwise use the second sorted class.
positive_class = 1 if 1 in classes else classes[-1]

# -----------------------------------------------------------------------------
# 2. Fixed hyperparameters
#    These are the modal (most frequently selected) settings from the previous
#    nested-CV run. Fold-specific parameters from the old run must not be
#    attached to newly generated repeated-CV folds.
# -----------------------------------------------------------------------------

def make_models():
    return {
        "Naive Bayes": Pipeline([
            ("scaler", StandardScaler()),
            ("classifier", GaussianNB(var_smoothing=1e-10)),
        ]),
        "KNN": Pipeline([
            ("scaler", StandardScaler()),
            ("classifier", KNeighborsClassifier(n_neighbors=3)),
        ]),
    }


# -----------------------------------------------------------------------------
# 3. Repeated stratified cross-validation
#    Both models are evaluated on exactly the same train/test splits.
# -----------------------------------------------------------------------------
cv = RepeatedStratifiedKFold(
    n_splits=N_SPLITS,
    n_repeats=N_REPEATS,
    random_state=RANDOM_STATE,
)

rows = []
total_splits = cv.get_n_splits(X, y)

for split_number, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
    repetition = (split_number - 1) // N_SPLITS + 1
    fold = (split_number - 1) % N_SPLITS + 1

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    for model_name, model in make_models().items():
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        probability_classes = model.classes_
        positive_index = int(np.where(probability_classes == positive_class)[0][0])
        y_prob = model.predict_proba(X_test)[:, positive_index]

        rows.append({
            "Model": model_name,
            "Repetition": repetition,
            "Fold": fold,
            "AUC": roc_auc_score(y_test, y_prob),
            "Accuracy": accuracy_score(y_test, y_pred),
            "Precision": precision_score(
                y_test, y_pred, pos_label=positive_class, zero_division=0
            ),
            "Recall": recall_score(
                y_test, y_pred, pos_label=positive_class, zero_division=0
            ),
            "F1": f1_score(
                y_test, y_pred, pos_label=positive_class, zero_division=0
            ),
        })

    print(
        f"Completed repetition {repetition}/{N_REPEATS}, "
        f"fold {fold}/{N_SPLITS} ({split_number}/{total_splits} splits)"
    )

fold_results = pd.DataFrame(rows)
fold_results.to_csv(
    os.path.join(RESULTS_DIR, "nb_knn_repeated_cv_fold_scores.csv"),
    index=False,
)

# -----------------------------------------------------------------------------
# 4. Average the five folds within each repetition
#    These repetition-level means are the paired observations used by Wilcoxon.
# -----------------------------------------------------------------------------
metric_columns = ["AUC", "Accuracy", "Precision", "Recall", "F1"]

repetition_results = (
    fold_results
    .groupby(["Model", "Repetition"], as_index=False)[metric_columns]
    .mean()
)

repetition_results.to_csv(
    os.path.join(RESULTS_DIR, "nb_knn_repetition_mean_scores.csv"),
    index=False,
)

summary = (
    repetition_results
    .groupby("Model")[metric_columns]
    .agg(["mean", "std"])
)

print("\n--- PERFORMANCE ACROSS REPETITIONS ---")
print(summary.round(4))
summary.to_csv(
    os.path.join(RESULTS_DIR, "nb_knn_repeated_cv_summary.csv")
)

# -----------------------------------------------------------------------------
# 5. Paired Wilcoxon signed-rank test on repetition-level mean AUC
# -----------------------------------------------------------------------------
auc_wide = (
    repetition_results
    .pivot(index="Repetition", columns="Model", values=PRIMARY_METRIC)
    .sort_index()
)

required_models = {"Naive Bayes", "KNN"}
if not required_models.issubset(auc_wide.columns):
    raise RuntimeError(
        f"Missing model results. Expected {required_models}, got {set(auc_wide.columns)}"
    )

nb_auc = auc_wide["Naive Bayes"].to_numpy()
knn_auc = auc_wide["KNN"].to_numpy()
differences = nb_auc - knn_auc

if np.allclose(differences, 0):
    statistic, p_value = 0.0, 1.0
else:
    result = wilcoxon(
        nb_auc,
        knn_auc,
        alternative="two-sided",
        zero_method="wilcox",
        method="auto",
    )
    statistic = float(result.statistic)
    p_value = float(result.pvalue)

wins_nb = int(np.sum(differences > 0))
wins_knn = int(np.sum(differences < 0))
ties = int(np.sum(np.isclose(differences, 0)))

wilcoxon_result = pd.DataFrame([{
    "Metric": PRIMARY_METRIC,
    "Unit of analysis": "Mean of 5 folds within each repetition",
    "Number of paired repetitions": N_REPEATS,
    "Naive Bayes mean": nb_auc.mean(),
    "Naive Bayes std": nb_auc.std(ddof=1),
    "KNN mean": knn_auc.mean(),
    "KNN std": knn_auc.std(ddof=1),
    "Mean difference (NB - KNN)": differences.mean(),
    "Median difference (NB - KNN)": np.median(differences),
    "NB better repetitions": wins_nb,
    "KNN better repetitions": wins_knn,
    "Tied repetitions": ties,
    "Wilcoxon statistic": statistic,
    "Two-sided p-value": p_value,
    "Significant at 0.05": p_value < 0.05,
}])

print("\n--- PAIRED WILCOXON TEST: NAIVE BAYES VS KNN REPETITION-MEAN AUC ---")
print(wilcoxon_result.round(6).to_string(index=False))

wilcoxon_result.to_csv(
    os.path.join(
        RESULTS_DIR,
        "wilcoxon_nb_vs_knn_repetition_mean_auc.csv",
    ),
    index=False,
)

print(
    "\nInterpretation note: the test uses one mean AUC per repetition, not all "
    "individual folds as independent observations. Repeated-CV estimates still "
    "reuse the same dataset, so report this comparison as exploratory rather "
    "than as a fully independent confirmatory experiment."
)
