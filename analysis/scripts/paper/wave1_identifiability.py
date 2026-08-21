"""Wave 1 representation comparison on the current GenePy matrices.

Random Forest nested CV with thesis outer partitions (seed 42). Does not
download or process BAM files. Computational settings are those frozen in
paper/protocol/CLAIMS_AND_PROTOCOL.md (10 TPE trials, 150 trees).
"""

from __future__ import annotations

from functools import partial
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

DATA_DIR = ROOT / "analysis/data/genepy_expanded"
VCF_FILE = ROOT / "analysis/02_targeting/uc_vqsr_rare_included.vcf.gz"
OUTPUT_DIR = ROOT / "paper/results/wave1"
N_OUTER_FOLDS = 5
N_OUTER_REPEATS = 10
N_INNER_FOLDS = 3
N_ESTIMATORS = 150
RANDOM_STATE = 42
RANDOM_SUBSET_SIZE = 30
RANDOM_SUBSET_SEEDS = list(range(3))
LOCKED_RF = {
    "max_depth": 5,
    "min_samples_split": 4,
    "min_samples_leaf": 2,
    "max_features": "sqrt",
    "class_weight": "balanced",
}


def continuous_mutual_information(X, y, random_state):
    return mutual_info_classif(
        X, y, discrete_features=False, random_state=random_state
    )


def feature_count_choices(n_features):
    choices = [value for value in [5, 10, 15, 20, 30, 45] if value < n_features]
    return choices + ["all"]


def make_pipeline(k, seed):
    steps = [
        ("imputer", SimpleImputer(
            strategy="median",
            add_indicator=False,
            keep_empty_features=True,
        )),
        ("scaler", StandardScaler()),
        ("classifier", RandomForestClassifier(
            n_estimators=N_ESTIMATORS,
            random_state=seed,
            n_jobs=1,
            **LOCKED_RF,
        )),
    ]
    if k != "all":
        mi_score = partial(continuous_mutual_information, random_state=seed)
        steps.insert(1, ("selector", SelectKBest(score_func=mi_score, k=k)))
    return Pipeline(steps)


def choose_k(X, y, seed):
    # Wave 1 compares whole representations. k-search is a thesis protocol detail.
    return "all", float("nan")


def calculate_metrics(y_true, prediction, probability):
    return {
        "Outer_AUROC": roc_auc_score(y_true, probability),
        "Outer_average_precision": average_precision_score(y_true, probability),
        "Outer_balanced_accuracy": balanced_accuracy_score(y_true, prediction),
        "Outer_precision": precision_score(y_true, prediction, zero_division=0),
        "Outer_sensitivity": recall_score(y_true, prediction, zero_division=0),
        "Outer_specificity": recall_score(
            y_true, prediction, pos_label=0, zero_division=0
        ),
        "Outer_F1": f1_score(y_true, prediction, zero_division=0),
    }


def evaluate_matrix(name, X, y, outer_splits):
    performance_rows = []
    prediction_rows = []
    print(f"\n=== {name} ({X.shape[1]} features) ===", flush=True)
    for split, (train_idx, test_idx) in enumerate(outer_splits):
        repetition = split // N_OUTER_FOLDS + 1
        fold = split % N_OUTER_FOLDS + 1
        seed = RANDOM_STATE + split
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        k, inner_auc = choose_k(X_train, y_train, seed)
        pipeline = make_pipeline(k, seed)
        pipeline.fit(X_train, y_train)
        probability = pipeline.predict_proba(X_test)[:, 1]
        prediction = pipeline.predict(X_test)
        metrics = calculate_metrics(y_test, prediction, probability)
        n_selected = X_train.shape[1]
        if "selector" in pipeline.named_steps:
            n_selected = int(pipeline["selector"].get_support().sum())
        performance_rows.append({
            "Representation": name,
            "Repetition": repetition,
            "Fold": fold,
            "Inner_AUROC": inner_auc,
            "Number_of_features_selected": n_selected,
            "Parameters": json.dumps({"k": k, **LOCKED_RF}, sort_keys=True, default=str),
            **metrics,
        })
        for sample, truth, prob, pred in zip(
            X_test.index, y_test, probability, prediction
        ):
            prediction_rows.append({
                "Representation": name,
                "Repetition": repetition,
                "Fold": fold,
                "Sample": sample,
                "True_label": int(truth),
                "Probability_UC": float(prob),
                "Predicted_label": int(pred),
            })
        print(
            f"  rep {repetition} fold {fold}: AUROC={metrics['Outer_AUROC']:.3f}",
            flush=True,
        )
    return pd.DataFrame(performance_rows), pd.DataFrame(prediction_rows)


def summarise(performance):
    repetition_means = (
        performance.groupby(["Representation", "Repetition"], as_index=False)
        .mean(numeric_only=True)
    )
    rows = []
    for name, group in repetition_means.groupby("Representation"):
        rows.append({
            "Representation": name,
            "AUROC_mean": group["Outer_AUROC"].mean(),
            "AUROC_std": group["Outer_AUROC"].std(ddof=1),
            "Balanced_accuracy_mean": group["Outer_balanced_accuracy"].mean(),
            "Balanced_accuracy_std": group["Outer_balanced_accuracy"].std(ddof=1),
            "Sensitivity_mean": group["Outer_sensitivity"].mean(),
            "Specificity_mean": group["Outer_specificity"].mean(),
            "Average_precision_mean": group["Outer_average_precision"].mean(),
            "F1_mean": group["Outer_F1"].mean(),
        })
    return repetition_means, pd.DataFrame(rows).sort_values(
        "AUROC_mean", ascending=False
    )


def load_labelled_matrix(path, drop_genes=None):
    data = pd.read_csv(path, index_col=0)
    y = data["target"].astype(int)
    X = data.drop(columns="target")
    if drop_genes:
        X = X.drop(columns=[gene for gene in drop_genes if gene in X.columns])
    return X, y


def parse_vcf_sample_qc(vcf_path, samples):
    """Per-sample missingness, mean DP, and Ti/Tv from the targeted VCF."""
    if not vcf_path.exists():
        return None
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    stats = {
        sample: {
            "vcf_missing_rate": 0,
            "vcf_called": 0,
            "vcf_depth_sum": 0,
            "vcf_depth_n": 0,
            "ti": 0,
            "tv": 0,
        }
        for sample in samples
    }
    n_sites = 0
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = line.rstrip().split("\t")[9:]
                sample_index = {
                    sample: i for i, sample in enumerate(header_samples)
                    if sample in stats
                }
                continue
            columns = line.rstrip().split("\t")
            ref, alt = columns[3], columns[4]
            if "," in alt:
                continue
            n_sites += 1
            is_ti = (ref, alt) in transitions
            fmt = columns[8].split(":")
            gt_i = fmt.index("GT") if "GT" in fmt else None
            dp_i = fmt.index("DP") if "DP" in fmt else None
            for sample, index in sample_index.items():
                fields = columns[9 + index].split(":")
                gt = fields[gt_i] if gt_i is not None and gt_i < len(fields) else "./."
                if gt in {"./.", ".|.", "."}:
                    stats[sample]["vcf_missing_rate"] += 1
                    continue
                stats[sample]["vcf_called"] += 1
                if dp_i is not None and dp_i < len(fields) and fields[dp_i] not in {".", ""}:
                    stats[sample]["vcf_depth_sum"] += float(fields[dp_i])
                    stats[sample]["vcf_depth_n"] += 1
                if is_ti:
                    stats[sample]["ti"] += 1
                else:
                    stats[sample]["tv"] += 1
    rows = []
    for sample, values in stats.items():
        called = values["vcf_called"]
        rows.append({
            "vcf_missing_rate": values["vcf_missing_rate"] / n_sites if n_sites else np.nan,
            "vcf_mean_depth": (
                values["vcf_depth_sum"] / values["vcf_depth_n"]
                if values["vcf_depth_n"] else np.nan
            ),
            "vcf_titv": (
                values["ti"] / values["tv"] if values["tv"] else np.nan
            ),
            "vcf_called_sites": called,
        })
    return pd.DataFrame(rows, index=list(stats))


def build_qc_matrix(genepy, callability, y):
    genes = [column for column in genepy.columns if column != "target"]
    aligned_call = callability.reindex(index=genepy.index, columns=genes)
    qc = pd.DataFrame(index=genepy.index)
    qc["mean_callability"] = aligned_call.mean(axis=1)
    qc["fully_callable_fraction"] = aligned_call.eq(1).mean(axis=1)
    qc["total_genepy_burden"] = genepy[genes].sum(axis=1)
    qc["nonzero_gene_count"] = (genepy[genes] > 0).sum(axis=1)
    qc["genepy_std"] = genepy[genes].std(axis=1)
    vcf_qc = parse_vcf_sample_qc(VCF_FILE, genepy.index)
    if vcf_qc is not None:
        qc = qc.join(vcf_qc)
    return qc, y


def univariate_ranks(X, y, callability=None):
    rows = []
    for gene in X.columns:
        uc = X.loc[y == 1, gene]
        jpt = X.loc[y == 0, gene]
        statistic, p_value = mannwhitneyu(
            uc, jpt, alternative="two-sided"
        )
        rows.append({
            "Gene": gene,
            "UC_mean": uc.mean(),
            "JPT_mean": jpt.mean(),
            "Mean_difference_UC_minus_JPT": uc.mean() - jpt.mean(),
            "MannWhitney_U": statistic,
            "MannWhitney_p": p_value,
            "Gene_AUROC": max(
                roc_auc_score(y, X[gene].fillna(X[gene].median())),
                1 - roc_auc_score(y, X[gene].fillna(X[gene].median())),
            ),
            "Signed_AUROC": roc_auc_score(
                y, X[gene].fillna(X[gene].median())
            ),
        })
    ranks = pd.DataFrame(rows).sort_values("MannWhitney_p")
    if callability is not None:
        gene_call_diff = (
            callability.reindex(index=X.index, columns=X.columns)
            .groupby(y.map({1: "UC", 0: "JPT"}))
            .mean()
            .T
        )
        if {"UC", "JPT"}.issubset(gene_call_diff.columns):
            ranks = ranks.merge(
                gene_call_diff.assign(
                    Callability_difference=lambda frame: frame["UC"] - frame["JPT"]
                )[["Callability_difference"]].rename_axis("Gene").reset_index(),
                on="Gene",
                how="left",
            )
            valid = ranks.dropna(subset=["Mean_difference_UC_minus_JPT", "Callability_difference"])
            if len(valid) > 5:
                slope, intercept = np.polyfit(
                    valid["Callability_difference"],
                    valid["Mean_difference_UC_minus_JPT"],
                    1,
                )
                ranks["Callability_residualized_mean_difference"] = (
                    ranks["Mean_difference_UC_minus_JPT"]
                    - intercept
                    - slope * ranks["Callability_difference"].fillna(0)
                )
    return ranks


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    X_a, y = load_labelled_matrix(
        DATA_DIR / "genepy_original_missing_as_zero_without_GJA3.csv"
    )
    X_b, y_b = load_labelled_matrix(
        DATA_DIR / "genepy_original_missing_as_zero.csv"
    )
    X_c, y_c = load_labelled_matrix(
        DATA_DIR / "genepy_original_callability_aware.csv",
        drop_genes=["GJA3"],
    )
    callability = pd.read_csv(DATA_DIR / "genepy_callability.csv", index_col=0)
    assert y.equals(y_b) and y.equals(y_c)

    X_qc, _ = build_qc_matrix(X_a.assign(target=y), callability, y)
    X_qc.to_csv(OUTPUT_DIR / "qc_feature_matrix.csv")

    ranks = univariate_ranks(X_a, y, callability)
    ranks.to_csv(OUTPUT_DIR / "univariate_gene_ranks.csv", index=False)

    outer_cv = RepeatedStratifiedKFold(
        n_splits=N_OUTER_FOLDS,
        n_repeats=N_OUTER_REPEATS,
        random_state=RANDOM_STATE,
    )
    outer_splits = list(outer_cv.split(X_a, y))
    with (OUTPUT_DIR / "outer_split_manifest.json").open("w") as handle:
        json.dump(
            {
                "n_splits": N_OUTER_FOLDS,
                "n_repeats": N_OUTER_REPEATS,
                "random_state": RANDOM_STATE,
                "tuning": "locked RF; k=all (whole-representation comparison)",
                "n_estimators": N_ESTIMATORS,
                "locked_rf": LOCKED_RF,
                "n_samples": int(len(y)),
                "n_uc": int((y == 1).sum()),
                "n_jpt": int((y == 0).sum()),
            },
            handle,
            indent=2,
        )

    experiments = [
        ("missing_as_zero_without_GJA3", X_a),
        ("missing_as_zero_with_GJA3", X_b),
        ("callability_aware_without_GJA3", X_c),
        ("qc_source_baseline", X_qc),
    ]
    rng = np.random.default_rng(RANDOM_STATE)
    genes = list(X_a.columns)
    subset_frames = []
    for seed in RANDOM_SUBSET_SEEDS:
        chosen = rng.choice(genes, size=RANDOM_SUBSET_SIZE, replace=False)
        subset_frames.append((f"random_30_genes_seed{seed}", X_a[chosen]))
        pd.Series(chosen, name="Gene").to_csv(
            OUTPUT_DIR / f"random_30_genes_seed{seed}.csv", index=False
        )
    experiments.extend(subset_frames)

    all_performance = []
    all_predictions = []
    for name, X in experiments:
        performance, predictions = evaluate_matrix(name, X, y, outer_splits)
        all_performance.append(performance)
        all_predictions.append(predictions)

    performance = pd.concat(all_performance, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)
    repetition_means, summary = summarise(performance)
    performance.to_csv(OUTPUT_DIR / "outer_fold_performance.csv", index=False)
    predictions.to_csv(OUTPUT_DIR / "out_of_fold_predictions.csv", index=False)
    repetition_means.to_csv(
        OUTPUT_DIR / "repetition_mean_performance.csv", index=False
    )
    summary.to_csv(OUTPUT_DIR / "representation_summary.csv", index=False)
    print("\nWave 1 representation summary")
    print(summary.round(4).to_string(index=False))
    print(f"\nWrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
