"""Wave 2 nested panel comparison — run only after joint-called GenePy matrices exist.

Compares UC panel vs autoimmune/IBD panel vs all callable genes on identical
outer CV partitions, with size-matched random panels. Does nothing if the
joint-called matrices are absent (Wave 1 / no BAMs).
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from analysis.scripts.paper.wave1_identifiability import (  # noqa: E402
    RANDOM_STATE,
    evaluate_matrix,
    summarise,
)
from sklearn.model_selection import RepeatedStratifiedKFold

JOINT_DIR = ROOT / "analysis/data/genepy_joint"
OUTPUT_DIR = ROOT / "paper/results/wave2"
UC_PANEL = ROOT / "uc-genepy-ml/config/genes_uc.txt"
AUTOIMMUNE_PANEL = ROOT / "paper/config/autoimmune_ibd_panel.txt"


def read_symbols(path):
    text = path.read_text().replace("\n", ",")
    return [symbol.strip() for symbol in text.split(",") if symbol.strip()]


def subset_matrix(matrix, symbols):
    present = [symbol for symbol in symbols if symbol in matrix.columns]
    missing = sorted(set(symbols) - set(present))
    return matrix[present], missing


def random_panel(matrix, n_genes, seed):
    rng = np.random.default_rng(seed)
    chosen = rng.choice(matrix.columns, size=min(n_genes, matrix.shape[1]), replace=False)
    return matrix[chosen], list(chosen)


def main():
    joint_path = JOINT_DIR / "genepy_all_callable.csv"
    if not joint_path.exists():
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        message = {
            "status": "deferred",
            "reason": (
                "Joint-called GenePy matrix not found. Download BAMs and run "
                "paper/docs/bam_gvcf_checklist.md before this comparison."
            ),
            "expected_input": str(joint_path),
        }
        with (OUTPUT_DIR / "nested_panels_status.json").open("w") as handle:
            json.dump(message, handle, indent=2)
        print(message["reason"])
        return

    data = pd.read_csv(joint_path, index_col=0)
    y = data["target"].astype(int)
    X_all = data.drop(columns="target")
    uc_symbols = read_symbols(UC_PANEL)
    autoimmune_symbols = read_symbols(AUTOIMMUNE_PANEL)
    X_uc, missing_uc = subset_matrix(X_all, uc_symbols)
    X_auto, missing_auto = subset_matrix(X_all, autoimmune_symbols)
    X_random_uc, genes_random_uc = random_panel(X_all, X_uc.shape[1], seed=0)
    X_random_auto, genes_random_auto = random_panel(X_all, X_auto.shape[1], seed=1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.Series(genes_random_uc).to_csv(OUTPUT_DIR / "random_uc_width_genes.csv", index=False)
    pd.Series(genes_random_auto).to_csv(
        OUTPUT_DIR / "random_autoimmune_width_genes.csv", index=False
    )
    with (OUTPUT_DIR / "panel_overlap.json").open("w") as handle:
        json.dump(
            {
                "n_all_callable": int(X_all.shape[1]),
                "n_uc_present": int(X_uc.shape[1]),
                "n_autoimmune_present": int(X_auto.shape[1]),
                "n_uc_missing_from_joint": len(missing_uc),
                "n_autoimmune_missing_from_joint": len(missing_auto),
            },
            handle,
            indent=2,
        )

    outer_splits = list(
        RepeatedStratifiedKFold(
            n_splits=5, n_repeats=10, random_state=RANDOM_STATE
        ).split(X_uc, y)
    )
    frames = []
    predictions = []
    for name, matrix in [
        ("uc_open_targets", X_uc),
        ("autoimmune_ibd", X_auto),
        ("all_callable_genes", X_all),
        ("random_uc_width", X_random_uc),
        ("random_autoimmune_width", X_random_auto),
    ]:
        performance, prediction = evaluate_matrix(name, matrix, y, outer_splits)
        frames.append(performance)
        predictions.append(prediction)
    performance = pd.concat(frames, ignore_index=True)
    prediction = pd.concat(predictions, ignore_index=True)
    repetition_means, summary = summarise(performance)
    performance.to_csv(OUTPUT_DIR / "nested_panel_outer_fold.csv", index=False)
    prediction.to_csv(OUTPUT_DIR / "nested_panel_predictions.csv", index=False)
    repetition_means.to_csv(OUTPUT_DIR / "nested_panel_repetitions.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "nested_panel_summary.csv", index=False)
    print(summary.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
