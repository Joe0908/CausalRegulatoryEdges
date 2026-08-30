"""Local, controls-only GRNBoost2-compatible edge inference.

The regression profile and early-stop rule match arboreto 0.1.6, while avoiding
its distributed runtime so the analysis remains reproducible in PID-isolated
containers.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import GradientBoostingRegressor
import yaml

from edge_causality.score_perturbations import log_normalize


class EarlyStopMonitor:
    def __init__(self, window_length: int = 25):
        self.window_length = window_length

    def __call__(self, current_round, regressor, _unused) -> bool:
        if current_round < self.window_length - 1:
            return False
        low = current_round - self.window_length + 1
        return float(np.mean(regressor.oob_improvement_[low : current_round + 1])) < 0


def fit_target(
    target_index: int,
    values: np.ndarray,
    tf_indices: np.ndarray,
    tf_names: list[str],
    seed: int,
) -> tuple[int, np.ndarray, int]:
    keep = tf_indices != target_index
    predictors = values[:, tf_indices[keep]]
    model = GradientBoostingRegressor(
        learning_rate=0.01,
        n_estimators=5000,
        max_features=0.1,
        subsample=0.9,
        random_state=seed,
    )
    model.fit(
        predictors,
        values[:, target_index],
        monitor=EarlyStopMonitor(window_length=25),
    )
    output = np.zeros(len(tf_names), dtype=np.float32)
    output[keep] = model.feature_importances_ * len(model.estimators_)
    return target_index, output, int(len(model.estimators_))


def infer(
    config_path: Path,
    input_dir: Path,
    output_dir: Path,
    workers: int,
) -> dict:
    config = load_config(config_path)
    counts = sparse.load_npz(input_dir / "rna_counts_cells_by_genes.npz").tocsr()
    features = pd.read_csv(input_dir / "gene_features.csv.gz")
    eligible = features.candidate_eligible.to_numpy(dtype=bool)
    selected = features.loc[eligible].reset_index(drop=True)
    values = log_normalize(counts[:, eligible]).toarray().astype(np.float32)
    tf_names = list(config["mvp"]["primary_tf_panel"])
    lookup = {
        name: int(index) for index, name in enumerate(selected.gene_name.astype(str))
    }
    missing = [tf for tf in tf_names if tf not in lookup]
    if missing:
        raise ValueError(f"Primary TFs absent after filtering: {missing}")
    tf_indices = np.array([lookup[tf] for tf in tf_names], dtype=int)
    seed = int(config["project"]["random_seed"])
    target_indices = list(range(values.shape[1]))

    def task(index: int):
        return fit_target(index, values, tf_indices, tf_names, seed + index)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(task, target_indices, chunksize=16))
    importance = np.zeros((len(tf_names), values.shape[1]), dtype=np.float32)
    trees = np.zeros(values.shape[1], dtype=np.int32)
    for target_index, target_importance, n_trees in results:
        importance[:, target_index] = target_importance
        trees[target_index] = n_trees

    rows = []
    for tf_index, tf in enumerate(tf_names):
        for target_index, target_row in selected.iterrows():
            if target_index == tf_indices[tf_index]:
                continue
            rows.append(
                {
                    "TF": tf,
                    "target": str(target_row.feature_key),
                    "target_symbol": str(target_row.gene_name),
                    "target_gene_id": str(target_row.gene_id),
                    "importance": float(importance[tf_index, target_index]),
                    "trees_fitted_for_target": int(trees[target_index]),
                }
            )
    edges = pd.DataFrame(rows)
    edges["rank_within_TF"] = edges.groupby("TF", observed=True).importance.rank(
        method="first", ascending=False
    )
    group_size = edges.groupby("TF", observed=True).target.transform("size")
    edges["percentile_within_TF"] = edges.rank_within_TF / group_size
    edges = edges.sort_values("importance", ascending=False)

    output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(output_dir / "grnboost2_edges.csv.gz", index=False, compression="gzip")
    for cutoff in config["edge_discovery"]["evaluation_cutoffs_percent"]:
        edges.loc[edges.percentile_within_TF <= float(cutoff) / 100].to_csv(
            output_dir / f"grnboost2_top{int(cutoff)}_percent.csv", index=False
        )
    summary = {
        "controls_only": True,
        "implementation": "arboreto_0.1.6_compatible_local_sklearn",
        "cells": int(values.shape[0]),
        "candidate_genes": int(values.shape[1]),
        "TFs": tf_names,
        "candidate_edges": int(len(edges)),
        "workers": int(workers),
        "mean_trees_per_target": float(trees.mean()),
        "median_trees_per_target": float(np.median(trees)),
        "maximum_trees_per_target": int(trees.max()),
        "bootstrap_complete": False,
    }
    with (output_dir / "grnboost2_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--input", type=Path, default=Path("data/processed/controls_day14")
    )
    parser.add_argument("--output", type=Path, default=Path("reports/grnboost2"))
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    print(json.dumps(infer(args.config, args.input, args.output, args.workers), indent=2))


if __name__ == "__main__":
    main()
