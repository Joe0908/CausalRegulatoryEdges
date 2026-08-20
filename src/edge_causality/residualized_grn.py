"""Controls-only residualized signed TF-target association network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as student_t
import yaml

from edge_causality.score_perturbations import bh_adjust, log_normalize


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def design_matrix(metadata: pd.DataFrame, library_sizes: np.ndarray) -> np.ndarray:
    categorical = pd.get_dummies(
        metadata[["replicate", "new_CellType"]].astype(str), drop_first=True
    ).to_numpy(dtype=float)
    log_depth = np.log1p(library_sizes).reshape(-1, 1)
    depth_z = (log_depth - log_depth.mean()) / max(float(log_depth.std()), 1e-8)
    return np.column_stack([np.ones(len(metadata)), depth_z, categorical])


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    q, _ = np.linalg.qr(design, mode="reduced")
    return values - q @ (q.T @ values)


def standardized_residuals(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    output = residualize(values, design)
    output -= output.mean(axis=0, keepdims=True)
    scale = output.std(axis=0, ddof=1, keepdims=True)
    scale[scale < 1e-8] = np.nan
    return output / scale


def stratified_bootstrap_indices(
    metadata: pd.DataFrame, rng: np.random.Generator
) -> np.ndarray:
    groups = metadata.groupby(["replicate", "new_CellType"], observed=True).indices
    sampled = [rng.choice(index, size=len(index), replace=True) for index in groups.values()]
    return np.concatenate(sampled)


def infer(config_path: Path, input_dir: Path, output_dir: Path) -> dict:
    config = load_config(config_path)
    settings = config["edge_discovery"]
    rng = np.random.default_rng(int(config["project"]["random_seed"]))
    counts = sparse.load_npz(input_dir / "rna_counts_cells_by_genes.npz").tocsr()
    metadata = pd.read_csv(input_dir / "cell_metadata.csv.gz", index_col=0)
    features = pd.read_csv(input_dir / "gene_features.csv.gz")
    eligible = features.candidate_eligible.to_numpy(dtype=bool)
    selected_features = features.loc[eligible].reset_index(drop=True)
    values = log_normalize(counts[:, eligible]).toarray().astype(np.float32)
    library_sizes = np.asarray(counts.sum(axis=1)).ravel()
    design = design_matrix(metadata, library_sizes)
    residuals = standardized_residuals(values, design)

    tf_names = list(config["mvp"]["primary_tf_panel"])
    feature_lookup = {
        name: int(index)
        for index, name in enumerate(selected_features.gene_name.astype(str))
    }
    missing = [tf for tf in tf_names if tf not in feature_lookup]
    if missing:
        raise ValueError(f"Primary TFs absent after feature filtering: {missing}")
    tf_indices = np.array([feature_lookup[tf] for tf in tf_names])
    correlations = residuals[:, tf_indices].T @ residuals / (len(metadata) - 1)
    correlations = np.clip(correlations, -0.999999, 0.999999)
    degrees_freedom = len(metadata) - np.linalg.matrix_rank(design) - 2
    t_stat = correlations * np.sqrt(
        degrees_freedom / np.maximum(1 - correlations**2, 1e-12)
    )
    p_values = 2 * student_t.sf(np.abs(t_stat), df=degrees_freedom)

    n_tf, n_targets = correlations.shape
    cutoffs = [int(x) for x in settings["evaluation_cutoffs_percent"]]
    selection_counts = {
        cutoff: np.zeros((n_tf, n_targets), dtype=np.int16) for cutoff in cutoffs
    }
    positive_counts = np.zeros((n_tf, n_targets), dtype=np.int16)
    iterations = int(settings["bootstrap_iterations"])
    for _ in range(iterations):
        sampled = stratified_bootstrap_indices(metadata, rng)
        boot_values = values[sampled]
        boot_design = design_matrix(metadata.iloc[sampled], library_sizes[sampled])
        boot_residuals = standardized_residuals(boot_values, boot_design)
        boot_corr = (
            boot_residuals[:, tf_indices].T @ boot_residuals / (len(sampled) - 1)
        )
        positive_counts += boot_corr > 0
        for tf_index, target_index in enumerate(tf_indices):
            boot_corr[tf_index, target_index] = np.nan
        for cutoff in cutoffs:
            top_k = max(1, int(np.ceil(n_targets * cutoff / 100)))
            for tf_index in range(n_tf):
                score = np.abs(boot_corr[tf_index])
                selected = np.argpartition(np.nan_to_num(score, nan=-np.inf), -top_k)[
                    -top_k:
                ]
                selection_counts[cutoff][tf_index, selected] += 1

    rows = []
    target_keys = selected_features.feature_key.astype(str).to_numpy()
    target_symbols = selected_features.gene_name.astype(str).to_numpy()
    target_ids = selected_features.gene_id.astype(str).to_numpy()
    for tf_index, tf in enumerate(tf_names):
        for target_index in range(n_targets):
            if target_index == tf_indices[tf_index]:
                continue
            positive_frequency = positive_counts[tf_index, target_index] / iterations
            row = {
                "TF": tf,
                "target": target_keys[target_index],
                "target_symbol": target_symbols[target_index],
                "target_gene_id": target_ids[target_index],
                "signed_association": float(correlations[tf_index, target_index]),
                "absolute_association": float(abs(correlations[tf_index, target_index])),
                "p_value": float(p_values[tf_index, target_index]),
                "positive_sign_frequency": float(positive_frequency),
                "bootstrap_sign_consistency": float(
                    max(positive_frequency, 1 - positive_frequency)
                ),
            }
            for cutoff in cutoffs:
                row[f"bootstrap_top{cutoff}_frequency"] = float(
                    selection_counts[cutoff][tf_index, target_index] / iterations
                )
            rows.append(row)
    edges = pd.DataFrame(rows)
    edges["fdr_global"] = bh_adjust(edges.p_value.to_numpy())
    selection_cutoff = int(settings["bootstrap_selection_percent"])
    frequency_column = f"bootstrap_top{selection_cutoff}_frequency"
    edges["stable_edge"] = (
        edges[frequency_column] >= float(settings["minimum_bootstrap_frequency"])
    )
    edges = edges.sort_values(
        ["stable_edge", "absolute_association"], ascending=[False, False]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(output_dir / "residualized_signed_edges.csv.gz", index=False, compression="gzip")
    stable = edges.loc[edges.stable_edge].copy()
    stable.to_csv(output_dir / "stable_residualized_edges.csv", index=False)
    summary = {
        "controls_only": True,
        "cells": int(len(metadata)),
        "candidate_genes": int(n_targets),
        "TFs": tf_names,
        "candidate_edges": int(len(edges)),
        "bootstrap_iterations": iterations,
        "stability_definition": (
            f"top {selection_cutoff}% within TF in >= "
            f"{float(settings['minimum_bootstrap_frequency']):.2f} bootstraps"
        ),
        "stable_edges": int(len(stable)),
        "stable_edges_by_TF": {
            str(k): int(v) for k, v in stable.TF.value_counts().items()
        },
        "fdr_0_05_edges": int((edges.fdr_global < 0.05).sum()),
    }
    with (output_dir / "residualized_grn_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--input", type=Path, default=Path("data/processed/controls_day14")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/residualized_grn")
    )
    args = parser.parse_args()
    print(json.dumps(infer(args.config, args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
