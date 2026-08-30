"""Cross-guide, nearest-control perturbation scoring for the day-14 MVP.

This implements the key logic of the source analysis (PCA-matched AAVS1_1
neighbors and a 100-gene perturbation signature) while cross-fitting signatures
across sgRNAs. A held-out guide is never used to define its own signature.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import yaml


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0, 1)
    return output


def log_normalize(counts: sparse.csr_matrix, scale: float = 10_000.0) -> sparse.csr_matrix:
    output = counts.astype(np.float32, copy=True)
    library_sizes = np.asarray(output.sum(axis=1)).ravel()
    factors = np.divide(
        scale,
        library_sizes,
        out=np.zeros_like(library_sizes, dtype=np.float32),
        where=library_sizes > 0,
    )
    output = output.multiply(factors[:, None]).tocsr()
    output.data = np.log1p(output.data)
    return output


def select_variable_genes(
    expression: sparse.csr_matrix, control_rows: np.ndarray, n_genes: int
) -> np.ndarray:
    controls = expression[control_rows]
    mean = np.asarray(controls.mean(axis=0)).ravel()
    mean_square = np.asarray(controls.power(2).mean(axis=0)).ravel()
    variance = np.maximum(mean_square - mean**2, 0)
    dispersion = variance / np.maximum(mean, 1e-6)
    return np.argsort(dispersion)[-min(n_genes, expression.shape[1]) :]


def reference_neighbors(
    pcs: np.ndarray,
    reference_rows: np.ndarray,
    n_neighbors: int,
) -> tuple[np.ndarray, np.ndarray]:
    reference_pcs = pcs[reference_rows]
    model = NearestNeighbors(n_neighbors=n_neighbors + 1, metric="euclidean")
    model.fit(reference_pcs)
    raw = model.kneighbors(pcs, return_distance=False)
    global_to_local = {int(row): i for i, row in enumerate(reference_rows)}
    neighbors = np.empty((len(pcs), n_neighbors), dtype=np.int32)
    for row, candidates in enumerate(raw):
        self_local = global_to_local.get(row)
        kept = [int(x) for x in candidates if int(x) != self_local][:n_neighbors]
        if len(kept) < n_neighbors:
            raise ValueError("Insufficient non-self reference neighbors")
        neighbors[row] = kept
    return neighbors, reference_rows


def mean_residual(
    expression: sparse.csr_matrix,
    reference_expression: sparse.csr_matrix,
    neighbors: np.ndarray,
    rows: np.ndarray,
) -> np.ndarray:
    observed = np.asarray(expression[rows].mean(axis=0)).ravel()
    weights = np.bincount(
        neighbors[rows].ravel(), minlength=reference_expression.shape[0]
    ).astype(np.float64)
    expected = np.asarray(weights @ reference_expression).ravel() / weights.sum()
    return observed - expected


def residual_block(
    expression: sparse.csr_matrix,
    reference_expression: sparse.csr_matrix,
    neighbors: np.ndarray,
    rows: np.ndarray,
    features: np.ndarray,
) -> np.ndarray:
    observed = expression[rows][:, features].toarray()
    reference = reference_expression[:, features].toarray()
    expected = reference[neighbors[rows]].mean(axis=1)
    return observed - expected


def score(config_path: Path, input_dir: Path, output_dir: Path) -> dict:
    config = load_config(config_path)
    settings = config["perturbation_score"]
    rng = np.random.default_rng(int(config["project"]["random_seed"]))

    counts = sparse.load_npz(input_dir / "rna_counts_cells_by_genes.npz").tocsr()
    metadata = pd.read_csv(input_dir / "cell_metadata.csv.gz", index_col=0)
    features = pd.read_csv(input_dir / "gene_features.csv.gz")
    eligible = features.candidate_eligible.to_numpy(dtype=bool)
    expression = log_normalize(counts[:, eligible])
    eligible_features = features.loc[eligible].reset_index(drop=True)

    guide = metadata.perturbation_name.astype(str).to_numpy()
    target = metadata.target.astype(str).to_numpy()
    discovery_controls = np.flatnonzero(
        np.isin(guide, config["data"]["discovery_control_guides"])
    )
    hvg = select_variable_genes(
        expression, discovery_controls, int(settings["variable_genes"])
    )
    dense_hvg = expression[:, hvg].toarray().astype(np.float32, copy=False)
    dense_hvg = StandardScaler(copy=False).fit_transform(dense_hvg)
    pcs = PCA(
        n_components=int(settings["principal_components"]),
        svd_solver="randomized",
        random_state=int(config["project"]["random_seed"]),
    ).fit_transform(dense_hvg)
    del dense_hvg

    reference_rows = np.flatnonzero(guide == settings["reference_guide"])
    neighbors, reference_rows = reference_neighbors(
        pcs, reference_rows, int(settings["nearest_reference_neighbors"])
    )
    reference_expression = expression[reference_rows]
    null_rows = np.flatnonzero(np.isin(guide, settings["null_guides"]))

    score_targets = list(config["mvp"]["primary_tf_panel"]) + list(
        config["mvp"]["state_dependence_positive_controls"]
    )
    guide_rows = []
    cell_rows = []
    signature_size = int(settings["signature_genes"])
    bootstrap_iterations = int(settings["bootstrap_iterations"])

    for tf in score_targets:
        tf_guides = sorted(np.unique(guide[target == tf]))
        vectors = {
            sg: mean_residual(
                expression,
                reference_expression,
                neighbors,
                np.flatnonzero(guide == sg),
            )
            for sg in tf_guides
        }
        for held_out in tf_guides:
            training_guides = [sg for sg in tf_guides if sg != held_out]
            training_sizes = np.array([np.sum(guide == sg) for sg in training_guides])
            training_vector = np.average(
                np.vstack([vectors[sg] for sg in training_guides]),
                axis=0,
                weights=training_sizes,
            )
            signature = np.argsort(np.abs(training_vector))[-signature_size:]
            held_rows = np.flatnonzero(guide == held_out)
            null_residual = residual_block(
                expression,
                reference_expression,
                neighbors,
                null_rows,
                signature,
            )
            difference_vector = (
                training_vector[signature] - null_residual.mean(axis=0)
            )
            vector_norm = float(np.linalg.norm(difference_vector))
            if vector_norm < 1e-8:
                raise ValueError(f"Degenerate perturbation vector for {tf}/{held_out}")
            null_projection = null_residual @ difference_vector / vector_norm
            projection_center = float(null_projection.mean())
            projection_scale = float(null_projection.std(ddof=1))
            if projection_scale < 1e-8:
                raise ValueError(f"Degenerate null projection for {tf}/{held_out}")
            null_scores = (
                null_projection - projection_center
            ) / projection_scale
            held_residual = residual_block(
                expression,
                reference_expression,
                neighbors,
                held_rows,
                signature,
            )
            held_projection = held_residual @ difference_vector / vector_norm
            held_scores = (
                held_projection - projection_center
            ) / projection_scale
            observed_median = float(np.median(held_scores))
            boot_medians = np.array(
                [
                    np.median(rng.choice(null_scores, size=len(held_scores), replace=True))
                    for _ in range(bootstrap_iterations)
                ]
            )
            null_threshold = float(
                np.quantile(boot_medians, float(settings["null_quantile"]))
            )
            p_value = float(
                (1 + np.sum(boot_medians >= observed_median))
                / (bootstrap_iterations + 1)
            )
            correlation = float(
                spearmanr(vectors[held_out], training_vector).statistic
            )
            guide_rows.append(
                {
                    "target": tf,
                    "guide": held_out,
                    "n_cells": int(len(held_rows)),
                    "crossfit_median_score": observed_median,
                    "bootstrap_null_q95": null_threshold,
                    "bootstrap_p_value": p_value,
                    "direction_spearman": correlation,
                    "signature_genes": ";".join(
                        eligible_features.iloc[signature].feature_key.astype(str)
                    ),
                }
            )
            cell_rows.extend(
                {
                    "cell_id": metadata.index[row],
                    "target": tf,
                    "guide": held_out,
                    "crossfit_perturbation_score": float(value),
                }
                for row, value in zip(held_rows, held_scores)
            )

    guide_scores = pd.DataFrame(guide_rows)
    guide_scores["bootstrap_fdr"] = bh_adjust(
        guide_scores.bootstrap_p_value.to_numpy()
    )
    guide_scores["effective_guide"] = (
        (guide_scores.crossfit_median_score > guide_scores.bootstrap_null_q95)
        & (guide_scores.bootstrap_fdr <= float(settings["guide_fdr_max"]))
        & (
            guide_scores.direction_spearman
            > float(settings["minimum_direction_correlation"])
        )
    )
    target_summary = (
        guide_scores.groupby("target", observed=True)
        .agg(
            guides=("guide", "size"),
            effective_guides=("effective_guide", "sum"),
            median_crossfit_score=("crossfit_median_score", "median"),
            minimum_direction_spearman=("direction_spearman", "min"),
        )
        .reset_index()
    )
    target_summary["passes_two_guide_rule"] = (
        target_summary.effective_guides
        >= int(config["causal_validation"]["minimum_effective_guides"])
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    guide_scores.to_csv(output_dir / "guide_crossfit_scores.csv", index=False)
    target_summary.to_csv(output_dir / "target_efficacy_summary.csv", index=False)
    pd.DataFrame(cell_rows).to_csv(
        output_dir / "cell_crossfit_scores.csv.gz", index=False, compression="gzip"
    )
    np.save(output_dir / "pca_coordinates.npy", pcs)
    pd.DataFrame(
        {
            "cell_id": metadata.index,
            **{f"PC{i + 1}": pcs[:, i] for i in range(pcs.shape[1])},
        }
    ).to_csv(output_dir / "pca_coordinates.csv.gz", index=False, compression="gzip")
    summary = {
        "cells": int(len(metadata)),
        "candidate_genes": int(expression.shape[1]),
        "variable_genes": int(len(hvg)),
        "reference_guide": settings["reference_guide"],
        "reference_cells": int(len(reference_rows)),
        "targets_scored": score_targets,
        "guides_scored": int(len(guide_scores)),
        "effective_guides": int(guide_scores.effective_guide.sum()),
        "targets_passing_two_guide_rule": target_summary.loc[
            target_summary.passes_two_guide_rule, "target"
        ].tolist(),
        "method_note": (
            "Nearest-AAVS1_1 residual signature with guide-held-out 100-gene "
            "cross-fitting and Mixscale-style projection/Z-scoring; efficacy "
            "screen, not an E1 causal-edge test."
        ),
    }
    with (output_dir / "score_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--input", type=Path, default=Path("data/processed/day14_all_guides")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/perturbation_score")
    )
    args = parser.parse_args()
    print(json.dumps(score(args.config, args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
