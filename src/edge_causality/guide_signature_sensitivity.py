"""Edge-specific leave-target-out sensitivity for guide efficacy signatures.

The committed guide screen is unchanged.  For an E0 edge whose target occurs
in a TF-matched held-out-guide signature, this audit removes that target,
recomputes the affected projection, and repeats the original 24-guide BH
family.  It tests target leakage; it does not define a new efficacy screen.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from edge_causality.score_perturbations import (
    bh_adjust,
    load_config,
    log_normalize,
    mean_residual,
    reference_neighbors,
    residual_block,
    select_variable_genes,
)


def effective_guide_flags(table: pd.DataFrame, settings: dict) -> pd.Series:
    """Reapply the committed efficacy rule after the 24-guide BH update."""
    return (
        table.crossfit_median_score.gt(table.bootstrap_null_q95)
        & table.bootstrap_fdr.le(float(settings["guide_fdr_max"]))
        & table.direction_spearman.gt(float(settings["minimum_direction_correlation"]))
    )


def _prepare(config: dict, input_dir: Path) -> dict:
    settings = config["perturbation_score"]
    counts = sparse.load_npz(input_dir / "rna_counts_cells_by_genes.npz").tocsr()
    metadata = pd.read_csv(input_dir / "cell_metadata.csv.gz", index_col=0)
    features = pd.read_csv(input_dir / "gene_features.csv.gz")
    eligible = features.candidate_eligible.to_numpy(bool)
    expression = log_normalize(counts[:, eligible])
    eligible_features = features.loc[eligible].reset_index(drop=True)
    guide = metadata.perturbation_name.astype(str).to_numpy()
    target = metadata.target.astype(str).to_numpy()
    controls = np.flatnonzero(np.isin(guide, config["data"]["discovery_control_guides"]))
    hvg = select_variable_genes(expression, controls, int(settings["variable_genes"]))
    dense = expression[:, hvg].toarray().astype(np.float32, copy=False)
    dense = StandardScaler(copy=False).fit_transform(dense)
    pcs = PCA(n_components=int(settings["principal_components"]), svd_solver="randomized",
              random_state=int(config["project"]["random_seed"])).fit_transform(dense)
    reference_rows = np.flatnonzero(guide == settings["reference_guide"])
    neighbors, reference_rows = reference_neighbors(
        pcs, reference_rows, int(settings["nearest_reference_neighbors"])
    )
    reference_expression = expression[reference_rows]
    null_rows = np.flatnonzero(np.isin(guide, settings["null_guides"]))
    return locals()


def _guide_context(prepared: dict, config: dict) -> list[dict]:
    expression = prepared["expression"]
    reference_expression = prepared["reference_expression"]
    neighbors = prepared["neighbors"]
    guide = prepared["guide"]
    target = prepared["target"]
    settings = config["perturbation_score"]
    contexts: list[dict] = []
    score_targets = list(config["mvp"]["primary_tf_panel"]) + list(
        config["mvp"]["state_dependence_positive_controls"]
    )
    for tf in score_targets:
        tf_guides = sorted(np.unique(guide[target == tf]))
        vectors = {
            sg: mean_residual(expression, reference_expression, neighbors, np.flatnonzero(guide == sg))
            for sg in tf_guides
        }
        for held_out in tf_guides:
            training_guides = [sg for sg in tf_guides if sg != held_out]
            sizes = np.array([np.sum(guide == sg) for sg in training_guides])
            training = np.average(np.vstack([vectors[sg] for sg in training_guides]), axis=0, weights=sizes)
            signature = np.argsort(np.abs(training))[-int(settings["signature_genes"]):]
            contexts.append({"target": tf, "guide": held_out, "tf_guides": tf_guides,
                             "vectors": vectors, "training_vector": training,
                             "signature": signature})
    return contexts


def _fit_context(context: dict, prepared: dict, config: dict,
                 rng_state: dict, remove_index: int | None = None) -> dict:
    settings = config["perturbation_score"]
    expression = prepared["expression"]
    reference_expression = prepared["reference_expression"]
    neighbors = prepared["neighbors"]
    guide = prepared["guide"]
    null_rows = prepared["null_rows"]
    signature = context["signature"]
    if remove_index is not None:
        signature = signature[signature != remove_index]
    training = context["training_vector"]
    null_residual = residual_block(expression, reference_expression, neighbors, null_rows, signature)
    difference = training[signature] - null_residual.mean(axis=0)
    vector_norm = float(np.linalg.norm(difference))
    null_projection = null_residual @ difference / vector_norm
    center = float(null_projection.mean())
    scale = float(null_projection.std(ddof=1))
    null_scores = (null_projection - center) / scale
    held_rows = np.flatnonzero(guide == context["guide"])
    held_residual = residual_block(expression, reference_expression, neighbors, held_rows, signature)
    held_scores = (held_residual @ difference / vector_norm - center) / scale
    median = float(np.median(held_scores))
    rng = np.random.default_rng()
    rng.bit_generator.state = deepcopy(rng_state)
    boot = np.array([
        np.median(rng.choice(null_scores, size=len(held_scores), replace=True))
        for _ in range(int(settings["bootstrap_iterations"]))
    ])
    return {
        "crossfit_median_score": median,
        "bootstrap_null_q95": float(np.quantile(boot, float(settings["null_quantile"]))),
        "bootstrap_p_value": float((1 + np.sum(boot >= median)) / (len(boot) + 1)),
        "direction_spearman": float(spearmanr(
            context["vectors"][context["guide"]], training
        ).statistic),
        "signature_size_after_removal": int(len(signature)),
    }


def audit(config_path: Path, input_dir: Path, e0_path: Path,
          guide_scores_path: Path, strict_time_path: Path, output_dir: Path) -> dict:
    config = load_config(config_path)
    settings = config["perturbation_score"]
    prepared = _prepare(config, input_dir)
    contexts = _guide_context(prepared, config)
    context_map = {(x["target"], x["guide"]): x for x in contexts}
    original = pd.read_csv(guide_scores_path)

    # Cache the exact RNG state used before every original held-out-guide bootstrap.
    rng = np.random.default_rng(int(config["project"]["random_seed"]))
    rng_states: dict[tuple[str, str], dict] = {}
    reconstructed = []
    for context in contexts:
        key = (context["target"], context["guide"])
        rng_states[key] = deepcopy(rng.bit_generator.state)
        fit = _fit_context(context, prepared, config, rng_states[key])
        reconstructed.append({"target": key[0], "guide": key[1], **fit})
        # Advance the shared generator exactly as in the committed implementation.
        rng.bit_generator.state = deepcopy(rng_states[key])
        null_length = len(prepared["null_rows"])
        held_length = int(np.sum(prepared["guide"] == key[1]))
        for _ in range(int(settings["bootstrap_iterations"])):
            rng.choice(null_length, size=held_length, replace=True)
    reconstructed = pd.DataFrame(reconstructed)
    check = original.merge(reconstructed, on=["target", "guide"], suffixes=("_original", "_recomputed"))
    numeric = ["crossfit_median_score", "bootstrap_null_q95", "bootstrap_p_value", "direction_spearman"]
    maximum_reconstruction_difference = float(max(
        np.nanmax(np.abs(check[f"{c}_original"] - check[f"{c}_recomputed"])) for c in numeric
    ))
    if maximum_reconstruction_difference > 1e-10:
        raise ValueError(f"Committed guide scores were not exactly reconstructed: {maximum_reconstruction_difference}")

    features = prepared["eligible_features"]
    symbol_to_indices = {
        symbol: group.index.to_numpy(dtype=int)
        for symbol, group in features.groupby("gene_name", observed=True)
    }
    e0 = pd.read_csv(e0_path)
    strict_time = pd.read_csv(strict_time_path) if strict_time_path.exists() else pd.DataFrame()
    records = []
    refits = []
    for edge in e0.itertuples(index=False):
        indices = symbol_to_indices.get(edge.target, np.array([], dtype=int))
        tf_contexts = [x for x in contexts if x["target"] == edge.TF]
        affected = [x for x in tf_contexts if np.intersect1d(x["signature"], indices).size]
        effective_affected = [
            x for x in affected
            if bool(original.loc[(original.target.eq(edge.TF)) & (original.guide.eq(x["guide"])), "effective_guide"].iloc[0])
        ]
        updated = original.copy()
        changed_guides = []
        for context in affected:
            remove = int(np.intersect1d(context["signature"], indices)[0])
            fit = _fit_context(context, prepared, config, rng_states[(edge.TF, context["guide"])], remove)
            mask = updated.target.eq(edge.TF) & updated.guide.eq(context["guide"])
            for column in ["crossfit_median_score", "bootstrap_null_q95", "bootstrap_p_value", "direction_spearman"]:
                updated.loc[mask, column] = fit[column]
            refits.append({"TF": edge.TF, "target": edge.target, "guide": context["guide"],
                           "removed_feature_key": features.iloc[remove].feature_key, **fit})
        updated["bootstrap_fdr"] = bh_adjust(updated.bootstrap_p_value.to_numpy())
        updated["effective_guide"] = effective_guide_flags(updated, settings)
        tf_original = original.loc[original.target.eq(edge.TF)].set_index("guide").effective_guide
        tf_updated = updated.loc[updated.target.eq(edge.TF)].set_index("guide").effective_guide
        changed_guides = sorted(tf_original.index[tf_original.ne(tf_updated)].tolist())
        original_testable = int(tf_original.sum()) >= int(config["causal_validation"]["minimum_effective_guides"])
        updated_testable = int(tf_updated.sum()) >= int(config["causal_validation"]["minimum_effective_guides"])
        strict_supported = bool(getattr(edge, "E1_supported", False))
        time_supported = False
        if not strict_time.empty:
            hit = strict_time.loc[strict_time.TF.eq(edge.TF) & strict_time.target.eq(edge.target)]
            time_supported = bool(hit.strict_any_timepoint_support.iloc[0]) if len(hit) else False
        records.append({
            "TF": edge.TF, "target": edge.target,
            "target_in_any_tf_signature": bool(affected),
            "target_in_effective_guide_signature": bool(effective_affected),
            "affected_guides": ";".join(x["guide"] for x in affected),
            "effective_guides_original": int(tf_original.sum()),
            "effective_guides_leave_target_out": int(tf_updated.sum()),
            "effective_guide_membership_changed": bool(changed_guides),
            "changed_guides": ";".join(changed_guides),
            "two_guide_testable_original": original_testable,
            "two_guide_testable_leave_target_out": updated_testable,
            "two_guide_testability_changed": original_testable != updated_testable,
            "strict_E1_supported": strict_supported,
            "strict_time_supported": time_supported,
        })

    edges = pd.DataFrame(records)
    refits = pd.DataFrame(refits)
    output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(output_dir / "leave_target_out_edge_sensitivity.csv.gz", index=False, compression="gzip")
    refits.to_csv(output_dir / "leave_target_out_guide_refits.csv.gz", index=False, compression="gzip")
    summary = {
        "E0_edges": int(len(edges)),
        "targets_overlapping_any_tf_signature": int(edges.target_in_any_tf_signature.sum()),
        "targets_overlapping_effective_guide_signature": int(edges.target_in_effective_guide_signature.sum()),
        "strict_E1_edges_affected": int((edges.strict_E1_supported & edges.target_in_any_tf_signature).sum()),
        "edge_guide_refits": int(len(refits)),
        "edges_with_effective_guide_membership_change": int(edges.effective_guide_membership_changed.sum()),
        "edges_with_two_guide_testability_change": int(edges.two_guide_testability_changed.sum()),
        "strict_E1_edges_with_membership_change": int((edges.strict_E1_supported & edges.effective_guide_membership_changed).sum()),
        "strict_time_supported_edges_with_membership_change": int((edges.strict_time_supported & edges.effective_guide_membership_changed).sum()),
        "maximum_original_score_reconstruction_difference": maximum_reconstruction_difference,
        "interpretation": "post-freeze target-leakage sensitivity; committed guide efficacy and E1 definitions unchanged",
    }
    with (output_dir / "leave_target_out_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument("--input", type=Path, default=Path("data/processed/day14_all_guides"))
    parser.add_argument("--e0", type=Path, default=Path("reports/validation/E0_to_E1_edge_matrix.csv.gz"))
    parser.add_argument("--guide-scores", type=Path, default=Path("reports/perturbation_score/guide_crossfit_scores.csv"))
    parser.add_argument("--strict-time", type=Path, default=Path("reports/time_resolved_support/strict_time_resolved_total_support.csv.gz"))
    parser.add_argument("--output", type=Path, default=Path("reports/perturbation_score"))
    args = parser.parse_args()
    print(json.dumps(audit(args.config, args.input, args.e0, args.guide_scores, args.strict_time, args.output), indent=2))


if __name__ == "__main__":
    main()
