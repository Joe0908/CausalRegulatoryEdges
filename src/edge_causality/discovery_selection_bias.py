"""Audit whether terminal-state E0 discovery misses preterminal regulation.

The frozen day-14 E0 and state-matched E1 sets are never redefined.  This
post-freeze audit reconstructs matched observational discovery sets at each
collection time, repeats discovery with a power-matched bootstrap, and refits
strict collection-time *total-response* support over the union of discovered
edges.  The latter is deliberately not labelled E1 because it is not the
original state-matched endpoint.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact
import yaml

from edge_causality.residualized_grn import (
    standardized_residuals,
    stratified_bootstrap_indices,
)
from edge_causality.score_perturbations import log_normalize
from edge_causality.state_dependence import feature_keys
from edge_causality.time_resolved_support import (
    _effective_guides,
    add_strict_support_calls,
    fit_tf_strict_total_support,
)


KEYS = ["TF", "target"]


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def slug(timepoint: str) -> str:
    return str(timepoint).replace(" ", "")


def _allocate_integer(total: int, weights: np.ndarray) -> np.ndarray:
    """Largest-remainder integer allocation with an exact total."""
    weights = np.asarray(weights, dtype=float)
    if total < 0 or not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Invalid allocation inputs")
    raw = total * weights / weights.sum()
    output = np.floor(raw).astype(int)
    remainder = total - int(output.sum())
    if remainder:
        order = np.argsort(-(raw - output), kind="stable")
        output[order[:remainder]] += 1
    return output


def balanced_stratified_bootstrap_indices(
    metadata: pd.DataFrame,
    rng: np.random.Generator,
    total_cells: int,
) -> np.ndarray:
    """Sample an exact common N with equal library contribution.

    Within each library, state quotas follow that timepoint's observed state
    mixture.  Sampling remains with replacement within library-by-state strata.
    This controls cell count and library weighting without pretending that the
    developmental state mixture can be made identical across time.
    """
    libraries = sorted(metadata.replicate.astype(str).unique())
    library_quota = _allocate_integer(total_cells, np.ones(len(libraries)))
    sampled: list[np.ndarray] = []
    for library, quota in zip(libraries, library_quota):
        block = metadata.loc[metadata.replicate.astype(str).eq(library)]
        groups = {
            str(name): np.asarray(index, dtype=int)
            for name, index in block.groupby(
                "new_CellType", observed=True
            ).groups.items()
        }
        names = sorted(groups)
        state_sizes = np.array([len(groups[name]) for name in names], dtype=float)
        state_quota = _allocate_integer(int(quota), state_sizes)
        for name, n_draw in zip(names, state_quota):
            if n_draw:
                sampled.append(rng.choice(groups[name], size=int(n_draw), replace=True))
    output = np.concatenate(sampled)
    if len(output) != total_cells:
        raise AssertionError("Power-matched bootstrap did not preserve common N")
    return output


def discovery_design_matrix(
    metadata: pd.DataFrame,
    library_sizes: np.ndarray,
    adjust_state: bool,
) -> np.ndarray:
    columns = ["replicate"] + (["new_CellType"] if adjust_state else [])
    categorical = pd.get_dummies(
        metadata[columns].astype(str), drop_first=True
    ).to_numpy(dtype=float)
    log_depth = np.log1p(library_sizes).reshape(-1, 1)
    depth_z = (log_depth - log_depth.mean()) / max(float(log_depth.std()), 1e-8)
    return np.column_stack([np.ones(len(metadata)), depth_z, categorical])


def _correlations(
    values: np.ndarray,
    metadata: pd.DataFrame,
    library_sizes: np.ndarray,
    tf_positions: np.ndarray,
    adjust_state: bool = True,
) -> np.ndarray:
    residuals = standardized_residuals(
        values, discovery_design_matrix(metadata, library_sizes, adjust_state)
    )
    output = residuals[:, tf_positions].T @ residuals / max(len(metadata) - 1, 1)
    return np.clip(output, -0.999999, 0.999999)


def _bootstrap_frequencies(
    values: np.ndarray,
    metadata: pd.DataFrame,
    library_sizes: np.ndarray,
    tf_positions: np.ndarray,
    iterations: int,
    top_fraction: float,
    rng: np.random.Generator,
    common_cells: int | None,
    adjust_state: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    n_tf = len(tf_positions)
    n_targets = values.shape[1]
    top_k = max(1, int(np.ceil(n_targets * top_fraction)))
    selected = np.zeros((n_tf, n_targets), dtype=np.int16)
    positive = np.zeros((n_tf, n_targets), dtype=np.int16)
    for _ in range(iterations):
        if common_cells is None:
            rows = stratified_bootstrap_indices(metadata, rng)
        else:
            rows = balanced_stratified_bootstrap_indices(
                metadata, rng, common_cells
            )
        boot_values = values[rows]
        boot_metadata = metadata.iloc[rows]
        boot_sizes = library_sizes[rows]
        corr = _correlations(
            boot_values,
            boot_metadata,
            boot_sizes,
            tf_positions,
            adjust_state=adjust_state,
        )
        positive += corr > 0
        for tf_index, target_index in enumerate(tf_positions):
            corr[tf_index, target_index] = np.nan
        for tf_index in range(n_tf):
            score = np.nan_to_num(np.abs(corr[tf_index]), nan=-np.inf)
            positions = np.argpartition(score, -top_k)[-top_k:]
            selected[tf_index, positions] += 1
    return selected / iterations, positive / iterations


def discover_timepoint(
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    feature_indices: np.ndarray,
    edge_universe: pd.DataFrame,
    fixed_targets: pd.DataFrame,
    tf_names: list[str],
    settings: dict,
    seed: int,
    common_cells: int,
) -> tuple[pd.DataFrame, dict]:
    """Reconstruct one timepoint E0 with full and power-matched bootstraps."""
    values = log_normalize(counts[:, feature_indices]).toarray().astype(np.float32)
    library_sizes = np.asarray(counts.sum(axis=1)).ravel()
    target_position = {
        target: index for index, target in enumerate(fixed_targets.target.astype(str))
    }
    tf_positions = np.array([target_position[tf] for tf in tf_names], dtype=int)
    pooled = _correlations(values, metadata, library_sizes, tf_positions)
    full_frequency, full_positive = _bootstrap_frequencies(
        values,
        metadata,
        library_sizes,
        tf_positions,
        int(settings["bootstrap_iterations"]),
        float(settings["top_fraction"]),
        np.random.default_rng(seed),
        None,
    )
    matched_frequency, matched_positive = _bootstrap_frequencies(
        values,
        metadata,
        library_sizes,
        tf_positions,
        int(settings["bootstrap_iterations"]),
        float(settings["top_fraction"]),
        np.random.default_rng(seed + 100_000),
        common_cells,
    )
    unadjusted_pooled = _correlations(
        values, metadata, library_sizes, tf_positions, adjust_state=False
    )
    unadjusted_frequency, unadjusted_positive = _bootstrap_frequencies(
        values,
        metadata,
        library_sizes,
        tf_positions,
        int(
            settings["state_residualization_sensitivity"][
                "bootstrap_iterations"
            ]
        ),
        float(settings["top_fraction"]),
        np.random.default_rng(seed + 200_000),
        None,
        adjust_state=False,
    )
    tf_position = {tf: index for index, tf in enumerate(tf_names)}
    output = edge_universe[["TF", "target", "target_symbol", "target_gene_id"]].copy()
    matrix_rows = np.array([tf_position[tf] for tf in output.TF], dtype=int)
    matrix_columns = np.array(
        [target_position[target] for target in output.target], dtype=int
    )
    output["signed_association"] = pooled[matrix_rows, matrix_columns]
    output["absolute_association"] = output.signed_association.abs()
    output["rank_fraction"] = output.groupby(
        "TF", observed=True
    ).absolute_association.rank(method="average", ascending=False, pct=True)
    output["bootstrap_top5_frequency"] = full_frequency[
        matrix_rows, matrix_columns
    ]
    output["bootstrap_sign_consistency"] = np.maximum(
        full_positive[matrix_rows, matrix_columns],
        1 - full_positive[matrix_rows, matrix_columns],
    )
    output["E0_timepoint"] = (
        output.bootstrap_top5_frequency
        >= float(settings["minimum_bootstrap_frequency"])
    )
    output["power_matched_bootstrap_top5_frequency"] = matched_frequency[
        matrix_rows, matrix_columns
    ]
    output["power_matched_bootstrap_sign_consistency"] = np.maximum(
        matched_positive[matrix_rows, matrix_columns],
        1 - matched_positive[matrix_rows, matrix_columns],
    )
    output["E0_timepoint_power_matched"] = (
        output.power_matched_bootstrap_top5_frequency
        >= float(settings["minimum_bootstrap_frequency"])
    )
    output["state_unadjusted_signed_association"] = unadjusted_pooled[
        matrix_rows, matrix_columns
    ]
    output["state_unadjusted_absolute_association"] = output[
        "state_unadjusted_signed_association"
    ].abs()
    output["state_unadjusted_rank_fraction"] = output.groupby(
        "TF", observed=True
    ).state_unadjusted_absolute_association.rank(
        method="average", ascending=False, pct=True
    )
    output["state_unadjusted_bootstrap_top5_frequency"] = unadjusted_frequency[
        matrix_rows, matrix_columns
    ]
    output["state_unadjusted_bootstrap_sign_consistency"] = np.maximum(
        unadjusted_positive[matrix_rows, matrix_columns],
        1 - unadjusted_positive[matrix_rows, matrix_columns],
    )
    output["state_unadjusted_E0_timepoint"] = (
        output.state_unadjusted_bootstrap_top5_frequency
        >= float(settings["minimum_bootstrap_frequency"])
    )
    summary = {
        "cells": int(len(metadata)),
        "power_matched_cells_per_bootstrap": int(common_cells),
        "libraries": int(metadata.replicate.nunique()),
        "library_counts": {
            str(key): int(value)
            for key, value in metadata.replicate.value_counts().sort_index().items()
        },
        "state_counts": {
            str(key): int(value)
            for key, value in metadata.new_CellType.value_counts().items()
        },
        "E0_edges": int(output.E0_timepoint.sum()),
        "power_matched_E0_edges": int(output.E0_timepoint_power_matched.sum()),
        "state_unadjusted_E0_edges": int(output.state_unadjusted_E0_timepoint.sum()),
    }
    return output, summary


def classify_discovery_pattern(pattern: str) -> str:
    pattern = str(pattern).zfill(4)
    if pattern == "1111":
        return "persistent_discovery"
    if pattern == "0001":
        return "day14_only_discovery"
    if pattern == "0000":
        return "never_selected"
    if pattern[-1] == "0":
        active = [index for index, value in enumerate(pattern[:3]) if value == "1"]
        if active and max(active) <= 1:
            return "early_only_discovery"
        if active == [1]:
            return "transient_midpoint_discovery"
        if active == [2]:
            return "day11_only_preterminal_discovery"
        return "preterminal_only_discovery"
    return "cross_time_discovery"


def build_discovery_master(
    edge_universe: pd.DataFrame,
    blocks: dict[str, pd.DataFrame],
    frozen_e0: pd.DataFrame,
    timepoints: list[str],
) -> pd.DataFrame:
    output = edge_universe[
        [
            "TF",
            "target",
            "target_symbol",
            "target_gene_id",
            "signed_association",
            "absolute_association",
            "stable_edge",
        ]
    ].rename(
        columns={
            "signed_association": "signed_association_day14_frozen",
            "absolute_association": "absolute_association_day14_frozen",
            "stable_edge": "frozen_day14_E0",
        }
    )
    frozen_keys = pd.MultiIndex.from_frame(frozen_e0[KEYS])
    output["frozen_day14_E0"] = pd.MultiIndex.from_frame(output[KEYS]).isin(
        frozen_keys
    )
    for timepoint in timepoints:
        local_slug = slug(timepoint)
        raw_block = blocks[timepoint]
        block = raw_block[
            KEYS
            + [
                "signed_association",
                "absolute_association",
                "rank_fraction",
                "bootstrap_top5_frequency",
                "bootstrap_sign_consistency",
                "E0_timepoint",
                "power_matched_bootstrap_top5_frequency",
                "power_matched_bootstrap_sign_consistency",
                "E0_timepoint_power_matched",
                "state_unadjusted_signed_association",
                "state_unadjusted_absolute_association",
                "state_unadjusted_rank_fraction",
                "state_unadjusted_bootstrap_top5_frequency",
                "state_unadjusted_bootstrap_sign_consistency",
                "state_unadjusted_E0_timepoint",
            ]
        ].rename(
            columns={
                column: f"{column}_{local_slug}"
                for column in raw_block.columns
                if column not in KEYS
            }
        )
        output = output.merge(block, on=KEYS, validate="one_to_one")
    pattern_columns = [
        f"E0_timepoint_{slug(timepoint)}" for timepoint in timepoints[:-1]
    ]
    output["discovery_pattern"] = output.apply(
        lambda row: "".join("1" if bool(row[column]) else "0" for column in pattern_columns)
        + ("1" if bool(row.frozen_day14_E0) else "0"),
        axis=1,
    )
    output["power_matched_discovery_pattern"] = output.apply(
        lambda row: "".join(
            "1"
            if bool(row[f"E0_timepoint_power_matched_{slug(timepoint)}"])
            else "0"
            for timepoint in timepoints[:-1]
        )
        + ("1" if bool(row.frozen_day14_E0) else "0"),
        axis=1,
    )
    output["discovery_class"] = output.discovery_pattern.map(
        classify_discovery_pattern
    )
    output["any_preterminal_E0"] = output[pattern_columns].any(axis=1)
    output["preterminal_E0_omitted_by_day14"] = (
        output.any_preterminal_E0 & ~output.frozen_day14_E0
    )
    matched_columns = [
        f"E0_timepoint_power_matched_{slug(timepoint)}"
        for timepoint in timepoints[:-1]
    ]
    output["power_matched_preterminal_E0"] = output[matched_columns].any(axis=1)
    output["power_robust_preterminal_omission"] = (
        output.preterminal_E0_omitted_by_day14
        & output.power_matched_preterminal_E0
    )
    unadjusted_columns = [
        f"state_unadjusted_E0_timepoint_{slug(timepoint)}"
        for timepoint in timepoints[:-1]
    ]
    output["state_unadjusted_preterminal_E0"] = output[
        unadjusted_columns
    ].any(axis=1)
    output["state_unadjusted_preterminal_omission"] = (
        output.state_unadjusted_preterminal_E0 & ~output.frozen_day14_E0
    )
    output["state_unadjusted_only_preterminal_omission"] = (
        output.state_unadjusted_preterminal_omission
        & ~output.any_preterminal_E0
    )
    return output


def prepare_union_candidates(
    master: pd.DataFrame,
    validation: pd.DataFrame,
    timepoints: list[str],
    selection_prefix: str = "E0_timepoint_",
    association_prefix: str = "signed_association_",
) -> pd.DataFrame:
    discovery_columns = [
        f"{selection_prefix}{slug(timepoint)}" for timepoint in timepoints[:-1]
    ]
    selected = master.loc[
        master[discovery_columns].any(axis=1) | master.frozen_day14_E0
    ].copy()
    guide_by_tf = (
        validation.groupby("TF", observed=True)[
            ["effective_guides_used", "effective_guide_names"]
        ]
        .first()
        .reset_index()
    )
    e1 = validation[KEYS + ["E1_supported", "E1_direction_concordant"]]
    selected = selected.merge(guide_by_tf, on="TF", validate="many_to_one")
    selected = selected.merge(e1, on=KEYS, how="left", validate="one_to_one")
    selected["E1_supported"] = selected.E1_supported.fillna(False)
    selected["E1_direction_concordant"] = selected.E1_direction_concordant.fillna(False)

    anchor_values = []
    anchor_timepoints = []
    for row in selected.itertuples(index=False):
        candidates: list[tuple[str, float]] = []
        row_dict = row._asdict()
        for timepoint in timepoints[:-1]:
            local_slug = slug(timepoint)
            if bool(row_dict[f"{selection_prefix}{local_slug}"]):
                candidates.append(
                    (
                        timepoint,
                        float(row_dict[f"{association_prefix}{local_slug}"]),
                    )
                )
        if bool(row_dict["frozen_day14_E0"]):
            candidates.append(
                ("day 14", float(row_dict["signed_association_day14_frozen"]))
            )
        anchor, value = max(candidates, key=lambda item: abs(item[1]))
        anchor_timepoints.append(anchor)
        anchor_values.append(value)
    selected["discovery_anchor_timepoint"] = anchor_timepoints
    selected["signed_association"] = anchor_values
    return selected


def fit_union_total_response(
    config: dict,
    pseudobulk_dir: Path,
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    timepoints = list(config["state_dependence"]["ordered_timepoints"])
    reference_guides = list(config["data"]["intervention_reference_guides"])
    z = np.load(pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_counts.npz")
    counts = z["counts"].astype(np.float64)
    keys = feature_keys(z["gene_name"].astype(str), z["gene_id"].astype(str))
    lookup = {key: index for index, key in enumerate(keys)}
    groups = pd.read_csv(
        pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_groups.csv"
    )
    library_total = counts.sum(axis=1, keepdims=True)
    log_cpm = np.log2(counts / np.maximum(library_total, 1) * 1_000_000 + 0.5)
    candidates = candidates.copy()
    candidates["feature_index"] = candidates.target.map(lookup)
    if candidates.feature_index.isna().any():
        raise ValueError("Union targets absent from pseudobulk features")
    blocks = []
    for _, block in candidates.groupby("TF", observed=True):
        blocks.append(
            fit_tf_strict_total_support(
                block,
                log_cpm,
                groups,
                _effective_guides(block),
                reference_guides,
                timepoints,
                int(config["state_dependence"]["minimum_cells_per_pseudobulk"]),
            )
        )
    return add_strict_support_calls(
        pd.concat(blocks, ignore_index=True),
        config["time_resolved_support"],
        timepoints,
    )


def add_timepoint_observational_direction(
    table: pd.DataFrame,
    timepoints: list[str],
    association_prefix: str,
    output_prefix: str,
) -> pd.DataFrame:
    output = table.copy()
    for timepoint in timepoints:
        local_slug = slug(timepoint)
        output[f"{output_prefix}{local_slug}"] = (
            np.sign(output[f"{association_prefix}{local_slug}"])
            == -np.sign(output[f"effect_{timepoint}"])
        )
    return output


def _fisher_rate(table: pd.DataFrame, exposed: pd.Series, outcome: pd.Series) -> dict:
    exposed = exposed.to_numpy(dtype=bool)
    outcome = outcome.to_numpy(dtype=bool)
    contingency = np.array(
        [
            [np.sum(exposed & outcome), np.sum(exposed & ~outcome)],
            [np.sum(~exposed & outcome), np.sum(~exposed & ~outcome)],
        ]
    )
    result = fisher_exact(contingency)
    return {
        "contingency": contingency.tolist(),
        "odds_ratio": float(result.statistic),
        "p_value": float(result.pvalue),
    }


def summarize_audit(
    master: pd.DataFrame,
    intervention: pd.DataFrame,
    sample_summary: dict,
    timepoints: list[str],
    frozen_time_support: pd.DataFrame,
) -> dict:
    preterminal = timepoints[:-1]
    testable = intervention.strict_time_resolved_testable.fillna(False)
    early_support = intervention[
        [f"strict_total_support_{timepoint}" for timepoint in preterminal]
    ].any(axis=1)
    omitted = ~intervention.frozen_day14_E0
    matched_time_support = np.zeros(len(intervention), dtype=bool)
    matched_time_direction = np.zeros(len(intervention), dtype=bool)
    for timepoint in preterminal:
        local_match = (
            intervention[f"E0_timepoint_{slug(timepoint)}"].to_numpy(dtype=bool)
            & intervention[f"strict_total_support_{timepoint}"].to_numpy(dtype=bool)
        )
        matched_time_support |= local_match
        matched_time_direction |= (
            local_match
            & intervention[
                f"observational_direction_concordant_{slug(timepoint)}"
            ].to_numpy(dtype=bool)
        )
    early_supported_n = int(early_support.sum())
    early_omitted_n = int((early_support & omitted).sum())
    matched_supported_n = int(matched_time_support.sum())
    matched_omitted_n = int((matched_time_support & omitted).sum())
    comparison = intervention.loc[testable].copy()
    comparison_outcome = comparison[
        [f"strict_total_support_{timepoint}" for timepoint in preterminal]
    ].any(axis=1)
    fisher = _fisher_rate(
        comparison,
        ~comparison.frozen_day14_E0,
        comparison_outcome,
    )
    return {
        "status": "completed",
        "original_E0_unchanged": True,
        "original_state_matched_E1_unchanged": True,
        "candidate_edges": int(len(master)),
        "sample_summary": sample_summary,
        "timepoint_E0_edges": {
            timepoint: int(master[f"E0_timepoint_{slug(timepoint)}"].sum())
            for timepoint in timepoints
        },
        "power_matched_timepoint_E0_edges": {
            timepoint: int(
                master[f"E0_timepoint_power_matched_{slug(timepoint)}"].sum()
            )
            for timepoint in timepoints
        },
        "day14_reproduction": {
            "frozen_edges": int(master.frozen_day14_E0.sum()),
            "reconstructed_edges": int(master.E0_timepoint_day14.sum()),
            "membership_disagreements": int(
                (master.frozen_day14_E0 != master.E0_timepoint_day14).sum()
            ),
        },
        "discovery_classes": {
            str(key): int(value)
            for key, value in master.discovery_class.value_counts().items()
        },
        "preterminal_E0_omitted_by_day14": int(
            master.preterminal_E0_omitted_by_day14.sum()
        ),
        "power_robust_preterminal_omission": int(
            master.power_robust_preterminal_omission.sum()
        ),
        "union_edges_tested": int(len(intervention)),
        "union_testable_edges": int(testable.sum()),
        "union_edge_time_hypotheses": int(testable.sum() * len(timepoints)),
        "union_support_by_timepoint": {
            timepoint: int(intervention[f"strict_total_support_{timepoint}"].sum())
            for timepoint in timepoints
        },
        "frozen_edges_within_union_support_by_timepoint": {
            timepoint: int(
                (
                    intervention.frozen_day14_E0
                    & intervention[f"strict_total_support_{timepoint}"]
                ).sum()
            )
            for timepoint in timepoints
        },
        "union_any_preterminal_support": early_supported_n,
        "early_supported_omitted_by_day14": early_omitted_n,
        "early_supported_omission_rate": (
            float(early_omitted_n / early_supported_n) if early_supported_n else None
        ),
        "matched_time_discovery_and_support": matched_supported_n,
        "matched_time_supported_omitted_by_day14": matched_omitted_n,
        "matched_time_supported_omission_rate": (
            float(matched_omitted_n / matched_supported_n)
            if matched_supported_n
            else None
        ),
        "matched_time_direction_concordant": int(matched_time_direction.sum()),
        "power_robust_matched_time_supported_omissions": int(
            (
                matched_time_support
                & omitted
                & intervention.power_matched_preterminal_E0
            ).sum()
        ),
        "power_robust_early_supported_omissions": int(
            (
                early_support
                & omitted
                & intervention.power_matched_preterminal_E0
            ).sum()
        ),
        "omitted_vs_day14_selected_early_support_fisher": fisher,
        "frozen_619_support_by_timepoint": {
            timepoint: int(
                frozen_time_support[f"strict_total_support_{timepoint}"].sum()
            )
            for timepoint in timepoints
        },
        "multiplicity": "BH across all finite union edge-by-time total-response tests",
        "interpretation_boundary": (
            "The audit estimates discovery-stage omission within six TFs and the "
            "guide-testable atlas. It does not prove direct TF-target binding or "
            "recover edges for TFs without two effective guides."
        ),
    }


def summarize_state_residualization_sensitivity(
    master: pd.DataFrame,
    intervention: pd.DataFrame,
    timepoints: list[str],
) -> dict:
    preterminal = timepoints[:-1]
    testable = intervention.strict_time_resolved_testable.fillna(False)
    early_support = intervention[
        [f"strict_total_support_{timepoint}" for timepoint in preterminal]
    ].any(axis=1)
    matched = np.zeros(len(intervention), dtype=bool)
    matched_direction = np.zeros(len(intervention), dtype=bool)
    for timepoint in preterminal:
        local_slug = slug(timepoint)
        local = (
            intervention[f"state_unadjusted_E0_timepoint_{local_slug}"].to_numpy(
                dtype=bool
            )
            & intervention[f"strict_total_support_{timepoint}"].to_numpy(dtype=bool)
        )
        matched |= local
        matched_direction |= (
            local
            & intervention[
                f"state_unadjusted_direction_concordant_{local_slug}"
            ].to_numpy(dtype=bool)
        )
    omitted = ~intervention.frozen_day14_E0
    state_only = intervention.state_unadjusted_only_preterminal_omission.to_numpy(
        dtype=bool
    )
    return {
        "timepoint_E0_edges": {
            timepoint: int(
                master[f"state_unadjusted_E0_timepoint_{slug(timepoint)}"].sum()
            )
            for timepoint in timepoints
        },
        "preterminal_E0_omitted_by_day14": int(
            master.state_unadjusted_preterminal_omission.sum()
        ),
        "state_unadjusted_only_preterminal_omission": int(
            master.state_unadjusted_only_preterminal_omission.sum()
        ),
        "union_edges_tested": int(len(intervention)),
        "union_testable_edges": int(testable.sum()),
        "union_support_by_timepoint": {
            timepoint: int(intervention[f"strict_total_support_{timepoint}"].sum())
            for timepoint in timepoints
        },
        "any_preterminal_support": int(early_support.sum()),
        "early_supported_omitted_by_day14": int((early_support & omitted).sum()),
        "matched_time_discovery_and_support": int(matched.sum()),
        "matched_time_supported_omitted_by_day14": int((matched & omitted).sum()),
        "matched_time_direction_concordant": int(matched_direction.sum()),
        "state_unadjusted_only_matched_time_supported": int(
            (matched & state_only).sum()
        ),
        "state_unadjusted_only_matched_time_direction_concordant": int(
            (matched_direction & state_only).sum()
        ),
        "interpretation": (
            "Removing state from the observational adjustment changes the estimand "
            "to include between-state developmental covariance; it is a bounded "
            "sensitivity analysis, not a replacement E0 definition."
        ),
    }


def run(config_path: Path, mvp_path: Path, output_dir: Path) -> dict:
    settings = load_yaml(config_path)
    mvp = load_yaml(mvp_path)
    inputs = settings["inputs"]
    controls = Path(inputs["controls"])
    counts = sparse.load_npz(controls / "rna_counts_cells_by_genes.npz").tocsr()
    metadata = pd.read_csv(controls / "cell_metadata.csv.gz", index_col=0)
    features = pd.read_csv(controls / "gene_features.csv.gz")
    edge_universe = pd.read_csv(inputs["candidate_universe"])
    frozen_e0 = pd.read_csv(inputs["frozen_E0"])
    validation = pd.read_csv(inputs["validation_matrix"])
    timepoints = list(settings["estimand"]["ordered_timepoints"])

    if len(edge_universe) != int(settings["estimand"]["candidate_edges"]):
        raise ValueError("Candidate edge universe differs from frozen design")
    fixed_targets = edge_universe[
        ["target", "target_symbol", "target_gene_id"]
    ].drop_duplicates("target")
    feature_lookup = dict(
        zip(features.feature_key.astype(str), features.index.to_numpy(dtype=int))
    )
    missing = sorted(set(fixed_targets.target.astype(str)) - set(feature_lookup))
    if missing:
        raise ValueError(f"{len(missing)} fixed targets absent from controls")
    feature_indices = np.array(
        [feature_lookup[target] for target in fixed_targets.target.astype(str)],
        dtype=int,
    )
    tf_names = list(mvp["mvp"]["primary_tf_panel"])
    allowed_states = set(settings["estimand"]["erythroid_cell_types"])
    controls_guides = set(settings["estimand"]["control_guides"])
    minimum_library = int(settings["estimand"]["minimum_cells_per_library"])

    selected_rows: dict[str, np.ndarray] = {}
    selected_metadata: dict[str, pd.DataFrame] = {}
    for timepoint in timepoints:
        states = (
            set(settings["estimand"]["original_day14_cell_types"])
            if timepoint == "day 14"
            else allowed_states
        )
        selected = (
            metadata.Timepoint.eq(timepoint)
            & metadata.new_CellType.isin(states)
            & metadata.perturbation_name.isin(controls_guides)
        )
        library_counts = metadata.loc[selected, "replicate"].value_counts()
        retained = set(library_counts.index[library_counts >= minimum_library])
        selected &= metadata.replicate.isin(retained)
        rows = np.flatnonzero(selected.to_numpy())
        selected_rows[timepoint] = rows
        selected_metadata[timepoint] = metadata.iloc[rows].reset_index(drop=True)
    common_cells = min(len(block) for block in selected_metadata.values())

    blocks: dict[str, pd.DataFrame] = {}
    sample_summary: dict[str, dict] = {}
    base_seed = int(settings["analysis"]["random_seed"])
    for index, timepoint in enumerate(timepoints):
        # Preserve the original day-14 RNG stream for an exact reproduction check.
        seed = base_seed if timepoint == "day 14" else base_seed + (index + 1) * 1000
        block, block_summary = discover_timepoint(
            counts[selected_rows[timepoint]],
            selected_metadata[timepoint],
            feature_indices,
            edge_universe,
            fixed_targets,
            tf_names,
            settings["discovery"],
            seed,
            common_cells,
        )
        blocks[timepoint] = block
        sample_summary[timepoint] = block_summary

    master = build_discovery_master(
        edge_universe, blocks, frozen_e0, timepoints
    )
    candidates = prepare_union_candidates(master, validation, timepoints)
    intervention = fit_union_total_response(
        mvp, Path(inputs["pseudobulk"]), candidates
    )
    intervention = add_timepoint_observational_direction(
        intervention,
        timepoints,
        association_prefix="signed_association_",
        output_prefix="observational_direction_concordant_",
    )
    state_candidates = prepare_union_candidates(
        master,
        validation,
        timepoints,
        selection_prefix="state_unadjusted_E0_timepoint_",
        association_prefix="state_unadjusted_signed_association_",
    )
    state_intervention = fit_union_total_response(
        mvp, Path(inputs["pseudobulk"]), state_candidates
    )
    state_intervention = add_timepoint_observational_direction(
        state_intervention,
        timepoints,
        association_prefix="state_unadjusted_signed_association_",
        output_prefix="state_unadjusted_direction_concordant_",
    )
    frozen_time_support = pd.read_csv(
        "reports/time_resolved_support/strict_time_resolved_total_support.csv.gz"
    )
    summary = summarize_audit(
        master,
        intervention,
        sample_summary,
        timepoints,
        frozen_time_support,
    )
    summary["state_residualization_sensitivity"] = (
        summarize_state_residualization_sensitivity(
            master, state_intervention, timepoints
        )
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    master.to_csv(
        output_dir / "all_candidate_timepoint_discovery.csv.gz",
        index=False,
        compression="gzip",
    )
    master.loc[
        master.any_preterminal_E0 | master.frozen_day14_E0
    ].to_csv(output_dir / "timepoint_E0_union.csv", index=False)
    intervention.to_csv(
        output_dir / "union_strict_total_response_support.csv.gz",
        index=False,
        compression="gzip",
    )
    intervention.loc[
        intervention[[f"strict_total_support_{timepoint}" for timepoint in timepoints]].any(axis=1)
    ].to_csv(output_dir / "union_supported_edges.csv", index=False)
    state_intervention.to_csv(
        output_dir / "state_unadjusted_union_strict_total_response_support.csv.gz",
        index=False,
        compression="gzip",
    )
    state_intervention.loc[
        state_intervention[
            [f"strict_total_support_{timepoint}" for timepoint in timepoints]
        ].any(axis=1)
    ].to_csv(
        output_dir / "state_unadjusted_union_supported_edges.csv", index=False
    )
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/discovery_selection_bias.yaml"),
    )
    parser.add_argument(
        "--mvp-config", type=Path, default=Path("config/mvp.yaml")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/discovery_selection_bias"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.mvp_config, args.output), indent=2))


if __name__ == "__main__":
    main()
