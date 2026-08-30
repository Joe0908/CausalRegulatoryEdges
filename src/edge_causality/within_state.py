"""Secondary perturbation-by-erythroid-state analysis.

Cell-state labels are observed after perturbation, so this model is explicitly a
secondary diagnostic and not a replacement for the primary exogenous-timepoint
interaction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import classify_edge, fit_interaction_model


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def state_interaction_design(
    groups: pd.DataFrame,
    target_guides: list[str],
    states: list[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    condition = groups.guide.isin(target_guides).astype(float).to_numpy()
    replicate = pd.get_dummies(
        groups.replicate.astype(str), drop_first=True
    ).to_numpy(dtype=float)
    state_values = pd.Categorical(groups.cell_type, categories=states, ordered=True)
    state_dummies = pd.get_dummies(state_values, drop_first=True).to_numpy(dtype=float)
    reduced = np.column_stack(
        [np.ones(len(groups)), condition, replicate, state_dummies]
    )
    interactions = [condition * state_dummies[:, index] for index in range(state_dummies.shape[1])]
    full = np.column_stack([reduced, *interactions])
    interaction_start = reduced.shape[1]
    contrasts = {}
    for state_index, state in enumerate(states):
        contrast = np.zeros(full.shape[1], dtype=float)
        contrast[1] = 1.0
        if state_index > 0:
            contrast[interaction_start + state_index - 1] = 1.0
        contrasts[state] = contrast
    return reduced, full, contrasts


def guide_state_effects(
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    indices: np.ndarray,
    target_guides: list[str],
    reference_guides: list[str],
    state: str,
) -> np.ndarray:
    output = []
    for guide in target_guides:
        differences = []
        rows = groups.index[groups.guide.eq(guide) & groups.cell_type.eq(state)]
        for row in rows:
            reference = groups.index[
                groups.guide.isin(reference_guides)
                & groups.replicate.eq(groups.loc[row, "replicate"])
                & groups.cell_type.eq(state)
            ]
            if len(reference):
                differences.append(
                    log_cpm[row, indices] - log_cpm[reference][:, indices].mean(axis=0)
                )
        output.append(
            np.mean(differences, axis=0)
            if differences
            else np.full(len(indices), np.nan)
        )
    return np.vstack(output)


def run(
    config_path: Path,
    pseudobulk_dir: Path,
    primary_edge_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["state_dependence"]
    # BFU-E lacks AAVS1 pseudobulks above the frozen 15-cell threshold, so it is
    # retained in the raw aggregation but excluded prospectively from this fit.
    requested_states = list(settings["secondary_erythroid_states"])
    reference_guides = list(config["data"]["intervention_reference_guides"])
    minimum_cells = int(settings["minimum_cells_per_pseudobulk"])

    z = np.load(pseudobulk_dir / "erythroid_state_candidate_counts.npz")
    counts = z["counts"].astype(np.float64)
    keys = z["feature_key"].astype(str)
    lookup = {key: index for index, key in enumerate(keys)}
    raw_groups = pd.read_csv(pseudobulk_dir / "erythroid_state_candidate_groups.csv")
    keep = raw_groups.n_cells.to_numpy() >= minimum_cells
    counts = counts[keep]
    groups = raw_groups.loc[keep].reset_index(drop=True)
    library_total = counts.sum(axis=1, keepdims=True)
    log_cpm = np.log2(counts / np.maximum(library_total, 1) * 1_000_000 + 0.5)

    state_coverage = (
        groups.loc[groups.guide.isin(reference_guides)]
        .groupby("cell_type", observed=True)
        .agg(reference_guides=("guide", "nunique"), reference_reps=("replicate", "nunique"))
    )
    states = [
        state
        for state in requested_states
        if state in state_coverage.index
        and state_coverage.loc[state, "reference_guides"] >= 2
        and state_coverage.loc[state, "reference_reps"] >= 2
    ]

    primary = pd.read_csv(primary_edge_path)
    blocks = []
    for tf, edges in primary.groupby("TF", observed=True):
        target_guides = sorted(
            groups.loc[groups.guide_target.eq(tf), "guide"].unique().tolist()
        )
        selected = groups.index[
            groups.guide.isin(target_guides + reference_guides)
            & groups.cell_type.isin(states)
        ].to_numpy()
        selected_groups = groups.loc[selected].reset_index(drop=True)
        reduced, full, contrasts = state_interaction_design(
            selected_groups, target_guides, states
        )
        indices = edges.target.map(lookup).to_numpy(dtype=int)
        weights = selected_groups.n_cells.to_numpy(dtype=float)
        weights /= np.median(weights)
        model = fit_interaction_model(
            log_cpm[selected][:, indices], reduced, full, weights, contrasts
        )
        block = edges[["TF", "target"]].copy().reset_index(drop=True)
        for name, value in model.items():
            block[name] = value
        for state in states:
            guide_matrix = guide_state_effects(
                log_cpm,
                groups,
                indices,
                target_guides,
                reference_guides,
                state,
            )
            pooled = block[f"effect_{state}"].to_numpy()
            block[f"consistent_guides_{state}"] = np.sum(
                np.sign(guide_matrix) == np.sign(pooled)[None, :], axis=0
            )
        blocks.append(block)

    results = pd.concat(blocks, ignore_index=True)
    results["interaction_fdr"] = bh_adjust(results.interaction_p_value.to_numpy())
    for state in states:
        results[f"effect_fdr_{state}"] = bh_adjust(
            results[f"effect_p_value_{state}"].to_numpy()
        )
    results["within_state_edge_class"] = results.apply(
        classify_edge, axis=1, settings=settings, timepoints=states
    )
    results["within_state_heterogeneous"] = results.interaction_fdr < float(
        settings["interaction_fdr_max"]
    )
    rename = {
        column: f"within_state_{column}"
        for column in results.columns
        if column not in {"TF", "target", "within_state_edge_class", "within_state_heterogeneous"}
    }
    results = results.rename(columns=rename)
    combined = primary.merge(results, on=["TF", "target"], validate="one_to_one")

    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir / "within_erythroid_state_interactions.csv.gz",
        index=False,
        compression="gzip",
    )
    combined.to_csv(
        output_dir / "timepoint_and_within_state_edges.csv.gz",
        index=False,
        compression="gzip",
    )
    primary_state = combined.state_dependent
    fate = combined.fate_shift_associated
    within = combined.within_state_heterogeneous
    summary = {
        "analysis_role": "secondary_post_treatment_state_diagnostic",
        "states_requested": requested_states,
        "states_modeled": states,
        "excluded_states": [state for state in requested_states if state not in states],
        "minimum_cells_per_profile": minimum_cells,
        "E0_edges_tested": int(len(combined)),
        "within_state_interaction_fdr_lt_0_05": int(within.sum()),
        "primary_timepoint_dependent_edges": int(primary_state.sum()),
        "primary_edges_also_within_state_heterogeneous": int((primary_state & within).sum()),
        "fate_shift_associated_primary_edges": int((primary_state & fate).sum()),
        "fate_shift_associated_also_within_state_heterogeneous": int(
            (primary_state & fate & within).sum()
        ),
        "non_fate_associated_primary_edges": int((primary_state & ~fate).sum()),
        "non_fate_associated_also_within_state_heterogeneous": int(
            (primary_state & ~fate & within).sum()
        ),
        "primary_and_within_state_candidates": combined.loc[
            primary_state & within, ["TF", "target"]
        ].to_dict(orient="records"),
        "within_state_classes": {
            str(key): int(value)
            for key, value in combined.within_state_edge_class.value_counts().items()
        },
        "caution": (
            "Author erythroid labels are measured after perturbation; these results diagnose persistence within labels but do not by themselves identify a causal mediator."
        ),
    }
    with (output_dir / "within_state_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--pseudobulk", type=Path, default=Path("data/processed/state_pseudobulk")
    )
    parser.add_argument(
        "--primary-edges",
        type=Path,
        default=Path(
            "reports/state_dependence/timepoint_edges_with_trajectory_diagnostic.csv.gz"
        ),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/state_dependence")
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.config, args.pseudobulk, args.primary_edges, args.output),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
