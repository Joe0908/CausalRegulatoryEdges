"""Diagnose perturbation-induced movement on a controls-only temporal axis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from scipy.stats import fisher_exact
import yaml

from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import (
    fit_interaction_model,
    interaction_design,
)


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def numeric_timepoint(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.extract(r"(\d+)")[0].astype(float).to_numpy()


def fit_control_trajectory(
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    eligible: np.ndarray,
    discovery_guides: list[str],
    variable_genes: int,
    principal_components: int,
    ridge_alphas: list[float],
    random_seed: int,
) -> tuple[np.ndarray, dict]:
    """Fit and project a temporal axis using discovery controls only."""
    control_rows = groups.guide.isin(discovery_guides).to_numpy()
    control_expression = log_cpm[control_rows]
    eligible_indices = np.flatnonzero(eligible)
    variance = control_expression[:, eligible_indices].var(axis=0)
    selected = eligible_indices[
        np.argsort(variance)[-min(variable_genes, len(eligible_indices)) :]
    ]
    scaler = StandardScaler()
    control_scaled = scaler.fit_transform(control_expression[:, selected])
    n_components = min(
        principal_components, control_scaled.shape[0] - 1, control_scaled.shape[1]
    )
    pca = PCA(n_components=n_components, random_state=random_seed)
    control_pc = pca.fit_transform(control_scaled)
    y = numeric_timepoint(groups.loc[control_rows, "timepoint"])
    replicate = groups.loc[control_rows, "replicate"].astype(str).to_numpy()

    # Leave one entire library out so the diagnostic cannot rely on a library-
    # specific signature to recover collection time.
    prediction = np.zeros_like(y)
    splitter = LeaveOneGroupOut()
    for train, test in splitter.split(control_pc, y, replicate):
        model = RidgeCV(alphas=np.asarray(ridge_alphas)).fit(
            control_pc[train], y[train]
        )
        prediction[test] = model.predict(control_pc[test])

    final_model = RidgeCV(alphas=np.asarray(ridge_alphas)).fit(control_pc, y)
    all_pc = pca.transform(scaler.transform(log_cpm[:, selected]))
    score = final_model.predict(all_pc)
    diagnostics = {
        "control_profiles": int(control_rows.sum()),
        "selected_variable_genes": int(len(selected)),
        "principal_components": int(n_components),
        "explained_variance_ratio_sum": float(pca.explained_variance_ratio_.sum()),
        "leave_one_library_out_r2": float(r2_score(y, prediction)),
        "leave_one_library_out_mae_days": float(mean_absolute_error(y, prediction)),
        "final_ridge_alpha": float(final_model.alpha_),
    }
    return score, diagnostics


def guide_score_differences(
    scored_groups: pd.DataFrame,
    target_guides: list[str],
    reference_guides: list[str],
    timepoint: str,
) -> np.ndarray:
    output = []
    for guide in target_guides:
        differences = []
        rows = scored_groups.index[
            scored_groups.guide.eq(guide)
            & scored_groups.timepoint.eq(timepoint)
        ]
        for row in rows:
            reference = scored_groups.loc[
                scored_groups.guide.isin(reference_guides)
                & scored_groups.replicate.eq(scored_groups.loc[row, "replicate"]),
                "trajectory_score",
            ]
            if len(reference):
                differences.append(
                    scored_groups.loc[row, "trajectory_score"] - reference.mean()
                )
        output.append(np.mean(differences) if differences else np.nan)
    return np.asarray(output)


def estimate_tf_shifts(
    scored_groups: pd.DataFrame,
    tfs: list[str],
    reference_guides: list[str],
    timepoints: list[str],
    minimum_cells: int,
) -> pd.DataFrame:
    rows = []
    for tf in tfs:
        target_guides = sorted(
            scored_groups.loc[scored_groups.target.eq(tf), "guide"].unique().tolist()
        )
        selected = scored_groups.index[
            scored_groups.guide.isin(target_guides + reference_guides)
            & (scored_groups.n_cells >= minimum_cells)
        ].to_numpy()
        selected_groups = scored_groups.loc[selected].reset_index(drop=True)
        reduced, full, contrasts, _ = interaction_design(
            selected_groups, target_guides, timepoints
        )
        weights = selected_groups.n_cells.to_numpy(dtype=float)
        weights /= np.median(weights)
        model = fit_interaction_model(
            selected_groups.trajectory_score.to_numpy()[:, None],
            reduced,
            full,
            weights,
            contrasts,
        )
        for timepoint in timepoints:
            guide_difference = guide_score_differences(
                scored_groups,
                target_guides,
                reference_guides,
                timepoint,
            )
            pooled = float(model[f"effect_{timepoint}"][0])
            rows.append(
                {
                    "TF": tf,
                    "timepoint": timepoint,
                    "trajectory_shift_days": pooled,
                    "trajectory_shift_se": float(model[f"effect_se_{timepoint}"][0]),
                    "trajectory_shift_p_value": float(
                        model[f"effect_p_value_{timepoint}"][0]
                    ),
                    "interaction_p_value": float(model["interaction_p_value"][0]),
                    "consistent_guides": int(
                        np.sum(np.sign(guide_difference) == np.sign(pooled))
                    ),
                    "minimum_guide_shift": float(np.nanmin(np.abs(guide_difference))),
                    "targeting_guides": ";".join(target_guides),
                }
            )
    output = pd.DataFrame(rows)
    output["trajectory_shift_fdr"] = bh_adjust(
        output.trajectory_shift_p_value.to_numpy()
    )
    interaction = (
        output[["TF", "interaction_p_value"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    interaction["interaction_fdr"] = bh_adjust(interaction.interaction_p_value.to_numpy())
    return output.merge(interaction[["TF", "interaction_fdr"]], on="TF")


def annotate_edges(
    edges: pd.DataFrame,
    shifts: pd.DataFrame,
    timepoints: list[str],
) -> pd.DataFrame:
    output = edges.copy()
    effect_matrix = output[[f"effect_{timepoint}" for timepoint in timepoints]].to_numpy()
    strongest = np.argmax(np.abs(effect_matrix), axis=1)
    output["strongest_effect_timepoint"] = [timepoints[index] for index in strongest]
    shift_columns = shifts[
        [
            "TF",
            "timepoint",
            "trajectory_shift_days",
            "trajectory_shift_fdr",
            "consistent_guides",
            "trajectory_shift_supported",
        ]
    ].rename(columns={"timepoint": "strongest_effect_timepoint"})
    output = output.merge(
        shift_columns,
        on=["TF", "strongest_effect_timepoint"],
        how="left",
        validate="many_to_one",
    )
    output["fate_shift_associated"] = (
        output.state_dependent & output.trajectory_shift_supported.fillna(False)
    )
    return output


def run(
    config_path: Path,
    pseudobulk_dir: Path,
    feature_path: Path,
    edge_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    trajectory_settings = config["control_trajectory"]
    state_settings = config["state_dependence"]
    timepoints = list(state_settings["ordered_timepoints"])
    reference_guides = list(config["data"]["intervention_reference_guides"])

    z = np.load(pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_counts.npz")
    counts = z["counts"].astype(np.float64)
    groups = pd.read_csv(
        pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_groups.csv"
    )
    features = pd.read_csv(feature_path)
    if len(features) != counts.shape[1]:
        raise ValueError("Control feature table and pseudobulk matrix do not align")
    library_total = counts.sum(axis=1, keepdims=True)
    log_cpm = np.log2(counts / np.maximum(library_total, 1) * 1_000_000 + 0.5)
    score, diagnostics = fit_control_trajectory(
        log_cpm,
        groups,
        features.candidate_eligible.to_numpy(dtype=bool),
        list(trajectory_settings["discovery_guides"]),
        int(trajectory_settings["variable_genes"]),
        int(trajectory_settings["principal_components"]),
        [float(value) for value in trajectory_settings["ridge_alphas"]],
        int(config["project"]["random_seed"]),
    )
    scored = groups.copy()
    scored["trajectory_score"] = score

    null_differences = []
    for timepoint in timepoints:
        null_differences.extend(
            guide_score_differences(
                scored,
                list(trajectory_settings["discovery_guides"]),
                reference_guides,
                timepoint,
            ).tolist()
        )
    null_threshold = float(
        np.quantile(
            np.abs(null_differences),
            float(trajectory_settings["null_absolute_quantile"]),
        )
    )
    tfs = list(config["mvp"]["primary_tf_panel"]) + list(
        state_settings["positive_control_tfs"]
    )
    shifts = estimate_tf_shifts(
        scored,
        tfs,
        reference_guides,
        timepoints,
        int(state_settings["minimum_cells_per_pseudobulk"]),
    )
    shifts["trajectory_shift_supported"] = (
        (shifts.trajectory_shift_fdr < 0.05)
        & (shifts.trajectory_shift_days.abs() >= null_threshold)
        & (shifts.consistent_guides >= int(state_settings["minimum_consistent_guides"]))
    )

    edges = pd.read_csv(edge_path)
    annotated = annotate_edges(edges, shifts, timepoints)
    non_fate_state_dependence = (
        annotated.state_dependent & ~annotated.fate_shift_associated
    )
    contingency = pd.crosstab(
        annotated.atlas_supported_day14_unsupported,
        non_fate_state_dependence,
    ).reindex(index=[False, True], columns=[False, True], fill_value=0)
    non_fate_odds_ratio, non_fate_p_value = fisher_exact(contingency.to_numpy())
    output_dir.mkdir(parents=True, exist_ok=True)
    scored.to_csv(
        output_dir / "pseudobulk_control_trajectory_scores.csv.gz",
        index=False,
        compression="gzip",
    )
    shifts.to_csv(output_dir / "tf_trajectory_shifts.csv", index=False)
    annotated.to_csv(
        output_dir / "timepoint_edges_with_trajectory_diagnostic.csv.gz",
        index=False,
        compression="gzip",
    )
    summary = {
        "trajectory_definition": (
            "controls-only RNA temporal axis; perturbed guide pseudobulks projected after fitting"
        ),
        **diagnostics,
        "null_absolute_shift_95_percent_days": null_threshold,
        "supported_tf_timepoint_shifts": int(shifts.trajectory_shift_supported.sum()),
        "state_dependent_edges": int(annotated.state_dependent.sum()),
        "state_dependent_edges_fate_shift_associated": int(
            annotated.fate_shift_associated.sum()
        ),
        "state_dependent_edges_not_fate_shift_associated": int(
            (annotated.state_dependent & ~annotated.fate_shift_associated).sum()
        ),
        "atlas_supported_day14_unsupported_state_dependent": int(
            (
                annotated.atlas_supported_day14_unsupported
                & annotated.state_dependent
            ).sum()
        ),
        "atlas_supported_day14_unsupported_state_dependent_fate_associated": int(
            (
                annotated.atlas_supported_day14_unsupported
                & annotated.state_dependent
                & annotated.fate_shift_associated
            ).sum()
        ),
        "non_fate_state_dependence_enrichment_odds_ratio": float(non_fate_odds_ratio),
        "non_fate_state_dependence_enrichment_fisher_p_value": float(non_fate_p_value),
        "interpretation": (
            "A fate-shift association is a diagnostic flag, not proof that the edge effect is mediated by composition."
        ),
    }
    with (output_dir / "trajectory_shift_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--pseudobulk", type=Path, default=Path("data/processed/pseudobulk")
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/processed/controls_all_timepoints/gene_features.csv.gz"),
    )
    parser.add_argument(
        "--edges",
        type=Path,
        default=Path("reports/state_dependence/timepoint_interaction_edges.csv.gz"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/state_dependence")
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.config, args.pseudobulk, args.features, args.edges, args.output),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
