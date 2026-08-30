"""Preterminal transportability of frozen day-14 observational edges.

The confirmatory predictors in this module use controls from days 7, 9 and 11.
Day 14 supplies only the already-frozen E0 membership and association direction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import chi2, fisher_exact, spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import yaml


KEY_COLUMNS = ["TF", "target"]


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def timepoint_slug(timepoint: str) -> str:
    return str(timepoint).replace(" ", "")


def design_matrix(metadata: pd.DataFrame, library_sizes: np.ndarray) -> np.ndarray:
    categorical = pd.get_dummies(
        metadata[["replicate", "new_CellType"]].astype(str), drop_first=True
    ).to_numpy(dtype=float)
    log_depth = np.log1p(library_sizes).reshape(-1, 1)
    depth_scale = max(float(log_depth.std()), 1e-8)
    depth_z = (log_depth - log_depth.mean()) / depth_scale
    return np.column_stack([np.ones(len(metadata)), depth_z, categorical])


def standardized_residuals(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    q, _ = np.linalg.qr(design, mode="reduced")
    output = values - q @ (q.T @ values)
    output -= output.mean(axis=0, keepdims=True)
    scale = output.std(axis=0, ddof=1, keepdims=True)
    scale[scale < 1e-8] = np.nan
    return output / scale


def log_normalize_selected(
    counts: sparse.csr_matrix, feature_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    factors = np.divide(
        10_000.0,
        library_sizes,
        out=np.zeros_like(library_sizes),
        where=library_sizes > 0,
    )
    selected = counts[:, feature_indices].astype(np.float32)
    selected = selected.multiply(factors[:, None]).tocsr()
    selected.data = np.log1p(selected.data)
    return selected.toarray().astype(np.float32), library_sizes


def correlation_matrix(
    values: np.ndarray,
    metadata: pd.DataFrame,
    library_sizes: np.ndarray,
    tf_positions: np.ndarray,
) -> np.ndarray:
    residuals = standardized_residuals(
        values, design_matrix(metadata, library_sizes)
    )
    correlations = residuals[:, tf_positions].T @ residuals / max(len(values) - 1, 1)
    return np.clip(correlations, -0.999999, 0.999999)


def edge_block_from_matrix(
    edge_pairs: pd.DataFrame,
    correlations: np.ndarray,
    tf_names: list[str],
    target_positions: dict[str, int],
) -> pd.DataFrame:
    tf_positions = {tf: i for i, tf in enumerate(tf_names)}
    output = edge_pairs.copy()
    output["signed_association"] = [
        correlations[tf_positions[tf], target_positions[target]]
        for tf, target in output[KEY_COLUMNS].itertuples(index=False, name=None)
    ]
    output["absolute_association"] = output.signed_association.abs()
    output["rank_fraction"] = output.groupby(
        "TF", observed=True
    ).absolute_association.rank(method="average", ascending=False, pct=True)
    return output


def estimate_timepoint(
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    feature_indices: np.ndarray,
    fixed_targets: pd.DataFrame,
    edge_pairs: pd.DataFrame,
    tf_names: list[str],
) -> tuple[pd.DataFrame, dict]:
    values, library_sizes = log_normalize_selected(counts, feature_indices)
    target_positions = {
        target: i for i, target in enumerate(fixed_targets.target.astype(str))
    }
    missing_tfs = [tf for tf in tf_names if tf not in target_positions]
    if missing_tfs:
        raise ValueError(f"TFs absent from fixed target universe: {missing_tfs}")
    tf_positions = np.array([target_positions[tf] for tf in tf_names], dtype=int)
    pooled = correlation_matrix(values, metadata, library_sizes, tf_positions)
    output = edge_block_from_matrix(
        edge_pairs, pooled, tf_names, target_positions
    )

    loo_values = []
    replicates = sorted(metadata.replicate.astype(str).unique())
    for omitted in replicates:
        keep = ~metadata.replicate.astype(str).eq(omitted).to_numpy()
        loo = correlation_matrix(
            values[keep], metadata.loc[keep].reset_index(drop=True), library_sizes[keep], tf_positions
        )
        loo_block = edge_block_from_matrix(
            edge_pairs, loo, tf_names, target_positions
        )
        loo_values.append(loo_block.signed_association.to_numpy(dtype=float))
    loo_matrix = np.vstack(loo_values)
    pooled_sign = np.sign(output.signed_association.to_numpy(dtype=float))
    output["loo_sign_fraction"] = np.mean(
        np.sign(loo_matrix) == pooled_sign[None, :], axis=0
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        output["loo_fisher_z_sd"] = np.nanstd(
            np.arctanh(np.clip(loo_matrix, -0.999999, 0.999999)), axis=0, ddof=1
        )

    raw_selected = counts[:, feature_indices]
    detection = np.asarray((raw_selected > 0).mean(axis=0)).ravel()
    cpm_factors = np.divide(
        1_000_000.0,
        library_sizes,
        out=np.zeros_like(library_sizes),
        where=library_sizes > 0,
    )
    mean_cpm = np.asarray(
        raw_selected.multiply(cpm_factors[:, None]).mean(axis=0)
    ).ravel()
    target_detection = dict(zip(fixed_targets.target.astype(str), detection))
    target_mean_cpm = dict(zip(fixed_targets.target.astype(str), mean_cpm))
    output["detection_fraction"] = output.target.map(target_detection)
    output["mean_cpm"] = output.target.map(target_mean_cpm)

    summary = {
        "cells": int(len(metadata)),
        "libraries": int(len(replicates)),
        "replicate_counts": {
            str(k): int(v) for k, v in metadata.replicate.value_counts().items()
        },
        "cell_type_counts": {
            str(k): int(v) for k, v in metadata.new_CellType.value_counts().items()
        },
        "estimable_edges": int(output.signed_association.notna().sum()),
        "correlation_fisher_z_weight": int(
            max(
                len(metadata)
                - np.linalg.matrix_rank(design_matrix(metadata, library_sizes))
                - 3,
                1,
            )
        ),
    }
    return output, summary


def build_transport_features(
    day14_edges: pd.DataFrame,
    timepoint_edges: dict[str, pd.DataFrame],
    e1_matrix: pd.DataFrame,
    perturbation_trajectory: pd.DataFrame,
    settings: dict,
    timepoint_fisher_z_weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    base_columns = [
        "TF",
        "target",
        "target_symbol",
        "target_gene_id",
        "signed_association",
        "absolute_association",
        "stable_edge",
        "detection_fraction",
        "mean_cpm",
        "residual_rank_fraction",
        "author_TF_sensitive",
        "author_log2fc",
        "author_beta_weight",
        "author_p_weight",
        "author_effect_025",
        "author_direction_concordant",
        "author_supported_concordant",
    ]
    available_base_columns = [
        column for column in base_columns if column in day14_edges.columns
    ]
    frozen = day14_edges.loc[
        day14_edges.stable_edge, available_base_columns
    ].copy()
    frozen = frozen.rename(
        columns={
            "signed_association": "signed_association_day14",
            "absolute_association": "absolute_association_day14",
            "residual_rank_fraction": "rank_fraction_day14",
        }
    )
    # Preserve the frozen config's baseline feature name while retaining the
    # explicit day-14 label used by the trajectory columns.
    frozen["absolute_association"] = frozen.absolute_association_day14
    wide = frozen
    timepoints = list(settings["estimand"]["confirmatory_timepoints"])
    for timepoint in timepoints:
        slug = timepoint_slug(timepoint)
        block = timepoint_edges[timepoint][
            KEY_COLUMNS
            + [
                "signed_association",
                "absolute_association",
                "rank_fraction",
                "loo_sign_fraction",
                "loo_fisher_z_sd",
                "detection_fraction",
                "mean_cpm",
            ]
        ].rename(
            columns={
                column: f"{column}_{slug}"
                for column in [
                    "signed_association",
                    "absolute_association",
                    "rank_fraction",
                    "loo_sign_fraction",
                    "loo_fisher_z_sd",
                    "detection_fraction",
                    "mean_cpm",
                ]
            }
        )
        wide = wide.merge(block, on=KEY_COLUMNS, validate="one_to_one")

    association_columns = [
        f"signed_association_{timepoint_slug(timepoint)}" for timepoint in timepoints
    ]
    rank_columns = [
        f"rank_fraction_{timepoint_slug(timepoint)}" for timepoint in timepoints
    ]
    loo_columns = [
        f"loo_sign_fraction_{timepoint_slug(timepoint)}" for timepoint in timepoints
    ]
    association = wide[association_columns].to_numpy(dtype=float)
    day14_sign = np.sign(wide.signed_association_day14.to_numpy(dtype=float))
    signs = np.sign(association)
    finite = np.isfinite(association)
    sign_matches = (signs == day14_sign[:, None]) & finite
    wide["estimable_prior_timepoints"] = finite.sum(axis=1)
    wide["prior_sign_concordance_fraction"] = np.divide(
        sign_matches.sum(axis=1),
        finite.sum(axis=1),
        out=np.full(len(wide), np.nan),
        where=finite.sum(axis=1) > 0,
    )
    wide["prior_median_rank_fraction"] = np.nanmedian(
        wide[rank_columns].to_numpy(dtype=float), axis=1
    )
    wide["prior_median_rank_score"] = 1 - wide.prior_median_rank_fraction
    fisher_z = np.arctanh(np.clip(association, -0.999999, 0.999999))
    wide["prior_median_aligned_fisher_z"] = np.nanmedian(
        fisher_z * day14_sign[:, None], axis=1
    )
    wide["prior_fisher_z_sd"] = np.nanstd(fisher_z, axis=1, ddof=1)
    wide["prior_association_sd"] = np.nanstd(association, axis=1, ddof=1)
    wide["prior_association_range"] = np.nanmax(
        association, axis=1
    ) - np.nanmin(association, axis=1)
    wide["prior_rank_sd"] = np.nanstd(
        wide[rank_columns].to_numpy(dtype=float), axis=1, ddof=1
    )
    wide["minimum_loo_sign_fraction"] = np.nanmin(
        wide[loo_columns].to_numpy(dtype=float), axis=1
    )

    wide["prior_association_heterogeneity_q"] = np.nan
    wide["prior_association_heterogeneity_p"] = np.nan
    if timepoint_fisher_z_weights is not None:
        weights = np.array(
            [timepoint_fisher_z_weights[timepoint] for timepoint in timepoints],
            dtype=float,
        )
        for row_index in range(len(wide)):
            valid = np.isfinite(fisher_z[row_index]) & (weights > 0)
            if valid.sum() < 2:
                continue
            local_z = fisher_z[row_index, valid]
            local_weights = weights[valid]
            mean_z = np.average(local_z, weights=local_weights)
            statistic = float(np.sum(local_weights * (local_z - mean_z) ** 2))
            wide.loc[row_index, "prior_association_heterogeneity_q"] = statistic
            wide.loc[row_index, "prior_association_heterogeneity_p"] = float(
                chi2.sf(statistic, df=int(valid.sum() - 1))
            )
    wide["prior_association_heterogeneity_fdr"] = benjamini_hochberg(
        wide.prior_association_heterogeneity_p.to_numpy(dtype=float)
    )

    rule = settings["transportability"]["descriptive_transportable_rule"]
    wide["observationally_transportable"] = (
        (wide.estimable_prior_timepoints == len(timepoints))
        & (wide.prior_sign_concordance_fraction == 1.0)
        & (
            wide.prior_median_rank_fraction
            <= float(rule["maximum_median_rank_fraction"])
        )
        & (
            wide.minimum_loo_sign_fraction
            >= float(rule["minimum_leave_one_library_out_sign_fraction"])
        )
    )

    e1_columns = [
        "TF",
        "target",
        "effective_guides_used",
        "perturbation_log2fc",
        "perturbation_fdr",
        "E1_supported",
        "E1_direction_concordant",
    ]
    wide = wide.merge(
        e1_matrix[e1_columns], on=KEY_COLUMNS, how="left", validate="one_to_one"
    )
    wide["strict_E1_testable"] = wide.effective_guides_used >= 2

    trajectory_columns = KEY_COLUMNS + [
        f"effect_{timepoint}" for timepoint in ["day 7", "day 9", "day 11", "day 14"]
    ]
    trajectory_columns += [
        column
        for column in [
            "interaction_fdr",
            "state_dependent",
            "trajectory_shift_supported",
            "fate_shift_associated",
        ]
        if column in perturbation_trajectory.columns
    ]
    wide = wide.merge(
        perturbation_trajectory[trajectory_columns],
        on=KEY_COLUMNS,
        how="left",
        validate="one_to_one",
    )
    observational_columns = association_columns + ["signed_association_day14"]
    perturbation_columns = [
        "effect_day 7",
        "effect_day 9",
        "effect_day 11",
        "effect_day 14",
    ]
    magnitude_correlations = []
    oriented_correlations = []
    direction_concordance = []
    for _, row in wide.iterrows():
        observational = row[observational_columns].to_numpy(dtype=float)
        perturbation = row[perturbation_columns].to_numpy(dtype=float)
        if (
            np.isfinite(observational).all()
            and np.isfinite(perturbation).all()
            and np.ptp(np.abs(observational)) > 0
            and np.ptp(np.abs(perturbation)) > 0
        ):
            magnitude_correlations.append(
                float(spearmanr(np.abs(observational), np.abs(perturbation)).statistic)
            )
        else:
            magnitude_correlations.append(np.nan)
        oriented_observational = np.sign(row.signed_association_day14) * observational
        oriented_perturbation = -np.sign(row.signed_association_day14) * perturbation
        if (
            np.isfinite(oriented_observational).all()
            and np.isfinite(oriented_perturbation).all()
            and np.ptp(oriented_observational) > 0
            and np.ptp(oriented_perturbation) > 0
        ):
            oriented_correlations.append(
                float(spearmanr(oriented_observational, oriented_perturbation).statistic)
            )
        else:
            oriented_correlations.append(np.nan)
        direction_valid = (
            np.isfinite(observational)
            & np.isfinite(perturbation)
            & (observational != 0)
            & (perturbation != 0)
        )
        direction_concordance.append(
            float(
                np.mean(
                    np.sign(perturbation[direction_valid])
                    == -np.sign(observational[direction_valid])
                )
            )
            if direction_valid.any()
            else np.nan
        )
    wide["observational_perturbation_magnitude_spearman"] = magnitude_correlations
    wide["observational_perturbation_oriented_spearman"] = oriented_correlations
    wide["observational_perturbation_direction_concordance_fraction"] = (
        direction_concordance
    )
    wide["log1p_mean_cpm"] = np.log1p(wide.mean_cpm)

    support = wide.author_effect_025.fillna(False)
    wide["descriptive_edge_class"] = np.select(
        [
            wide.observationally_transportable & support,
            wide.observationally_transportable & ~support,
            ~wide.observationally_transportable & support,
        ],
        [
            "transportable_and_atlas_supported",
            "transportable_but_atlas_unsupported",
            "variable_but_atlas_supported",
        ],
        default="variable_and_atlas_unsupported",
    )
    wide["day14_selected_not_preterminal_transportable"] = (
        (wide.prior_sign_concordance_fraction < 1.0)
        | (
            wide.prior_median_rank_fraction
            > float(rule["maximum_median_rank_fraction"])
        )
    )
    if "state_dependent" in wide:
        wide["observationally_variable_and_perturbationally_variable"] = (
            ~wide.observationally_transportable & wide.state_dependent.fillna(False)
        )
    wide["transportable_but_strict_E1_unsupported"] = (
        wide.observationally_transportable
        & wide.strict_E1_testable
        & ~wide.E1_supported.fillna(False)
    )
    direction_threshold = float(
        settings.get("diagnostics", {}).get(
            "minimum_direction_concordance_fraction", 0.5
        )
    )
    wide["observational_perturbation_direction_conflict"] = (
        wide.observational_perturbation_direction_concordance_fraction
        < direction_threshold
    )
    return wide


def model_pipeline(numeric_features: list[str]) -> Pipeline:
    numeric = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ]
    )
    categorical = OneHotEncoder(
        handle_unknown="ignore", drop="first", sparse_output=False
    )
    transform = ColumnTransformer(
        [
            ("numeric", numeric, numeric_features),
            ("TF", categorical, ["TF"]),
        ],
        sparse_threshold=0,
    )
    return Pipeline(
        [
            ("transform", transform),
            (
                "model",
                LogisticRegression(
                    C=1.0, solver="liblinear", max_iter=5000
                ),
            ),
        ]
    )


def group_bootstrap_metric_difference(
    frame: pd.DataFrame,
    metric: str,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    groups = frame.target.astype(str).unique()
    group_rows = {
        group: frame.index[frame.target.astype(str).eq(group)].to_numpy()
        for group in groups
    }
    differences = []
    for _ in range(iterations):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([group_rows[group] for group in sampled])
        sample = frame.loc[indices]
        y = sample.outcome.to_numpy(dtype=int)
        if len(np.unique(y)) < 2:
            continue
        if metric == "AUPRC":
            baseline = average_precision_score(y, sample.baseline_prediction)
            extended = average_precision_score(y, sample.extended_prediction)
        elif metric == "AUROC":
            baseline = roc_auc_score(y, sample.baseline_prediction)
            extended = roc_auc_score(y, sample.extended_prediction)
        else:
            raise ValueError(f"Unknown metric: {metric}")
        differences.append(float(extended - baseline))
    return np.asarray(differences, dtype=float)


def cross_validated_incremental_prediction(
    edges: pd.DataFrame, settings: dict
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prediction = settings["prediction"]
    cv = prediction["cross_validation"]
    baseline_numeric = [
        feature for feature in prediction["baseline_features"] if feature != "TF"
    ]
    extended_numeric = baseline_numeric + list(prediction["extended_features"])
    outcomes = [
        prediction["primary_outcome"],
        "author_supported_concordant",
        "author_TF_sensitive",
        "strict_day14_E1_among_testable_TFs",
    ]
    metric_rows = []
    prediction_rows = []
    random_seed = int(settings["analysis"]["random_seed"])
    for outcome_name in outcomes:
        if outcome_name == "strict_day14_E1_among_testable_TFs":
            frame = edges.loc[edges.strict_E1_testable].copy()
            outcome_column = "E1_supported"
        else:
            frame = edges.copy()
            outcome_column = outcome_name
        frame = frame.loc[frame[outcome_column].notna()].reset_index(drop=True)
        y = frame[outcome_column].astype(int).to_numpy()
        groups = frame.target.astype(str).to_numpy()
        if len(np.unique(y)) < 2:
            continue
        for repeat in range(int(cv["repeats"])):
            splitter = StratifiedGroupKFold(
                n_splits=int(cv["splits"]),
                shuffle=True,
                random_state=random_seed + repeat,
            )
            baseline_prediction = np.full(len(frame), np.nan)
            extended_prediction = np.full(len(frame), np.nan)
            for train, test in splitter.split(frame, y, groups):
                baseline = model_pipeline(baseline_numeric)
                extended = model_pipeline(extended_numeric)
                baseline.fit(frame.iloc[train], y[train])
                extended.fit(frame.iloc[train], y[train])
                baseline_prediction[test] = baseline.predict_proba(frame.iloc[test])[:, 1]
                extended_prediction[test] = extended.predict_proba(frame.iloc[test])[:, 1]
            for model_name, values in [
                ("baseline", baseline_prediction),
                ("extended", extended_prediction),
            ]:
                metric_rows.append(
                    {
                        "outcome": outcome_name,
                        "repeat": repeat,
                        "model": model_name,
                        "positives": int(y.sum()),
                        "edges": int(len(y)),
                        "AUPRC": float(average_precision_score(y, values)),
                        "AUROC": float(roc_auc_score(y, values)),
                    }
                )
            for row, baseline_value, extended_value in zip(
                frame.itertuples(index=False), baseline_prediction, extended_prediction
            ):
                prediction_rows.append(
                    {
                        "outcome_name": outcome_name,
                        "repeat": repeat,
                        "TF": row.TF,
                        "target": row.target,
                        "outcome": int(getattr(row, outcome_column)),
                        "baseline_prediction": float(baseline_value),
                        "extended_prediction": float(extended_value),
                    }
                )

    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    summary_rows = []
    rng = np.random.default_rng(random_seed)
    iterations = int(prediction["target_cluster_bootstrap_iterations"])
    for outcome_name, group in metrics.groupby("outcome", observed=True):
        pivot = group.pivot(index="repeat", columns="model", values=["AUPRC", "AUROC"])
        first_predictions = predictions.loc[
            predictions.outcome_name.eq(outcome_name) & predictions.repeat.eq(0)
        ].reset_index(drop=True)
        for metric in ["AUPRC", "AUROC"]:
            difference = pivot[(metric, "extended")] - pivot[(metric, "baseline")]
            bootstrap = group_bootstrap_metric_difference(
                first_predictions, metric, iterations, rng
            )
            summary_rows.append(
                {
                    "outcome": outcome_name,
                    "metric": metric,
                    "edges": int(group.edges.iloc[0]),
                    "positives": int(group.positives.iloc[0]),
                    "baseline_mean": float(pivot[(metric, "baseline")].mean()),
                    "extended_mean": float(pivot[(metric, "extended")].mean()),
                    "mean_difference_across_repeats": float(difference.mean()),
                    "minimum_difference_across_repeats": float(difference.min()),
                    "maximum_difference_across_repeats": float(difference.max()),
                    "repeat_positive_fraction": float((difference > 0).mean()),
                    "target_bootstrap_difference_q025": float(
                        np.quantile(bootstrap, 0.025)
                    ),
                    "target_bootstrap_difference_q975": float(
                        np.quantile(bootstrap, 0.975)
                    ),
                }
            )
    return metrics, predictions, pd.DataFrame(summary_rows)


def feature_ablation_prediction(edges: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """Post-freeze diagnostic identifying which temporal feature moves ranking."""
    prediction = settings["prediction"]
    cv = prediction["cross_validation"]
    baseline = [
        feature for feature in prediction["baseline_features"] if feature != "TF"
    ]
    temporal = list(prediction["extended_features"])
    specifications = {"baseline": baseline, "extended_all": baseline + temporal}
    specifications.update(
        {f"plus_{feature}": baseline + [feature] for feature in temporal}
    )
    rows = []
    seed = int(settings["analysis"]["random_seed"])
    for outcome in ["author_effect_025", "author_supported_concordant"]:
        y = edges[outcome].astype(int).to_numpy()
        groups = edges.target.astype(str).to_numpy()
        for repeat in range(int(cv["repeats"])):
            splitter = StratifiedGroupKFold(
                n_splits=int(cv["splits"]),
                shuffle=True,
                random_state=seed + repeat,
            )
            predictions = {
                name: np.full(len(edges), np.nan) for name in specifications
            }
            for train, test in splitter.split(edges, y, groups):
                for name, features in specifications.items():
                    model = model_pipeline(features)
                    model.fit(edges.iloc[train], y[train])
                    predictions[name][test] = model.predict_proba(edges.iloc[test])[:, 1]
            for name, values in predictions.items():
                rows.append(
                    {
                        "outcome": outcome,
                        "repeat": repeat,
                        "model": name,
                        "AUPRC": float(average_precision_score(y, values)),
                        "AUROC": float(roc_auc_score(y, values)),
                    }
                )
    detail = pd.DataFrame(rows)
    baseline_means = (
        detail.loc[detail.model.eq("baseline")]
        .groupby("outcome", observed=True)[["AUPRC", "AUROC"]]
        .mean()
    )
    summary = (
        detail.groupby(["outcome", "model"], observed=True)[["AUPRC", "AUROC"]]
        .agg(["mean", "min", "max"])
        .reset_index()
    )
    summary.columns = [
        "_".join(str(value) for value in column if value).rstrip("_")
        if isinstance(column, tuple)
        else str(column)
        for column in summary.columns
    ]
    summary["AUPRC_delta_from_baseline"] = [
        row.AUPRC_mean - baseline_means.loc[row.outcome, "AUPRC"]
        for row in summary.itertuples(index=False)
    ]
    summary["AUROC_delta_from_baseline"] = [
        row.AUROC_mean - baseline_means.loc[row.outcome, "AUROC"]
        for row in summary.itertuples(index=False)
    ]
    return summary


def tf_transportability_enrichment(edges: pd.DataFrame) -> pd.DataFrame:
    """Exact within-TF diagnostic for the prespecified binary display rule."""
    rows = []
    for outcome in ["author_effect_025", "author_supported_concordant", "E1_supported"]:
        frame = edges
        if outcome == "E1_supported":
            frame = edges.loc[edges.strict_E1_testable]
        for tf, group in frame.groupby("TF", observed=True):
            table = (
                pd.crosstab(group.observationally_transportable, group[outcome])
                .reindex(index=[False, True], columns=[False, True], fill_value=0)
            )
            odds_ratio, p_value = fisher_exact(table.to_numpy(), alternative="two-sided")
            rows.append(
                {
                    "outcome": outcome,
                    "TF": str(tf),
                    "variable_edges": int(table.loc[False].sum()),
                    "variable_supported": int(table.loc[False, True]),
                    "transportable_edges": int(table.loc[True].sum()),
                    "transportable_supported": int(table.loc[True, True]),
                    "odds_ratio": float(odds_ratio),
                    "fisher_p_value": float(p_value),
                }
            )
    output = pd.DataFrame(rows)
    output["fisher_fdr_within_outcome"] = np.nan
    for _, index in output.groupby("outcome", observed=True).groups.items():
        index = np.asarray(list(index), dtype=int)
        p = output.loc[index, "fisher_p_value"].to_numpy(dtype=float)
        order = np.argsort(p)
        adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        restored = np.empty_like(adjusted)
        restored[order] = np.clip(adjusted, 0, 1)
        output.loc[index, "fisher_fdr_within_outcome"] = restored
    return output


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment that preserves missing values."""
    p_values = np.asarray(p_values, dtype=float)
    output = np.full(len(p_values), np.nan)
    finite = np.flatnonzero(np.isfinite(p_values))
    if not len(finite):
        return output
    p = p_values[finite]
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    output[finite] = restored
    return output


def summarize(
    edges: pd.DataFrame,
    sample_summary: dict,
    prediction_summary: pd.DataFrame,
    tf_enrichment: pd.DataFrame,
) -> dict:
    class_counts = {
        str(k): int(v) for k, v in edges.descriptive_edge_class.value_counts().items()
    }
    transport_rates = {}
    for outcome in [
        "author_effect_025",
        "author_supported_concordant",
        "author_TF_sensitive",
        "E1_supported",
    ]:
        block = edges
        if outcome == "E1_supported":
            block = edges.loc[edges.strict_E1_testable]
        transport_rates[outcome] = {
            str(bool(k)): {
                "edges": int(len(group)),
                "supported": int(group[outcome].sum()),
                "support_rate": float(group[outcome].mean()),
            }
            for k, group in block.groupby("observationally_transportable")
        }
    return {
        "design_status": "frozen_before_analysis",
        "confirmatory_timepoints": list(sample_summary),
        "sample_summary": sample_summary,
        "frozen_E0_edges": int(len(edges)),
        "observationally_transportable_edges": int(
            edges.observationally_transportable.sum()
        ),
        "preterminal_association_heterogeneity_fdr_005_edges": int(
            (edges.prior_association_heterogeneity_fdr < 0.05).sum()
        ),
        "diagnostic_edge_flags": {
            column: int(edges[column].fillna(False).sum())
            for column in [
                "day14_selected_not_preterminal_transportable",
                "observationally_variable_and_perturbationally_variable",
                "transportable_but_strict_E1_unsupported",
                "observational_perturbation_direction_conflict",
            ]
            if column in edges
        },
        "descriptive_edge_classes": class_counts,
        "support_rates_by_transportability": transport_rates,
        "median_observational_perturbation_magnitude_spearman": float(
            edges.observational_perturbation_magnitude_spearman.median()
        ),
        "median_observational_perturbation_oriented_spearman": float(
            edges.observational_perturbation_oriented_spearman.median()
        ),
        "incremental_prediction": prediction_summary.to_dict(orient="records"),
        "within_TF_transportability_enrichment": tf_enrichment.to_dict(
            orient="records"
        ),
        "interpretation_boundary": (
            "Transportability is a prediction feature within this atlas, not a "
            "formal invariant-causal test or proof of direct regulation."
        ),
    }


def run(
    config_path: Path,
    output_dir: Path,
) -> dict:
    settings = load_yaml(config_path)
    inputs = settings["inputs"]
    controls = Path(inputs["controls"])
    counts = sparse.load_npz(controls / "rna_counts_cells_by_genes.npz").tocsr()
    metadata = pd.read_csv(controls / "cell_metadata.csv.gz", index_col=0)
    features = pd.read_csv(controls / "gene_features.csv.gz")
    day14_edges = pd.read_csv(inputs["day14_edges"])
    e1_matrix = pd.read_csv(inputs["e1_matrix"])
    perturbation_trajectory = pd.read_csv(inputs["perturbation_trajectory"])

    target_metadata = (
        day14_edges[["target", "target_symbol", "target_gene_id"]]
        .drop_duplicates("target")
        .reset_index(drop=True)
    )
    feature_lookup = dict(
        zip(features.feature_key.astype(str), features.index.to_numpy(dtype=int))
    )
    missing = sorted(set(target_metadata.target.astype(str)) - set(feature_lookup))
    if missing:
        raise ValueError(f"{len(missing)} frozen target features are absent")
    feature_indices = np.array(
        [feature_lookup[target] for target in target_metadata.target.astype(str)],
        dtype=int,
    )
    edge_pairs = day14_edges[
        ["TF", "target", "target_symbol", "target_gene_id"]
    ].copy()
    tf_names = list(day14_edges.TF.drop_duplicates().astype(str))

    allowed_states = set(settings["estimand"]["erythroid_cell_types"])
    allowed_guides = set(settings["estimand"]["control_guides"])
    timepoint_edges = {}
    sample_summary = {}
    for timepoint in settings["estimand"]["confirmatory_timepoints"]:
        selected = (
            metadata.Timepoint.eq(timepoint)
            & metadata.new_CellType.isin(allowed_states)
            & metadata.perturbation_name.isin(allowed_guides)
        )
        selected_metadata = metadata.loc[selected].copy()
        library_counts = selected_metadata.replicate.value_counts()
        retained_libraries = set(
            library_counts.index[
                library_counts >= int(settings["estimand"]["minimum_cells_per_library"])
            ]
        )
        selected &= metadata.replicate.isin(retained_libraries)
        selected_rows = np.flatnonzero(selected.to_numpy())
        selected_metadata = metadata.iloc[selected_rows].reset_index(drop=True)
        block, block_summary = estimate_timepoint(
            counts[selected_rows],
            selected_metadata,
            feature_indices,
            target_metadata,
            edge_pairs,
            tf_names,
        )
        timepoint_edges[timepoint] = block
        sample_summary[timepoint] = block_summary

    wide_all = None
    for timepoint, block in timepoint_edges.items():
        slug = timepoint_slug(timepoint)
        renamed = block.rename(
            columns={
                column: f"{column}_{slug}"
                for column in block.columns
                if column not in KEY_COLUMNS + ["target_symbol", "target_gene_id"]
            }
        )
        keep_columns = KEY_COLUMNS + [
            column
            for column in renamed.columns
            if column not in KEY_COLUMNS + ["target_symbol", "target_gene_id"]
        ]
        wide_all = (
            renamed[keep_columns]
            if wide_all is None
            else wide_all.merge(renamed[keep_columns], on=KEY_COLUMNS, validate="one_to_one")
        )

    fisher_z_weights = {
        timepoint: summary["correlation_fisher_z_weight"]
        for timepoint, summary in sample_summary.items()
    }
    frozen = build_transport_features(
        day14_edges,
        timepoint_edges,
        e1_matrix,
        perturbation_trajectory,
        settings,
        fisher_z_weights,
    )
    metrics, predictions, prediction_summary = cross_validated_incremental_prediction(
        frozen, settings
    )
    ablation = feature_ablation_prediction(frozen, settings)
    tf_enrichment = tf_transportability_enrichment(frozen)
    summary = summarize(frozen, sample_summary, prediction_summary, tf_enrichment)

    output_dir.mkdir(parents=True, exist_ok=True)
    wide_all.to_csv(
        output_dir / "all_candidate_preterminal_associations.csv.gz",
        index=False,
        compression="gzip",
    )
    frozen.to_csv(
        output_dir / "frozen_E0_temporal_transportability.csv",
        index=False,
    )
    metrics.to_csv(output_dir / "prediction_metrics_by_repeat.csv", index=False)
    predictions.to_csv(
        output_dir / "prediction_oof.csv.gz", index=False, compression="gzip"
    )
    prediction_summary.to_csv(
        output_dir / "incremental_prediction_summary.csv", index=False
    )
    ablation.to_csv(output_dir / "temporal_feature_ablation.csv", index=False)
    tf_enrichment.to_csv(
        output_dir / "within_TF_transportability_enrichment.csv", index=False
    )
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/temporal_transportability.yaml"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/temporal_transportability")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
