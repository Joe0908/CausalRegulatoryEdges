"""Audit intervention support at each collection time without redefining E1.

Two deliberately separate modes are provided:

``strict``
    Refit the total-RNA pseudobulk model with the frozen effective guides,
    require complete guide and leave-one-effective-guide-out sign agreement,
    and control FDR across the full edge-by-time family.

``provisional``
    Reuse the committed all-guide timepoint summary when the large pseudobulk
    matrix is unavailable.  These calls are explicitly screening labels and
    are never named E1 because timepoint-specific leave-one-guide-out results
    cannot be reconstructed from the aggregate table.

The original state-matched day-14 E1 result is read only for comparison and is
not modified by either mode.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import (
    feature_keys,
    fit_interaction_model,
    interaction_design,
    load_config,
    per_guide_context_effects,
)


def add_edge_time_fdr(
    table: pd.DataFrame,
    timepoints: list[str],
    p_value_prefix: str = "effect_p_value_",
    output_prefix: str = "effect_fdr_edge_time_",
) -> pd.DataFrame:
    """BH-adjust the single prespecified edge-by-time hypothesis family."""
    output = table.copy()
    locations: list[tuple[str, np.ndarray]] = []
    pooled: list[np.ndarray] = []
    for timepoint in timepoints:
        values = output[f"{p_value_prefix}{timepoint}"].to_numpy(dtype=float)
        finite = np.flatnonzero(np.isfinite(values))
        locations.append((timepoint, finite))
        pooled.append(values[finite])
        output[f"{output_prefix}{timepoint}"] = np.nan
    if not pooled or sum(len(values) for values in pooled) == 0:
        return output
    adjusted = bh_adjust(np.concatenate(pooled))
    offset = 0
    for timepoint, finite in locations:
        width = len(finite)
        output.loc[finite, f"{output_prefix}{timepoint}"] = adjusted[
            offset : offset + width
        ]
        offset += width
    return output


def support_pattern(row: pd.Series, timepoints: list[str], prefix: str) -> str:
    return "".join(
        "1" if bool(row[f"{prefix}{timepoint}"]) else "0"
        for timepoint in timepoints
    )


def classify_temporal_support(
    row: pd.Series,
    timepoints: list[str],
    support_prefix: str,
    testable_column: str,
    interaction_fdr_max: float,
) -> str:
    """Assign cautious support-shape labels without treating non-support as zero."""
    if not bool(row[testable_column]):
        return "guide_efficacy_limited"
    support = np.array(
        [bool(row[f"{support_prefix}{timepoint}"]) for timepoint in timepoints]
    )
    if support.all():
        return "persistent_support_detected"
    if not support.any():
        return "no_timepoint_support_detected"
    if not bool(row.interaction_fdr < interaction_fdr_max):
        return "localized_support_no_detected_heterogeneity"

    indices = np.flatnonzero(support)
    last = len(timepoints) - 1
    if indices.max() <= 1 and not support[last]:
        return "early_window_candidate_with_heterogeneity"
    if indices.min() >= 1 and support[last] and not support[0]:
        return "late_onset_candidate_with_heterogeneity"
    if not support[0] and not support[last] and np.any(support[1:last]):
        return "transient_candidate_with_heterogeneity"
    return "time_varying_support_candidate"


def contrast_name(first: str, second: str) -> str:
    return f"{first}_vs_{second}"


def fit_pairwise_time_contrasts(
    response: np.ndarray,
    full: np.ndarray,
    weights: np.ndarray,
    contrasts: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Estimate every prespecified pairwise collection-time contrast."""
    root_weight = np.sqrt(weights).reshape(-1, 1)
    x = full * root_weight
    y = response * root_weight
    inverse = np.linalg.pinv(x.T @ x)
    beta = inverse @ x.T @ y
    residual = response - full @ beta
    degrees_freedom = max(1, len(response) - np.linalg.matrix_rank(full))
    sigma2 = (weights[:, None] * residual**2).sum(axis=0) / degrees_freedom
    output: dict[str, np.ndarray] = {}
    for first, second in combinations(contrasts, 2):
        name = contrast_name(first, second)
        contrast = contrasts[first] - contrasts[second]
        difference = contrast @ beta
        multiplier = float(contrast @ inverse @ contrast)
        standard_error = np.sqrt(np.maximum(sigma2 * multiplier, 0))
        statistic = np.divide(
            difference,
            standard_error,
            out=np.zeros_like(difference),
            where=standard_error > 0,
        )
        output[f"effect_difference_{name}"] = difference
        output[f"effect_difference_se_{name}"] = standard_error
        output[f"effect_difference_p_value_{name}"] = 2 * student_t.sf(
            np.abs(statistic), degrees_freedom
        )
    return output


def _pair_name(timepoints: list[str], first: str, second: str) -> str:
    first_index = timepoints.index(first)
    second_index = timepoints.index(second)
    if first_index < second_index:
        return contrast_name(first, second)
    return contrast_name(second, first)


def classify_strict_temporal_shape(
    row: pd.Series,
    timepoints: list[str],
    settings: dict,
) -> str:
    """Require contrasts/equivalence before assigning resolved temporal shapes."""
    if not bool(row.strict_time_resolved_testable):
        return "guide_efficacy_limited"
    support = np.array(
        [bool(row[f"strict_total_support_{timepoint}"]) for timepoint in timepoints]
    )
    if support.all():
        return "persistent_support_detected"
    if not support.any():
        return "no_timepoint_support_detected"

    def contrast_supported(first: str, second: str) -> bool:
        name = _pair_name(timepoints, first, second)
        return bool(
            row[f"effect_difference_fdr_edge_pair_{name}"]
            < float(settings["contrast_fdr_max"])
        )

    early = timepoints[:2]
    middle = timepoints[1:-1]
    first = timepoints[0]
    last = timepoints[-1]

    for timepoint in middle:
        if (
            bool(row[f"strict_total_support_{timepoint}"])
            and not support[0]
            and not support[-1]
            and bool(row[f"negligible_effect_supported_{first}"])
            and bool(row[f"negligible_effect_supported_{last}"])
            and contrast_supported(timepoint, first)
            and contrast_supported(timepoint, last)
        ):
            return "transient_support_detected"

    supported_early = [
        timepoint
        for timepoint in early
        if bool(row[f"strict_total_support_{timepoint}"])
    ]
    if supported_early and not support[-1]:
        strongest_early = max(
            supported_early, key=lambda tp: abs(float(row[f"effect_{tp}"]))
        )
        if (
            contrast_supported(strongest_early, last)
            and abs(float(row[f"effect_{strongest_early}"]))
            > abs(float(row[f"effect_{last}"]))
        ):
            if bool(row[f"negligible_effect_supported_{last}"]):
                return "early_specific_support_detected"
            return "attenuating_support_detected"

    if support[-1] and not support[0] and contrast_supported(last, first):
        if abs(float(row[f"effect_{last}"])) > abs(float(row[f"effect_{first}"])):
            if bool(row[f"negligible_effect_supported_{first}"]):
                return "late_onset_support_detected"
            return "late_amplification_support_detected"

    if bool(row.interaction_fdr < float(settings["interaction_fdr_max"])):
        return "time_varying_support_candidate"
    return "localized_support_no_detected_heterogeneity"


def _effective_guides(block: pd.DataFrame) -> list[str]:
    values = block.effective_guide_names.dropna().astype(str).unique().tolist()
    if not values:
        return []
    if len(values) != 1:
        raise ValueError(f"Conflicting effective-guide sets for {block.TF.iloc[0]}")
    return [guide for guide in values[0].split(";") if guide]


def _untestable_block(
    candidates: pd.DataFrame,
    timepoints: list[str],
    reason: str,
    effective_guides: list[str],
) -> pd.DataFrame:
    output = candidates.copy().reset_index(drop=True)
    output["effective_guides_used_time_resolved"] = len(effective_guides)
    output["effective_guide_names_time_resolved"] = ";".join(effective_guides)
    output["strict_time_resolved_testable"] = False
    output["strict_time_resolved_unresolved_reason"] = reason
    for column in [
        "interaction_f_stat",
        "interaction_p_value",
        "model_df_numerator",
        "model_df_denominator",
    ]:
        output[column] = np.nan
    for timepoint in timepoints:
        for stem in [
            "effect_",
            "effect_se_",
            "effect_p_value_",
            "consistent_effective_guides_",
            "effective_guide_direction_consistent_",
            "leave_one_effective_guide_out_direction_consistent_",
        ]:
            output[f"{stem}{timepoint}"] = np.nan
    for first, second in combinations(timepoints, 2):
        name = contrast_name(first, second)
        for stem in [
            "effect_difference_",
            "effect_difference_se_",
            "effect_difference_p_value_",
        ]:
            output[f"{stem}{name}"] = np.nan
    return output


def fit_tf_strict_total_support(
    candidates: pd.DataFrame,
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    effective_guides: list[str],
    reference_guides: list[str],
    timepoints: list[str],
    minimum_cells: int,
) -> pd.DataFrame:
    """Refit one TF using only frozen effective guides and calculate LOO signs."""
    if len(effective_guides) < 2:
        return _untestable_block(
            candidates,
            timepoints,
            "fewer_than_two_effective_guides",
            effective_guides,
        )

    selected = groups.index[
        groups.guide.isin(effective_guides + reference_guides)
        & (groups.n_cells >= minimum_cells)
    ].to_numpy()
    selected_groups = groups.loc[selected].reset_index(drop=True)
    if selected_groups.timepoint.nunique() < len(timepoints):
        return _untestable_block(
            candidates,
            timepoints,
            "insufficient_timepoint_pseudobulk_coverage",
            effective_guides,
        )

    reduced, full, contrasts, _ = interaction_design(
        selected_groups, effective_guides, timepoints
    )
    weights = selected_groups.n_cells.to_numpy(dtype=float)
    weights /= np.median(weights)
    indices = candidates.feature_index.to_numpy(dtype=int)
    model = fit_interaction_model(
        log_cpm[selected][:, indices], reduced, full, weights, contrasts
    )

    output = candidates.copy().reset_index(drop=True)
    output["effective_guides_used_time_resolved"] = len(effective_guides)
    output["effective_guide_names_time_resolved"] = ";".join(effective_guides)
    output["strict_time_resolved_testable"] = True
    output["strict_time_resolved_unresolved_reason"] = "not_applicable"
    for name, values in model.items():
        output[name] = values
    pairwise = fit_pairwise_time_contrasts(
        log_cpm[selected][:, indices], full, weights, contrasts
    )
    for name, values in pairwise.items():
        output[name] = values

    eligible_groups = groups.loc[groups.n_cells >= minimum_cells]
    per_guide = per_guide_context_effects(
        log_cpm,
        eligible_groups,
        indices,
        effective_guides,
        reference_guides,
        timepoints,
    )
    loo_by_timepoint: dict[str, list[np.ndarray]] = {
        timepoint: [] for timepoint in timepoints
    }
    for omitted in effective_guides:
        retained = [guide for guide in effective_guides if guide != omitted]
        loo_rows = groups.index[
            groups.guide.isin(retained + reference_guides)
            & (groups.n_cells >= minimum_cells)
        ].to_numpy()
        loo_groups = groups.loc[loo_rows].reset_index(drop=True)
        loo_reduced, loo_full, loo_contrasts, _ = interaction_design(
            loo_groups, retained, timepoints
        )
        loo_weights = loo_groups.n_cells.to_numpy(dtype=float)
        loo_weights /= np.median(loo_weights)
        loo_model = fit_interaction_model(
            log_cpm[loo_rows][:, indices],
            loo_reduced,
            loo_full,
            loo_weights,
            loo_contrasts,
        )
        for timepoint in timepoints:
            loo_by_timepoint[timepoint].append(loo_model[f"effect_{timepoint}"])

    for timepoint in timepoints:
        effect = output[f"effect_{timepoint}"].to_numpy(dtype=float)
        guide_matrix = per_guide[timepoint]
        same_guide_sign = np.sign(guide_matrix) == np.sign(effect)[None, :]
        output[f"consistent_effective_guides_{timepoint}"] = np.sum(
            same_guide_sign, axis=0
        )
        output[f"effective_guide_direction_consistent_{timepoint}"] = np.all(
            same_guide_sign, axis=0
        )
        loo_matrix = np.vstack(loo_by_timepoint[timepoint])
        output[f"leave_one_effective_guide_out_direction_consistent_{timepoint}"] = (
            np.all(np.sign(loo_matrix) == np.sign(effect)[None, :], axis=0)
        )
    return output


def add_strict_support_calls(
    table: pd.DataFrame,
    settings: dict,
    timepoints: list[str],
) -> pd.DataFrame:
    output = add_edge_time_fdr(table, timepoints)
    tested = output.strict_time_resolved_testable.fillna(False)
    output["interaction_fdr"] = np.nan
    finite_interaction = tested & output.interaction_p_value.notna()
    output.loc[finite_interaction, "interaction_fdr"] = bh_adjust(
        output.loc[finite_interaction, "interaction_p_value"].to_numpy()
    )
    for timepoint in timepoints:
        output[f"effect_fdr_within_timepoint_{timepoint}"] = np.nan
        finite = tested & output[f"effect_p_value_{timepoint}"].notna()
        output.loc[finite, f"effect_fdr_within_timepoint_{timepoint}"] = bh_adjust(
            output.loc[finite, f"effect_p_value_{timepoint}"].to_numpy()
        )
        output[f"strict_total_support_{timepoint}"] = (
            tested
            & (
                output[f"effect_fdr_edge_time_{timepoint}"]
                < float(settings["global_fdr_max"])
            )
            & (
                output[f"effect_{timepoint}"].abs()
                >= float(settings["minimum_absolute_effect"])
            )
            & output[f"effective_guide_direction_consistent_{timepoint}"].fillna(False)
            & output[
                f"leave_one_effective_guide_out_direction_consistent_{timepoint}"
            ].fillna(False)
        )
        output[f"strict_direction_concordant_{timepoint}"] = (
            np.sign(output.signed_association)
            == -np.sign(output[f"effect_{timepoint}"])
        )
        degrees_freedom = output.model_df_denominator.to_numpy(dtype=float)
        critical = student_t.ppf(
            1 - float(settings["equivalence_alpha"]), degrees_freedom
        )
        effect = output[f"effect_{timepoint}"].to_numpy(dtype=float)
        standard_error = output[f"effect_se_{timepoint}"].to_numpy(dtype=float)
        margin = float(settings["negligible_effect_margin"])
        output[f"negligible_effect_supported_{timepoint}"] = (
            tested
            & ((effect - critical * standard_error) > -margin)
            & ((effect + critical * standard_error) < margin)
        )

    pairwise_p_values = []
    pairwise_locations: list[tuple[str, np.ndarray]] = []
    for first, second in combinations(timepoints, 2):
        name = contrast_name(first, second)
        values = output[f"effect_difference_p_value_{name}"].to_numpy(dtype=float)
        finite = np.flatnonzero(tested.to_numpy() & np.isfinite(values))
        output[f"effect_difference_fdr_edge_pair_{name}"] = np.nan
        pairwise_locations.append((name, finite))
        pairwise_p_values.append(values[finite])
    if sum(len(values) for values in pairwise_p_values):
        adjusted = bh_adjust(np.concatenate(pairwise_p_values))
        offset = 0
        for name, finite in pairwise_locations:
            width = len(finite)
            output.loc[finite, f"effect_difference_fdr_edge_pair_{name}"] = adjusted[
                offset : offset + width
            ]
            offset += width
    support_columns = [f"strict_total_support_{timepoint}" for timepoint in timepoints]
    output["strict_any_timepoint_support"] = output[support_columns].any(axis=1)
    output["strict_support_pattern"] = output.apply(
        support_pattern, axis=1, timepoints=timepoints, prefix="strict_total_support_"
    )
    output["strict_temporal_support_class"] = output.apply(
        classify_strict_temporal_shape,
        axis=1,
        timepoints=timepoints,
        settings=settings,
    )
    return output


def run_strict(
    config_path: Path,
    pseudobulk_dir: Path,
    edge_matrix_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["time_resolved_support"]
    timepoints = list(config["state_dependence"]["ordered_timepoints"])
    reference_guides = list(config["data"]["intervention_reference_guides"])

    z = np.load(pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_counts.npz")
    counts = z["counts"].astype(np.float64)
    gene_name = z["gene_name"].astype(str)
    gene_id = z["gene_id"].astype(str)
    keys = feature_keys(gene_name, gene_id)
    lookup = {key: index for index, key in enumerate(keys)}
    groups = pd.read_csv(
        pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_groups.csv"
    )
    library_total = counts.sum(axis=1, keepdims=True)
    log_cpm = np.log2(counts / np.maximum(library_total, 1) * 1_000_000 + 0.5)

    e0 = pd.read_csv(edge_matrix_path)
    e0["feature_index"] = e0.target.map(lookup)
    if e0.feature_index.isna().any():
        missing = e0.loc[e0.feature_index.isna(), "target"].tolist()
        raise ValueError(f"E0 targets absent from pseudobulk features: {missing[:5]}")

    blocks = []
    for _, candidates in e0.groupby("TF", observed=True):
        guides = _effective_guides(candidates)
        blocks.append(
            fit_tf_strict_total_support(
                candidates,
                log_cpm,
                groups,
                guides,
                reference_guides,
                timepoints,
                int(config["state_dependence"]["minimum_cells_per_pseudobulk"]),
            )
        )
    results = add_strict_support_calls(
        pd.concat(blocks, ignore_index=True), settings, timepoints
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir / "strict_time_resolved_total_support.csv.gz",
        index=False,
        compression="gzip",
    )
    results.loc[results.strict_any_timepoint_support].to_csv(
        output_dir / "strict_time_resolved_supported_edges.csv", index=False
    )
    summary = {
        "status": "strict_refit_completed",
        "estimand": "collection-time-specific total RNA intervention response",
        "original_day14_state_matched_E1_unchanged": True,
        "E0_edges_frozen": int(len(results)),
        "edge_time_hypotheses": int(
            results.strict_time_resolved_testable.sum() * len(timepoints)
        ),
        "testable_edges": int(results.strict_time_resolved_testable.sum()),
        "support_by_timepoint": {
            timepoint: int(results[f"strict_total_support_{timepoint}"].sum())
            for timepoint in timepoints
        },
        "any_timepoint_support": int(results.strict_any_timepoint_support.sum()),
        "support_patterns": {
            str(key): int(value)
            for key, value in results.strict_support_pattern.value_counts().items()
        },
        "temporal_support_classes": {
            str(key): int(value)
            for key, value in results.strict_temporal_support_class.value_counts().items()
        },
        "multiplicity": "BH across all finite frozen E0 edge-by-time tests",
        "limitations": [
            "total effects may include perturbation-induced lineage-composition shifts",
            "lack of support at a timepoint is not evidence of a zero effect",
            "resolved early/transient/late labels require globally adjusted contrasts and endpoint equivalence",
            "frozen day-14 guide efficacy is assumed to transport across collection time",
        ],
    }
    with (output_dir / "strict_time_resolved_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    return summary


def add_provisional_support_calls(
    table: pd.DataFrame,
    settings: dict,
    timepoints: list[str],
) -> pd.DataFrame:
    output = add_edge_time_fdr(table, timepoints)
    output["provisional_time_resolved_testable"] = (
        output.effective_guides_used >= int(settings["minimum_effective_guides"])
    )
    for timepoint in timepoints:
        output[f"provisional_total_support_{timepoint}"] = (
            output.provisional_time_resolved_testable
            & (
                output[f"effect_fdr_edge_time_{timepoint}"]
                < float(settings["global_fdr_max"])
            )
            & (
                output[f"effect_{timepoint}"].abs()
                >= float(settings["minimum_absolute_effect"])
            )
            & (
                output[f"consistent_guides_{timepoint}"]
                >= int(settings["minimum_consistent_guides"])
            )
        )
    support_columns = [
        f"provisional_total_support_{timepoint}" for timepoint in timepoints
    ]
    output["provisional_any_timepoint_support"] = output[support_columns].any(axis=1)
    output["provisional_support_pattern"] = output.apply(
        support_pattern,
        axis=1,
        timepoints=timepoints,
        prefix="provisional_total_support_",
    )
    output["provisional_temporal_support_class"] = output.apply(
        classify_temporal_support,
        axis=1,
        timepoints=timepoints,
        support_prefix="provisional_total_support_",
        testable_column="provisional_time_resolved_testable",
        interaction_fdr_max=float(settings["interaction_fdr_max"]),
    )
    return output


def run_provisional(
    config_path: Path,
    timepoint_edges_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["time_resolved_support"]
    timepoints = list(config["state_dependence"]["ordered_timepoints"])
    results = add_provisional_support_calls(
        pd.read_csv(timepoint_edges_path), settings, timepoints
    )
    preterminal = (
        results[
            [f"provisional_total_support_{timepoint}" for timepoint in timepoints[:-1]]
        ].any(axis=1)
        & ~results[f"provisional_total_support_{timepoints[-1]}"]
    )
    day14 = results[f"provisional_total_support_{timepoints[-1]}"]
    strict = results.E1_supported.fillna(False)
    summary = {
        "status": "provisional_aggregate_screen_only",
        "strict_time_resolved_E1_completed": False,
        "reason_strict_not_completed": (
            "timepoint-specific effective-guide and leave-one-guide-out fits cannot "
            "be reconstructed from the committed aggregate table"
        ),
        "original_day14_state_matched_E1_unchanged": True,
        "E0_edges_frozen": int(len(results)),
        "edge_time_hypotheses": int(len(results) * len(timepoints)),
        "provisional_support_by_timepoint": {
            timepoint: int(results[f"provisional_total_support_{timepoint}"].sum())
            for timepoint in timepoints
        },
        "provisional_any_timepoint_support": int(
            results.provisional_any_timepoint_support.sum()
        ),
        "preterminal_only_candidates": int(preterminal.sum()),
        "preterminal_only_with_detected_interaction": int(
            (preterminal & (results.interaction_fdr < float(settings["interaction_fdr_max"]))).sum()
        ),
        "provisional_any_timepoint_support_by_TF": {
            str(key): int(value)
            for key, value in results.loc[
                results.provisional_any_timepoint_support, "TF"
            ].value_counts().items()
        },
        "provisional_support_patterns": {
            str(key): int(value)
            for key, value in results.provisional_support_pattern.value_counts().items()
        },
        "provisional_temporal_support_classes": {
            str(key): int(value)
            for key, value in results.provisional_temporal_support_class.value_counts().items()
        },
        "day14_provisional_supported": int(day14.sum()),
        "day14_overlap_with_strict_E1": int((day14 & strict).sum()),
        "day14_provisional_only": int((day14 & ~strict).sum()),
        "day14_strict_only": int((~day14 & strict).sum()),
        "multiplicity": "BH across all 2476 frozen E0 edge-by-time tests",
        "limitations": [
            "all targeting guides contributed to the committed timepoint model",
            "timepoint-specific leave-one-effective-guide-out estimates are unavailable",
            "total effects may include perturbation-induced lineage-composition shifts",
            "lack of support at a timepoint is not evidence of a zero effect",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir / "provisional_time_resolved_total_support.csv.gz",
        index=False,
        compression="gzip",
    )
    with (output_dir / "provisional_time_resolved_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["provisional", "strict"], default="provisional")
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--timepoint-edges",
        type=Path,
        default=Path("reports/state_dependence/timepoint_interaction_edges.csv.gz"),
    )
    parser.add_argument(
        "--pseudobulk", type=Path, default=Path("data/processed/pseudobulk")
    )
    parser.add_argument(
        "--edge-matrix",
        type=Path,
        default=Path("reports/validation/E0_to_E1_edge_matrix.csv.gz"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/time_resolved_support")
    )
    args = parser.parse_args()
    if args.mode == "strict":
        summary = run_strict(
            args.config, args.pseudobulk, args.edge_matrix, args.output
        )
    else:
        summary = run_provisional(args.config, args.timepoint_edges, args.output)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
