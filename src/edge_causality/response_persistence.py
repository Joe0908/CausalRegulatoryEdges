"""Separate absent time heterogeneity from absent perturbational response.

The collection-time model already stores four effect estimates and their
standard errors for every frozen E0 edge.  This module reconstructs the WLS
contrast covariance from the committed metadata, estimates a common response,
and applies equivalence-aware diagnostics to the interaction-null edges.

"No detected interaction" is deliberately retained as a heterogeneity result.
It is never used by itself as evidence for either stability or a null response.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
import yaml

from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import interaction_design


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def build_group_table(metadata: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct replicate-by-guide cell counts used by the WLS model."""
    required = {"replicate", "Timepoint", "perturbation_name", "target"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {sorted(missing)}")
    selected = metadata.dropna(subset=["replicate"]).copy()
    return (
        selected.groupby(
            ["replicate", "Timepoint", "perturbation_name", "target"],
            dropna=False,
            observed=True,
        )
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(
            columns={"Timepoint": "timepoint", "perturbation_name": "guide"}
        )
    )


def contrast_covariance_multiplier(
    groups: pd.DataFrame,
    target_guides: list[str],
    reference_guides: list[str],
    timepoints: list[str],
    minimum_cells: int,
) -> np.ndarray:
    """Return C(X'WX)^-1C' for the four timepoint effect contrasts."""
    selected = groups.loc[
        groups.guide.isin(target_guides + reference_guides)
        & (groups.n_cells >= minimum_cells)
    ].reset_index(drop=True)
    if selected.timepoint.nunique() < len(timepoints):
        raise ValueError("Insufficient timepoint coverage to reconstruct covariance")
    _, full, contrasts, _ = interaction_design(selected, target_guides, timepoints)
    weights = selected.n_cells.to_numpy(dtype=float)
    weights /= np.median(weights)
    weighted_design = full * np.sqrt(weights).reshape(-1, 1)
    inverse = np.linalg.pinv(weighted_design.T @ weighted_design)
    contrast_matrix = np.vstack([contrasts[timepoint] for timepoint in timepoints])
    return contrast_matrix @ inverse @ contrast_matrix.T


def estimate_common_response(
    row: pd.Series,
    covariance_multiplier: np.ndarray,
    timepoints: list[str],
    temporal_equivalence_margin: float,
    negligible_effect_margin: float,
    equivalence_alpha: float,
) -> dict[str, float | bool]:
    """Estimate a common effect and simultaneous all-timepoint diagnostics.

    The common effect is the GLS average of the four stored timepoint effects.
    Gene-specific residual variance is recovered from the stored standard
    errors and the reconstructed design multiplier.  Equivalence uses 90% CIs
    when ``equivalence_alpha`` is 0.05 (TOST convention).
    """
    effects = np.array(
        [row[f"effect_{timepoint}"] for timepoint in timepoints], dtype=float
    )
    standard_errors = np.array(
        [row[f"effect_se_{timepoint}"] for timepoint in timepoints], dtype=float
    )
    diagonal = np.diag(covariance_multiplier)
    sigma2_by_timepoint = standard_errors**2 / diagonal
    sigma2 = float(np.median(sigma2_by_timepoint))
    if not np.allclose(sigma2_by_timepoint, sigma2, rtol=1e-6, atol=1e-12):
        raise ValueError("Stored standard errors do not match reconstructed design")

    inverse_multiplier = np.linalg.pinv(covariance_multiplier)
    ones = np.ones(len(timepoints))
    precision = float(ones @ inverse_multiplier @ ones)
    gls_weights = inverse_multiplier @ ones / precision
    common_effect = float(gls_weights @ effects)
    common_se = float(np.sqrt(sigma2 / precision))
    degrees_freedom = int(row.model_df_denominator)
    if common_se > 0:
        common_p = float(
            2
            * student_t.sf(
                np.abs(common_effect / common_se), degrees_freedom
            )
        )
    else:
        common_p = 1.0

    critical = float(
        student_t.ppf(1 - equivalence_alpha, degrees_freedom)
    )
    temporal_equivalence = True
    for left in range(len(timepoints)):
        for right in range(left + 1, len(timepoints)):
            difference = effects[left] - effects[right]
            variance = sigma2 * (
                covariance_multiplier[left, left]
                + covariance_multiplier[right, right]
                - 2 * covariance_multiplier[left, right]
            )
            difference_se = float(np.sqrt(max(variance, 0)))
            lower = difference - critical * difference_se
            upper = difference + critical * difference_se
            temporal_equivalence &= (
                lower > -temporal_equivalence_margin
                and upper < temporal_equivalence_margin
            )

    negligible_at_all_timepoints = bool(
        np.all(effects - critical * standard_errors > -negligible_effect_margin)
        and np.all(effects + critical * standard_errors < negligible_effect_margin)
    )
    return {
        "common_effect": common_effect,
        "common_effect_se": common_se,
        "common_effect_p_value": common_p,
        "temporal_equivalence_supported": bool(temporal_equivalence),
        "uniformly_negligible_supported": negligible_at_all_timepoints,
    }


def add_response_classification(
    results: pd.DataFrame, settings: dict, timepoints: list[str]
) -> pd.DataFrame:
    """Add mutually exclusive response patterns and transparent component flags."""
    output = results.copy()
    effects = output[[f"effect_{timepoint}" for timepoint in timepoints]].to_numpy()
    absolute = np.abs(effects)
    guides = output[
        [f"consistent_guides_{timepoint}" for timepoint in timepoints]
    ].to_numpy()
    timepoint_fdr = output[
        [f"effect_fdr_{timepoint}" for timepoint in timepoints]
    ].to_numpy()

    interaction_null = output.interaction_fdr >= float(
        settings["interaction_fdr_max"]
    )
    same_direction = np.all(effects > 0, axis=1) | np.all(effects < 0, axis=1)
    guide_consistent = np.all(
        guides >= int(settings["minimum_consistent_guides"]), axis=1
    )
    guide_testable = output.effective_guides_used >= int(
        settings["minimum_effective_guides_for_response"]
    )
    moderate = float(settings["persistent_response_minimum_effect"])
    strong = float(settings["minimum_on_effect"])
    common_supported = output.common_effect_fdr < float(
        settings["common_effect_fdr_max"]
    )

    output["interaction_null"] = interaction_null
    output["all_timepoint_same_direction"] = same_direction
    output["all_timepoint_guide_consistent"] = guide_consistent
    output["guide_testable"] = guide_testable
    output["uniformly_small_point_estimates"] = np.all(
        absolute < float(settings["maximum_off_effect"]), axis=1
    )
    output["any_timepoint_response_supported"] = np.any(
        (timepoint_fdr < float(settings["common_effect_fdr_max"]))
        & (absolute >= strong),
        axis=1,
    )
    output["persistent_response_supported"] = (
        interaction_null
        & common_supported
        & (np.abs(output.common_effect) >= moderate)
        & same_direction
        & np.all(absolute >= moderate, axis=1)
        & guide_consistent
        & guide_testable
    )
    output["persistent_strong_response_supported"] = (
        output.persistent_response_supported
        & (np.abs(output.common_effect) >= strong)
    )
    output["temporally_stable_response_supported"] = (
        output.persistent_response_supported
        & output.temporal_equivalence_supported
    )
    output["uniformly_weak_supported"] = (
        interaction_null
        & output.uniformly_negligible_supported
        & guide_testable
    )

    patterns = np.full(len(output), "unresolved_no_interaction", dtype=object)
    patterns[output.uniformly_small_point_estimates.to_numpy()] = (
        "uniformly_small_estimates_only"
    )
    patterns[output.any_timepoint_response_supported.to_numpy()] = (
        "localized_response_no_heterogeneity"
    )
    patterns[output.uniformly_weak_supported.to_numpy()] = "uniformly_weak_supported"
    patterns[output.persistent_response_supported.to_numpy()] = (
        "persistent_response_supported"
    )
    patterns[(~interaction_null).to_numpy()] = "time_dependent"
    output["response_pattern"] = patterns
    return output


def add_unresolved_reason_audit(
    results: pd.DataFrame, settings: dict, timepoints: list[str]
) -> pd.DataFrame:
    """Reason-code the unresolved class without changing its frozen label."""
    output = results.copy()
    effects = output[[f"effect_{timepoint}" for timepoint in timepoints]].to_numpy()
    moderate = float(settings["persistent_response_minimum_effect"])
    common_supported = output.common_effect_fdr < float(
        settings["common_effect_fdr_max"]
    )
    unresolved = output.response_pattern.eq("unresolved_no_interaction")

    stable_small = (
        unresolved
        & output.guide_testable
        & common_supported
        & (np.abs(output.common_effect) < moderate)
        & output.all_timepoint_same_direction
        & output.all_timepoint_guide_consistent
        & output.temporal_equivalence_supported
    )
    near_threshold_persistent = (
        unresolved
        & output.guide_testable
        & common_supported
        & (np.abs(output.common_effect) >= moderate)
        & output.all_timepoint_same_direction
        & output.all_timepoint_guide_consistent
        & output.temporal_equivalence_supported
    )
    output["stable_small_response_supported"] = stable_small
    output["near_threshold_persistent"] = near_threshold_persistent

    reasons = np.full(len(output), "not_applicable", dtype=object)
    for position, (_, row) in enumerate(output.iterrows()):
        if row.response_pattern != "unresolved_no_interaction":
            continue
        row_effects = effects[position]
        if not row.guide_testable:
            reason = "guide_efficacy_limited"
        elif row.common_effect_fdr >= float(settings["common_effect_fdr_max"]):
            if not row.all_timepoint_same_direction:
                meaningful_polarity = (
                    row_effects.max() >= moderate
                    and row_effects.min() <= -moderate
                )
                if meaningful_polarity:
                    reason = "meaningful_polarity_conflict_interaction_unresolved"
                else:
                    reason = "near_zero_sign_crossing_common_unsupported"
            elif np.abs(row.common_effect) >= moderate:
                reason = "moderate_common_estimate_imprecise"
            else:
                reason = "small_same_direction_not_equivalently_weak"
        elif np.abs(row.common_effect) < moderate:
            reason = "small_common_response_supported"
        elif not row.all_timepoint_same_direction:
            reason = "common_response_direction_inconsistent"
        elif not row.all_timepoint_guide_consistent:
            reason = "common_response_guide_inconsistent"
        elif row.temporal_equivalence_supported:
            reason = "near_threshold_persistent_response"
        else:
            reason = "common_response_temporal_pattern_unresolved"
        reasons[position] = reason
    output["unresolved_reason"] = reasons

    refined = output.response_pattern.astype(object).to_numpy(copy=True)
    remaining = unresolved & ~stable_small & ~near_threshold_persistent
    refined[remaining.to_numpy()] = "reason_coded_unresolved"
    small_common = unresolved & output.unresolved_reason.eq(
        "small_common_response_supported"
    )
    refined[small_common.to_numpy()] = "small_common_response_supported"
    refined[stable_small.to_numpy()] = "temporally_stable_small_response_supported"
    refined[near_threshold_persistent.to_numpy()] = "near_threshold_persistent"
    output["refined_response_pattern"] = refined
    return output


def run(
    config_path: Path,
    interaction_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["state_dependence"]
    timepoints = list(settings["ordered_timepoints"])
    reference_guides = list(config["data"]["intervention_reference_guides"])
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0)
    groups = build_group_table(metadata)
    interactions = pd.read_csv(interaction_path)

    blocks = []
    for tf, block in interactions.groupby("TF", observed=True, sort=False):
        guide_sets = block.targeting_guides.drop_duplicates()
        if len(guide_sets) != 1:
            raise ValueError(f"Inconsistent targeting guide sets for {tf}")
        target_guides = str(guide_sets.iloc[0]).split(";")
        multiplier = contrast_covariance_multiplier(
            groups,
            target_guides,
            reference_guides,
            timepoints,
            int(settings["minimum_cells_per_pseudobulk"]),
        )
        estimates = block.apply(
            estimate_common_response,
            axis=1,
            result_type="expand",
            covariance_multiplier=multiplier,
            timepoints=timepoints,
            temporal_equivalence_margin=float(settings["minimum_effect_range"]),
            negligible_effect_margin=float(settings["maximum_off_effect"]),
            equivalence_alpha=float(settings["equivalence_alpha"]),
        ).reset_index(drop=True)
        blocks.append(pd.concat([block.reset_index(drop=True), estimates], axis=1))

    results = pd.concat(blocks, ignore_index=True)
    results["common_effect_fdr"] = bh_adjust(
        results.common_effect_p_value.to_numpy()
    )
    results = add_response_classification(results, settings, timepoints)
    results = add_unresolved_reason_audit(results, settings, timepoints)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir / "noninteraction_response_patterns.csv.gz",
        index=False,
        compression="gzip",
    )
    audit_columns = [
        "TF",
        "target",
        "response_pattern",
        "refined_response_pattern",
        "unresolved_reason",
        "stable_small_response_supported",
        "near_threshold_persistent",
        "common_effect",
        "common_effect_se",
        "common_effect_fdr",
        "interaction_fdr",
        "guide_testable",
        "effective_guides_used",
        "all_timepoint_same_direction",
        "all_timepoint_guide_consistent",
        "temporal_equivalence_supported",
        "uniformly_negligible_supported",
        *[f"effect_{timepoint}" for timepoint in timepoints],
        *[f"effect_fdr_{timepoint}" for timepoint in timepoints],
        *[f"consistent_guides_{timepoint}" for timepoint in timepoints],
        "E1_supported",
    ]
    unresolved_audit = results.loc[
        results.response_pattern.eq("unresolved_no_interaction"), audit_columns
    ]
    unresolved_audit.to_csv(
        output_dir / "unresolved_reason_audit.csv", index=False
    )

    interaction_null = results.interaction_null
    pattern_counts = {
        str(key): int(value)
        for key, value in results.response_pattern.value_counts().items()
    }
    unresolved_reason_counts = {
        str(key): int(value)
        for key, value in unresolved_audit.unresolved_reason.value_counts().items()
    }
    refined_pattern_counts = {
        str(key): int(value)
        for key, value in results.refined_response_pattern.value_counts().items()
    }
    summary = {
        "E0_edges_tested": int(len(results)),
        "time_dependent": int((~interaction_null).sum()),
        "no_detected_interaction": int(interaction_null.sum()),
        "response_patterns": pattern_counts,
        "refined_response_patterns": refined_pattern_counts,
        "unresolved_reason_counts": unresolved_reason_counts,
        "no_interaction_uniformly_small_point_estimates": int(
            (interaction_null & results.uniformly_small_point_estimates).sum()
        ),
        "no_interaction_any_timepoint_response_supported": int(
            (interaction_null & results.any_timepoint_response_supported).sum()
        ),
        "no_interaction_common_effect_fdr_lt_0_05": int(
            (interaction_null & (results.common_effect_fdr < 0.05)).sum()
        ),
        "persistent_response_supported": int(
            results.persistent_response_supported.sum()
        ),
        "persistent_strong_response_supported": int(
            results.persistent_strong_response_supported.sum()
        ),
        "temporally_stable_response_supported": int(
            results.temporally_stable_response_supported.sum()
        ),
        "uniformly_negligible_supported_all_TFs": int(
            (interaction_null & results.uniformly_negligible_supported).sum()
        ),
        "uniformly_weak_supported_guide_testable": int(
            results.uniformly_weak_supported.sum()
        ),
        "stable_small_response_supported": int(
            results.stable_small_response_supported.sum()
        ),
        "near_threshold_persistent": int(
            results.near_threshold_persistent.sum()
        ),
        "reason_coded_unresolved_remaining": int(
            results.refined_response_pattern.eq("reason_coded_unresolved").sum()
        ),
        "thresholds": {
            "common_effect_fdr_max": float(settings["common_effect_fdr_max"]),
            "persistent_response_minimum_effect": float(
                settings["persistent_response_minimum_effect"]
            ),
            "strong_response_minimum_effect": float(settings["minimum_on_effect"]),
            "negligible_effect_margin": float(settings["maximum_off_effect"]),
            "temporal_equivalence_margin": float(settings["minimum_effect_range"]),
            "equivalence_alpha": float(settings["equivalence_alpha"]),
        },
        "limitations": [
            "interaction non-significance is not temporal equivalence",
            "uniformly small point estimates are not evidence of negligible effects",
            "weak-response claims require at least two effective guides",
            "reason codes identify the evidence blocker, not a causal mechanism",
            "stable small responses are not proven negligible at every timepoint",
            "replicate labels are library/batch strata, not independent donors",
        ],
    }
    with (output_dir / "noninteraction_response_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--interactions",
        type=Path,
        default=Path("reports/state_dependence/timepoint_interaction_edges.csv.gz"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/state_dependence")
    )
    args = parser.parse_args()
    summary = run(args.config, args.interactions, args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
