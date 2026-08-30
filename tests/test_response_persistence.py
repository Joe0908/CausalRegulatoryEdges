import numpy as np
import pandas as pd

from edge_causality.response_persistence import (
    add_response_classification,
    add_unresolved_reason_audit,
    contrast_covariance_multiplier,
    estimate_common_response,
)


TIMEPOINTS = ["day 7", "day 9", "day 11", "day 14"]


def synthetic_groups() -> pd.DataFrame:
    rows = []
    for index, timepoint in enumerate(TIMEPOINTS):
        for replicate in [f"rep{2 * index + 1}", f"rep{2 * index + 2}"]:
            for guide in ["AAVS1_1", "AAVS1_2", "TF_1", "TF_2"]:
                rows.append(
                    {
                        "replicate": replicate,
                        "timepoint": timepoint,
                        "guide": guide,
                        "n_cells": 50,
                    }
                )
    return pd.DataFrame(rows)


def response_row(
    effects: list[float], multiplier: np.ndarray, residual_variance: float = 0.01
) -> pd.Series:
    row: dict[str, float] = {"model_df_denominator": 20}
    for index, timepoint in enumerate(TIMEPOINTS):
        row[f"effect_{timepoint}"] = effects[index]
        row[f"effect_se_{timepoint}"] = np.sqrt(
            residual_variance * multiplier[index, index]
        )
    return pd.Series(row)


def test_common_response_recovers_persistent_effect_and_equivalence() -> None:
    multiplier = contrast_covariance_multiplier(
        synthetic_groups(),
        ["TF_1", "TF_2"],
        ["AAVS1_1", "AAVS1_2"],
        TIMEPOINTS,
        15,
    )
    estimate = estimate_common_response(
        response_row([0.2, 0.2, 0.2, 0.2], multiplier),
        multiplier,
        TIMEPOINTS,
        temporal_equivalence_margin=0.25,
        negligible_effect_margin=0.125,
        equivalence_alpha=0.05,
    )
    assert np.isclose(estimate["common_effect"], 0.2)
    assert estimate["common_effect_p_value"] < 0.05
    assert estimate["temporal_equivalence_supported"]
    assert not estimate["uniformly_negligible_supported"]


def test_zero_effect_can_support_uniform_negligibility() -> None:
    multiplier = contrast_covariance_multiplier(
        synthetic_groups(),
        ["TF_1", "TF_2"],
        ["AAVS1_1", "AAVS1_2"],
        TIMEPOINTS,
        15,
    )
    estimate = estimate_common_response(
        response_row([0.0, 0.0, 0.0, 0.0], multiplier, residual_variance=0.0001),
        multiplier,
        TIMEPOINTS,
        temporal_equivalence_margin=0.25,
        negligible_effect_margin=0.125,
        equivalence_alpha=0.05,
    )
    assert estimate["uniformly_negligible_supported"]


def test_classification_does_not_equate_interaction_null_with_no_effect() -> None:
    rows = []
    for target, effect, common_fdr, negligible in [
        ("persistent", 0.2, 0.001, False),
        ("weak", 0.0, 1.0, True),
    ]:
        row = {
            "TF": "TF",
            "target": target,
            "interaction_fdr": 0.8,
            "common_effect": effect,
            "common_effect_fdr": common_fdr,
            "temporal_equivalence_supported": target == "persistent",
            "uniformly_negligible_supported": negligible,
            "effective_guides_used": 2,
        }
        for timepoint in TIMEPOINTS:
            row[f"effect_{timepoint}"] = effect
            row[f"effect_fdr_{timepoint}"] = common_fdr
            row[f"consistent_guides_{timepoint}"] = 2
        rows.append(row)
    settings = {
        "interaction_fdr_max": 0.05,
        "minimum_consistent_guides": 2,
        "minimum_effective_guides_for_response": 2,
        "persistent_response_minimum_effect": 0.125,
        "minimum_on_effect": 0.25,
        "maximum_off_effect": 0.125,
        "common_effect_fdr_max": 0.05,
    }
    classified = add_response_classification(
        pd.DataFrame(rows), settings, TIMEPOINTS
    ).set_index("target")
    assert classified.loc["persistent", "response_pattern"] == (
        "persistent_response_supported"
    )
    assert classified.loc["weak", "response_pattern"] == "uniformly_weak_supported"


def test_unresolved_audit_preserves_original_label_and_adds_reasons() -> None:
    rows = []
    specifications = [
        ("guide_limited", [0.2, 0.2, 0.2, 0.2], 0.001, False, False, 1),
        ("stable_small", [0.08, 0.09, 0.07, 0.08], 0.001, True, True, 2),
        ("near_threshold", [0.13, 0.14, 0.11, 0.15], 0.001, True, True, 2),
        ("polarity", [0.2, -0.2, 0.1, -0.1], 0.8, False, False, 2),
    ]
    for target, effects, common_fdr, same, equivalent, effective_guides in specifications:
        row = {
            "TF": "TF",
            "target": target,
            "response_pattern": "unresolved_no_interaction",
            "guide_testable": effective_guides >= 2,
            "effective_guides_used": effective_guides,
            "common_effect": float(np.mean(effects)),
            "common_effect_fdr": common_fdr,
            "all_timepoint_same_direction": same,
            "all_timepoint_guide_consistent": target != "guide_inconsistent",
            "temporal_equivalence_supported": equivalent,
        }
        for timepoint, effect in zip(TIMEPOINTS, effects):
            row[f"effect_{timepoint}"] = effect
        rows.append(row)
    settings = {
        "persistent_response_minimum_effect": 0.125,
        "common_effect_fdr_max": 0.05,
    }
    audited = add_unresolved_reason_audit(
        pd.DataFrame(rows), settings, TIMEPOINTS
    ).set_index("target")
    assert audited.response_pattern.eq("unresolved_no_interaction").all()
    assert audited.loc["guide_limited", "unresolved_reason"] == (
        "guide_efficacy_limited"
    )
    assert audited.loc["stable_small", "refined_response_pattern"] == (
        "temporally_stable_small_response_supported"
    )
    assert audited.loc["near_threshold", "refined_response_pattern"] == (
        "near_threshold_persistent"
    )
    assert audited.loc["polarity", "unresolved_reason"] == (
        "meaningful_polarity_conflict_interaction_unresolved"
    )
