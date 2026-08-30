import numpy as np
import pandas as pd

from edge_causality.time_resolved_support import (
    add_edge_time_fdr,
    add_provisional_support_calls,
    add_strict_support_calls,
    classify_temporal_support,
    fit_tf_strict_total_support,
)


TIMEPOINTS = ["day 7", "day 9", "day 11", "day 14"]


def test_edge_time_fdr_pools_all_timepoints() -> None:
    table = pd.DataFrame(
        {
            f"effect_p_value_{timepoint}": [0.001, 0.9]
            for timepoint in TIMEPOINTS
        }
    )
    observed = add_edge_time_fdr(table, TIMEPOINTS)
    for timepoint in TIMEPOINTS:
        assert np.isclose(observed.loc[0, f"effect_fdr_edge_time_{timepoint}"], 0.002)
        assert np.isclose(observed.loc[1, f"effect_fdr_edge_time_{timepoint}"], 0.9)


def test_provisional_calls_require_effect_size_guides_and_global_fdr() -> None:
    rows = []
    for target, effective_guides, effects, p_values, consistent in [
        ("supported", 2, [0.4, 0.0, 0.0, 0.0], [1e-6, 1.0, 1.0, 1.0], 2),
        ("guide_limited", 1, [0.4, 0.0, 0.0, 0.0], [1e-6, 1.0, 1.0, 1.0], 2),
        ("discordant", 2, [0.4, 0.0, 0.0, 0.0], [1e-6, 1.0, 1.0, 1.0], 1),
    ]:
        row = {
            "TF": "TF",
            "target": target,
            "effective_guides_used": effective_guides,
            "interaction_fdr": 0.01,
        }
        for timepoint, effect, p_value in zip(TIMEPOINTS, effects, p_values):
            row[f"effect_{timepoint}"] = effect
            row[f"effect_p_value_{timepoint}"] = p_value
            row[f"consistent_guides_{timepoint}"] = consistent
        rows.append(row)
    settings = {
        "minimum_effective_guides": 2,
        "minimum_consistent_guides": 2,
        "minimum_absolute_effect": 0.25,
        "global_fdr_max": 0.05,
        "interaction_fdr_max": 0.05,
    }
    observed = add_provisional_support_calls(
        pd.DataFrame(rows), settings, TIMEPOINTS
    ).set_index("target")
    assert observed.loc["supported", "provisional_total_support_day 7"]
    assert observed.loc["supported", "provisional_support_pattern"] == "1000"
    assert not observed.loc["guide_limited", "provisional_any_timepoint_support"]
    assert not observed.loc["discordant", "provisional_any_timepoint_support"]


def test_temporal_labels_do_not_call_non_supported_days_null() -> None:
    row = {"testable": True, "interaction_fdr": 0.8}
    for timepoint, supported in zip(TIMEPOINTS, [True, False, False, False]):
        row[f"support_{timepoint}"] = supported
    assert classify_temporal_support(
        pd.Series(row), TIMEPOINTS, "support_", "testable", 0.05
    ) == "localized_support_no_detected_heterogeneity"

    row["interaction_fdr"] = 0.01
    assert classify_temporal_support(
        pd.Series(row), TIMEPOINTS, "support_", "testable", 0.05
    ) == "early_window_candidate_with_heterogeneity"


def test_strict_refit_uses_effective_guides_and_loo() -> None:
    rows = []
    values = []
    guide_offsets = {
        "AAVS1_1": -0.02,
        "AAVS1_2": 0.02,
        "TF_1": 0.48,
        "TF_2": 0.50,
        "TF_3": 0.52,
        "TF_bad": -2.0,
    }
    for day_index, timepoint in enumerate(TIMEPOINTS):
        for replicate_index in range(2):
            replicate = f"rep{2 * day_index + replicate_index + 1}"
            baseline = 2.0 + 0.1 * day_index + 0.01 * replicate_index
            for guide, offset in guide_offsets.items():
                rows.append(
                    {
                        "guide": guide,
                        "replicate": replicate,
                        "timepoint": timepoint,
                        "n_cells": 50,
                    }
                )
                values.append([baseline + offset])
    candidates = pd.DataFrame(
        {
            "TF": ["TF"],
            "target": ["gene"],
            "feature_index": [0],
            "signed_association": [-1.0],
        }
    )
    fitted = fit_tf_strict_total_support(
        candidates,
        np.asarray(values),
        pd.DataFrame(rows),
        ["TF_1", "TF_2", "TF_3"],
        ["AAVS1_1", "AAVS1_2"],
        TIMEPOINTS,
        15,
    )
    settings = {
        "global_fdr_max": 0.05,
        "minimum_absolute_effect": 0.25,
        "interaction_fdr_max": 0.05,
        "contrast_fdr_max": 0.05,
        "negligible_effect_margin": 0.125,
        "equivalence_alpha": 0.05,
    }
    observed = add_strict_support_calls(fitted, settings, TIMEPOINTS).iloc[0]
    assert observed.strict_time_resolved_testable
    assert observed.strict_support_pattern == "1111"
    assert observed.strict_temporal_support_class == "persistent_support_detected"
    for timepoint in TIMEPOINTS:
        assert np.isclose(observed[f"effect_{timepoint}"], 0.5)
        assert observed[f"consistent_effective_guides_{timepoint}"] == 3
        assert observed[
            f"leave_one_effective_guide_out_direction_consistent_{timepoint}"
        ]
