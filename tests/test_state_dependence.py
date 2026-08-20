import numpy as np
import pandas as pd

from edge_causality.state_dependence import (
    classify_edge,
    fit_interaction_model,
    interaction_design,
)


TIMEPOINTS = ["day 7", "day 9", "day 11", "day 14"]


def test_interaction_model_recovers_context_effects() -> None:
    rows = []
    response = []
    true_effects = {"day 7": 0.0, "day 9": 1.0, "day 11": 2.0, "day 14": 3.0}
    for index, timepoint in enumerate(TIMEPOINTS):
        for replicate in [f"rep{2 * index + 1}", f"rep{2 * index + 2}"]:
            for guide, condition in [("AAVS1_1", 0), ("TF_1", 1)]:
                rows.append(
                    {"guide": guide, "replicate": replicate, "timepoint": timepoint}
                )
                response.append(10 + index + condition * true_effects[timepoint])
    groups = pd.DataFrame(rows)
    reduced, full, contrasts, _ = interaction_design(groups, ["TF_1"], TIMEPOINTS)
    fitted = fit_interaction_model(
        np.asarray(response)[:, None],
        reduced,
        full,
        np.ones(len(response)),
        contrasts,
    )
    for timepoint, expected in true_effects.items():
        assert np.allclose(fitted[f"effect_{timepoint}"], expected)


def test_classify_gated_edge() -> None:
    row = pd.Series(
        {
            "interaction_fdr": 0.001,
            "effect_day 7": 0.02,
            "effect_day 9": 0.10,
            "effect_day 11": 0.35,
            "effect_day 14": 0.50,
            "consistent_guides_day 7": 1,
            "consistent_guides_day 9": 2,
            "consistent_guides_day 11": 3,
            "consistent_guides_day 14": 3,
        }
    )
    settings = {
        "interaction_fdr_max": 0.05,
        "minimum_on_effect": 0.25,
        "maximum_off_effect": 0.125,
        "minimum_effect_range": 0.25,
        "minimum_consistent_guides": 2,
    }
    assert classify_edge(row, settings, TIMEPOINTS) == "gated"
