import pandas as pd

from edge_causality.guide_signature_sensitivity import effective_guide_flags


def test_effective_guide_rule_requires_all_frozen_components():
    table = pd.DataFrame({
        "crossfit_median_score": [1.0, 0.1, 1.0, 1.0],
        "bootstrap_null_q95": [0.2, 0.2, 0.2, 0.2],
        "bootstrap_fdr": [0.01, 0.01, 0.20, 0.01],
        "direction_spearman": [0.5, 0.5, 0.5, -0.1],
    })
    flags = effective_guide_flags(
        table, {"guide_fdr_max": 0.05, "minimum_direction_correlation": 0.0}
    )
    assert flags.tolist() == [True, False, False, False]


def test_effective_guide_rule_accepts_fdr_boundary():
    table = pd.DataFrame({
        "crossfit_median_score": [0.3],
        "bootstrap_null_q95": [0.2],
        "bootstrap_fdr": [0.05],
        "direction_spearman": [0.01],
    })
    assert effective_guide_flags(
        table, {"guide_fdr_max": 0.05, "minimum_direction_correlation": 0.0}
    ).iloc[0]
