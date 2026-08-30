import numpy as np

from edge_causality.e2_library_robustness import _random_effects


def test_random_effects_preserves_concordant_link_direction():
    result = _random_effects(np.array([0.20, 0.22, 0.18]), np.array([100, 120, 90]))
    assert result["meta_link_correlation"] > 0
    assert result["meta_link_p_value"] < 0.05
    assert result["meta_tau2"] >= 0


def test_random_effects_returns_missing_for_unevaluable_strata():
    result = _random_effects(np.array([np.nan, 0.2]), np.array([100, 3]))
    assert np.isnan(result["meta_link_correlation"])
