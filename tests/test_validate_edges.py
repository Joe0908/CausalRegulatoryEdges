import numpy as np

from edge_causality.validate_edges import weighted_ols


def test_weighted_ols_recovers_condition_effect() -> None:
    condition = np.array([0, 0, 1, 1], dtype=float)
    design = np.column_stack([np.ones(4), condition])
    response = np.column_stack([2 + 3 * condition, 5 - condition])
    beta, p_value = weighted_ols(response, design, np.ones(4), coefficient=1)
    assert np.allclose(beta, [3, -1])
    assert np.all(p_value <= 1)
