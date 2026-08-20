import numpy as np
import pandas as pd

from edge_causality.within_state_chromatin import fit_wls_condition


def test_fit_wls_condition_recovers_known_effect():
    groups = pd.DataFrame(
        {
            "guide": ["ctrl", "target"] * 4,
            "replicate": ["r1", "r1", "r2", "r2", "r3", "r3", "r4", "r4"],
            "cell_state": ["a", "a", "b", "b", "a", "a", "b", "b"],
            "n_cells": [20] * 8,
        }
    )
    response = np.array([[0, 1, 0, 1, 0, 1, 0, 1]], dtype=float).T
    effect, p_value = fit_wls_condition(response, groups, ["target"], True)
    assert np.isclose(effect[0], 1.0)
    assert p_value[0] < 0.05
