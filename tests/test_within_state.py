import numpy as np
import pandas as pd

from edge_causality.within_state import state_interaction_design


def test_state_design_produces_state_specific_contrasts() -> None:
    groups = pd.DataFrame(
        {
            "guide": ["A", "T", "A", "T"],
            "replicate": ["r1", "r1", "r1", "r1"],
            "cell_type": ["early", "early", "late", "late"],
        }
    )
    _, full, contrasts = state_interaction_design(groups, ["T"], ["early", "late"])
    assert full.shape == (4, 4)
    assert np.allclose(contrasts["early"], [0, 1, 0, 0])
    assert np.allclose(contrasts["late"], [0, 1, 0, 1])
