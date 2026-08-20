import numpy as np
import pandas as pd

from edge_causality.evaluate_author_truth import matched_null_rates


def test_matched_null_returns_requested_iterations() -> None:
    edges = pd.DataFrame(
        {
            "TF": ["A"] * 4,
            "target": list("wxyz"),
            "expression_bin": [0, 0, 0, 0],
            "detection_bin": [0, 0, 0, 0],
            "truth": [1, 0, 1, 0],
        }
    )
    observed = matched_null_rates(
        edges,
        np.array([True, False, False, False]),
        "truth",
        25,
        np.random.default_rng(1),
    )
    assert len(observed) == 25
    assert np.all((observed >= 0) & (observed <= 1))
