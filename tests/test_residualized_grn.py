import numpy as np
import pandas as pd

from edge_causality.residualized_grn import residualize, stratified_bootstrap_indices


def test_residualize_removes_linear_design_signal() -> None:
    design = np.column_stack([np.ones(10), np.arange(10)])
    values = (3 + 2 * np.arange(10)).reshape(-1, 1).astype(float)
    observed = residualize(values, design)
    assert np.max(np.abs(observed)) < 1e-10


def test_stratified_bootstrap_preserves_stratum_sizes() -> None:
    metadata = pd.DataFrame(
        {
            "replicate": ["a", "a", "b", "b", "b"],
            "new_CellType": ["x", "x", "y", "y", "y"],
        }
    )
    sampled = stratified_bootstrap_indices(metadata, np.random.default_rng(1))
    assert len(sampled) == len(metadata)
    assert np.sum(sampled < 2) == 2
    assert np.sum(sampled >= 2) == 3
