import numpy as np
import pandas as pd

from edge_causality.discovery_selection_bias import (
    _allocate_integer,
    balanced_stratified_bootstrap_indices,
    classify_discovery_pattern,
)


def test_allocate_integer_preserves_total() -> None:
    observed = _allocate_integer(11, np.array([1, 1, 1]))
    assert observed.sum() == 11
    assert observed.max() - observed.min() <= 1


def test_power_matched_bootstrap_equalizes_libraries_and_preserves_n() -> None:
    metadata = pd.DataFrame(
        {
            "replicate": ["a"] * 8 + ["b"] * 4,
            "new_CellType": ["x"] * 6 + ["y"] * 2 + ["x"] * 2 + ["y"] * 2,
        }
    )
    rows = balanced_stratified_bootstrap_indices(
        metadata, np.random.default_rng(3), 10
    )
    sampled = metadata.iloc[rows]
    assert len(sampled) == 10
    assert sampled.replicate.value_counts().to_dict() == {"a": 5, "b": 5}
    assert sampled.loc[sampled.replicate.eq("a"), "new_CellType"].value_counts()["x"] == 4


def test_discovery_patterns_are_cautious() -> None:
    assert classify_discovery_pattern("1111") == "persistent_discovery"
    assert classify_discovery_pattern("0001") == "day14_only_discovery"
    assert classify_discovery_pattern("1000") == "early_only_discovery"
    assert classify_discovery_pattern("0010") == "day11_only_preterminal_discovery"
    assert classify_discovery_pattern("1011") == "cross_time_discovery"
