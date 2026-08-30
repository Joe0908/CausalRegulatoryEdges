import pandas as pd

from edge_causality.trajectory_shift import numeric_timepoint


def test_numeric_timepoint_parses_labels() -> None:
    observed = numeric_timepoint(pd.Series(["day 7", "day 9", "day 11", "day 14"]))
    assert observed.tolist() == [7.0, 9.0, 11.0, 14.0]
