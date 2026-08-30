import numpy as np

from edge_causality.build_targeted_multiome import (
    local_peak_rows,
    merge_intervals,
    parse_interval,
)


def test_parse_interval():
    assert parse_interval("chr22:10-25") == ("chr22", 10, 25)


def test_merge_intervals_unites_overlaps_but_not_adjacent():
    intervals = [
        ("chr1", 30, 40),
        ("chr1", 10, 20),
        ("chr1", 18, 25),
        ("chr1", 25, 27),
        ("chr2", 1, 3),
    ]
    assert merge_intervals(intervals) == [
        ("chr1", 10, 25),
        ("chr1", 25, 27),
        ("chr1", 30, 40),
        ("chr2", 1, 3),
    ]


def test_local_peak_rows_uses_half_open_overlap():
    types = np.array(["Gene Expression", "Peaks", "Peaks", "Peaks"])
    intervals = np.array(["chr1:1-2", "chr1:10-20", "chr1:20-30", "chr2:10-20"])
    assert local_peak_rows(types, intervals, "chr1", 15, 20) == [
        (1, "chr1", 10, 20)
    ]
