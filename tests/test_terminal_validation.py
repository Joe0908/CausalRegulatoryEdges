import numpy as np
import pandas as pd

from edge_causality.terminal_validation import (
    activation_stage,
    bootstrap_lead_fraction,
    map_ludwig_peaks,
    parse_library_name,
)


def test_parse_library_name_recovers_design_fields():
    result = parse_library_name("Donor3_P7_Rep4_ATAC")
    assert result["donor"] == "Donor3"
    assert result["population"] == "P7"
    assert result["replicate"] == "Rep4"


def test_activation_stage_uses_half_dynamic_range():
    onset, dynamic = activation_stage(np.array([0, 0.2, 0.8, 2.0]), 0.5, 1.0)
    assert onset == 3
    assert dynamic == 2.0


def test_activation_stage_rejects_small_or_baseline_maximum():
    assert np.isnan(activation_stage(np.array([1.0, 1.1, 1.0]), 0.5, 1.0)[0])
    assert np.isnan(activation_stage(np.array([3.0, 2.0, 1.0]), 0.5, 1.0)[0])


def test_map_ludwig_peaks_requires_interval_overlap():
    candidates = pd.DataFrame(
        {
            "peak_id": ["a"],
            "hg19_chromosome": ["chr1"],
            "hg19_start": [100],
            "hg19_end": [200],
            "liftover_pass": [True],
        }
    )
    peaks = pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "start": [150, 200],
            "end": [175, 225],
        }
    )
    result = map_ludwig_peaks(candidates, peaks)
    assert result.ludwig_peak_index.tolist() == [0]


def test_bootstrap_lead_fraction_is_high_for_clear_lead():
    ordered = ["P1", "P2", "P3", "P4"]
    populations = np.repeat(ordered, 3)
    atac = np.repeat([0.0, 2.0, 3.0, 3.0], 3)
    rna = np.repeat([0.0, 0.0, 0.2, 3.0], 3)
    fraction = bootstrap_lead_fraction(
        atac,
        populations,
        rna,
        populations,
        ordered,
        0.5,
        1.0,
        1.0,
        1,
        100,
        np.random.default_rng(3),
    )
    assert fraction == 1.0
