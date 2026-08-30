import numpy as np
import pandas as pd

from edge_causality.external_validation import (
    commitment_metrics,
    frozen_peak_sets,
    interval_distance,
    interval_overlap,
    parse_interval,
    standardized_state_effect,
    summarize_results,
)


def test_parse_interval_accepts_both_project_and_external_formats():
    assert parse_interval("chr1:10-20") == ("chr1", 10, 20)
    assert parse_interval("chrX-30-45") == ("chrX", 30, 45)


def test_interval_geometry_uses_half_open_coordinates():
    assert interval_overlap(10, 20, 15, 25) == 5
    assert interval_overlap(10, 20, 20, 25) == 0
    assert interval_distance(10, 20, 22, 25) == 2


def test_frozen_peak_sets_separates_E2_from_linked_comparisons():
    evidence = pd.DataFrame(
        {
            "TF": ["GATA1", "GATA1", "NFE2"],
            "target": ["A", "B", "C"],
            "peak_id": ["p1", "p2", "p3"],
            "chromosome": ["chr1"] * 3,
            "start": [1, 2, 3],
            "end": [2, 3, 4],
            "distance_to_tss": [0, 0, 0],
            "candidate_role": ["x", "x", "x"],
            "E2_peak": [True, False, True],
            "link_pass": [True, True, True],
        }
    )
    result = frozen_peak_sets(evidence, "GATA1")
    assert result.set_index("peak_id").external_set.to_dict() == {
        "p1": "targeted_E2",
        "p2": "linked_non_E2",
    }


def test_standardized_state_effect_is_positive_for_establishment():
    values = np.array([0.0, 0.1, 1.0, 1.1])
    states = np.array(["MPP", "MPP", "Comm-Prog", "Comm-Prog"])
    assert standardized_state_effect(values, states, "MPP", "Comm-Prog") > 1


def test_commitment_metrics_identifies_late_induction():
    metrics = commitment_metrics(1.0, 2.0, 9.0)
    assert metrics["rna_commitment_fraction"] == 0.125
    assert metrics["rna_late_log2_fold_change"] > 1


def test_summary_requires_both_atac_establishment_and_rna_delay():
    peaks = pd.DataFrame(
        {
            "TF": ["GATA1", "GATA1"],
            "target": ["A", "B"],
            "external_set": ["targeted_E2", "linked_non_E2"],
            "peak_id": ["p1", "p2"],
            "mapped_external_peaks": [1, 1],
            "atac_establishment_pass": [True, False],
        }
    )
    rna = pd.DataFrame(
        {"target": ["A", "B"], "rna_delay_pass": [True, True]}
    )
    edges, summary = summarize_results(peaks, rna)
    assert edges.set_index("target").external_edge_pass.to_dict() == {
        "A": True,
        "B": False,
    }
    assert summary["external_edges_passed"] == 1
