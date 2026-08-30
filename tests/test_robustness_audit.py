import pandas as pd

from edge_causality.robustness_audit import call_e1, call_e2, e0_stability_grid


def test_e0_grid_reproduces_frozen_rule():
    edges = pd.DataFrame(
        {
            "TF": ["A", "A"],
            "target": ["x", "y"],
            "absolute_association": [0.2, 0.1],
            "bootstrap_top1_frequency": [0.8, 0.1],
            "bootstrap_top5_frequency": [0.9, 0.6],
            "bootstrap_top10_frequency": [1.0, 0.8],
        }
    )
    grid = {"selection_percent": [5], "minimum_bootstrap_frequency": [0.7]}
    result = e0_stability_grid(edges, grid)
    assert result.loc[0, "selected_edges"] == 1
    assert result.loc[0, "jaccard_with_frozen"] == 1


def test_e1_guide_robustness_is_enforced():
    edges = pd.DataFrame(
        {
            "effective_guides_used": [2, 2],
            "perturbation_fdr": [0.01, 0.01],
            "perturbation_log2fc": [0.4, 0.4],
            "guide_direction_consistent": [True, False],
            "leave_one_guide_out_direction_consistent": [True, True],
        }
    )
    assert call_e1(edges).tolist() == [True, False]


def test_e2_reclassification_uses_all_layers():
    evidence = pd.DataFrame(
        {
            "link_fdr": [0.01, 0.01],
            "link_correlation": [0.1, 0.1],
            "link_bootstrap_sign_fraction": [0.9, 0.9],
            "libraries_present": [5, 5],
            "strongest_effect_fdr": [0.01, 0.01],
            "strongest_effect": [-0.5, -0.5],
            "strongest_consistent_guides": [3, 3],
            "rna_effect_at_strongest_atac_timepoint": [-0.4, -0.4],
            "motif_best_relative_score": [0.9, 0.7],
        }
    )
    assert call_e2(evidence).tolist() == [True, False]
