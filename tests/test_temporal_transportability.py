import numpy as np
import pandas as pd

from edge_causality.temporal_transportability import (
    build_transport_features,
    edge_block_from_matrix,
    group_bootstrap_metric_difference,
)


def test_edge_ranks_use_full_fixed_target_universe() -> None:
    pairs = pd.DataFrame(
        {
            "TF": ["A", "A", "A", "B", "B", "B"],
            "target": ["x", "y", "z", "x", "y", "z"],
        }
    )
    correlations = np.array(
        [
            [0.9, 0.1, -0.4],
            [0.2, -0.8, 0.3],
        ]
    )
    block = edge_block_from_matrix(
        pairs,
        correlations,
        ["A", "B"],
        {"x": 0, "y": 1, "z": 2},
    )
    a = block.loc[block.TF.eq("A")].set_index("target")
    b = block.loc[block.TF.eq("B")].set_index("target")
    assert a.loc["x", "rank_fraction"] == 1 / 3
    assert a.loc["y", "rank_fraction"] == 1
    assert b.loc["y", "rank_fraction"] == 1 / 3


def test_target_cluster_bootstrap_preserves_finite_metric_differences() -> None:
    frame = pd.DataFrame(
        {
            "target": ["a", "b", "c", "d", "e", "f"],
            "outcome": [1, 0, 1, 0, 1, 0],
            "baseline_prediction": [0.7, 0.6, 0.65, 0.55, 0.6, 0.5],
            "extended_prediction": [0.9, 0.2, 0.85, 0.15, 0.8, 0.1],
        }
    )
    differences = group_bootstrap_metric_difference(
        frame, "AUPRC", 50, np.random.default_rng(3)
    )
    assert len(differences) > 0
    assert np.isfinite(differences).all()
    assert np.median(differences) >= 0


def test_transport_features_keep_frozen_baseline_name() -> None:
    day14 = pd.DataFrame(
        {
            "TF": ["A"],
            "target": ["x"],
            "target_symbol": ["x"],
            "target_gene_id": ["id"],
            "signed_association": [0.4],
            "absolute_association": [0.4],
            "stable_edge": [True],
            "detection_fraction": [0.5],
            "mean_cpm": [10.0],
            "residual_rank_fraction": [0.01],
            "author_TF_sensitive": [True],
            "author_effect_025": [True],
            "author_supported_concordant": [True],
        }
    )
    block = pd.DataFrame(
        {
            "TF": ["A"],
            "target": ["x"],
            "signed_association": [0.3],
            "absolute_association": [0.3],
            "rank_fraction": [0.02],
            "loo_sign_fraction": [1.0],
            "loo_fisher_z_sd": [0.01],
            "detection_fraction": [0.5],
            "mean_cpm": [9.0],
        }
    )
    e1 = pd.DataFrame(
        {
            "TF": ["A"],
            "target": ["x"],
            "effective_guides_used": [2],
            "perturbation_log2fc": [-0.4],
            "perturbation_fdr": [0.01],
            "E1_supported": [True],
            "E1_direction_concordant": [True],
        }
    )
    trajectory = pd.DataFrame(
        {
            "TF": ["A"],
            "target": ["x"],
            "effect_day 7": [-0.5],
            "effect_day 9": [-0.4],
            "effect_day 11": [-0.3],
            "effect_day 14": [-0.2],
        }
    )
    settings = {
        "estimand": {"confirmatory_timepoints": ["day 7", "day 9", "day 11"]},
        "transportability": {
            "descriptive_transportable_rule": {
                "maximum_median_rank_fraction": 0.10,
                "minimum_leave_one_library_out_sign_fraction": 0.67,
            }
        },
    }
    observed = build_transport_features(
        day14,
        {"day 7": block, "day 9": block, "day 11": block},
        e1,
        trajectory,
        settings,
        {"day 7": 100, "day 9": 100, "day 11": 100},
    )
    assert observed.loc[0, "absolute_association"] == 0.4
    assert observed.loc[0, "prior_association_range"] == 0.0
    assert observed.loc[0, "prior_association_heterogeneity_p"] == 1.0
    assert (
        observed.loc[
            0, "observational_perturbation_direction_concordance_fraction"
        ]
        == 1.0
    )
    assert not observed.loc[0, "observational_perturbation_direction_conflict"]
