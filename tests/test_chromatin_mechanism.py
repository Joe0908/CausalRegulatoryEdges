import numpy as np
import pandas as pd
from scipy import sparse

from edge_causality.chromatin_mechanism import (
    add_baseline_gating_diagnostics,
    bootstrap_sign_fraction,
    covariate_design,
    residualize,
)


def test_residualize_removes_linear_covariate():
    x = np.arange(10, dtype=float)
    design = np.column_stack([np.ones(10), x])
    residual = residualize(3 + 2 * x, design)
    assert np.max(np.abs(residual)) < 1e-10


def test_covariate_design_has_intercept_and_correct_rows():
    metadata = pd.DataFrame(
        {
            "nCount_RNA": [10, 20, 30],
            "nCount_atac": [100, 200, 300],
            "replicate": ["rep1", "rep1", "rep2"],
            "new_CellType": ["A", "B", "B"],
        }
    )
    design = covariate_design(metadata)
    assert design.shape[0] == 3
    assert np.all(design[:, 0] == 1)


def test_bootstrap_sign_fraction_is_high_for_strong_signal():
    rng = np.random.default_rng(3)
    x = np.arange(20, dtype=float)
    y = 2 * x + rng.normal(0, 0.1, 20)
    value = bootstrap_sign_fraction(
        x, y, np.repeat(["a", "b"], 10), 1.0, 50, rng
    )
    assert value > 0.95


def test_baseline_gating_diagnostic_tracks_accessibility_and_effects():
    metadata = pd.DataFrame(
        {
            "perturbation_name": ["NT_1"] * 4,
            "Timepoint": ["d1", "d2", "d3", "d4"],
            "nCount_RNA": [100] * 4,
        }
    )
    genes = sparse.csr_matrix([[1], [1], [1], [1]])
    atac = sparse.csr_matrix([[0], [1], [1], [2]])
    candidates = pd.DataFrame({"target": ["G"], "gene_index": [0]})
    evidence = pd.DataFrame(
        {
            "target": ["G"],
            "peak_index": [0],
            "effect_d1": [0.1],
            "effect_d2": [0.2],
            "effect_d3": [0.3],
            "effect_d4": [0.4],
            "rna_effect_d1": [0.1],
            "rna_effect_d2": [0.2],
            "rna_effect_d3": [0.3],
            "rna_effect_d4": [0.4],
        }
    )
    result = add_baseline_gating_diagnostics(
        evidence, genes, atac, metadata, candidates, ["NT_1"], ["d1", "d2", "d3", "d4"]
    )
    assert result.loc[0, "baseline_access_vs_absolute_atac_effect_correlation"] > 0
