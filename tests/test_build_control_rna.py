import numpy as np
from scipy import sparse

from edge_causality.build_control_rna import select_candidate_genes


def test_select_candidate_genes_applies_detection_cpm_and_prefix_filters() -> None:
    counts = sparse.csr_matrix(
        np.array(
            [
                [10, 0, 5],
                [10, 0, 5],
                [10, 1, 5],
                [10, 0, 5],
            ]
        )
    )
    observed = select_candidate_genes(
        counts,
        np.array(["GATA1", "RARE", "MT-CO1"]),
        minimum_detection_fraction=0.5,
        minimum_mean_cpm=1.0,
        excluded_prefixes=["MT-"],
    )
    assert observed.candidate_eligible.tolist() == [True, False, False]
