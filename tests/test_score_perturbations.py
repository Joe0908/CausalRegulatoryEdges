import numpy as np
from scipy import sparse

from edge_causality.score_perturbations import bh_adjust, log_normalize


def test_bh_adjust_is_monotone_in_rank() -> None:
    p = np.array([0.04, 0.001, 0.02])
    q = bh_adjust(p)
    ranked = q[np.argsort(p)]
    assert np.all(np.diff(ranked) >= 0)
    assert np.all((q >= 0) & (q <= 1))


def test_log_normalize_equalizes_nonzero_library_totals_before_log() -> None:
    counts = sparse.csr_matrix(np.array([[1, 1], [2, 2]], dtype=np.int32))
    observed = log_normalize(counts, scale=10.0).toarray()
    assert np.allclose(observed[0], observed[1])
