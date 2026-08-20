from pathlib import Path

import pandas as pd

from edge_causality.audit_metadata import guide_coverage, inspect_h5, load_metadata


def test_load_metadata_removes_blank_rows(tmp_path: Path) -> None:
    path = tmp_path / "metadata.csv.gz"
    frame = pd.DataFrame(
        {"replicate": ["rep13", None], "target": ["NT", None]},
        index=["rep13_AAAC-1", "blank"],
    )
    frame.to_csv(path, compression="gzip")
    observed, blanks = load_metadata(path)
    assert blanks == 1
    assert observed.index.tolist() == ["rep13_AAAC-1"]


def test_guide_coverage_fills_absent_replicates_with_zero() -> None:
    frame = pd.DataFrame(
        {
            "target": ["GATA1", "GATA1"],
            "perturbation_name": ["GATA1_a", "GATA1_a"],
            "replicate": ["rep13", "rep14"],
        }
    )
    observed = guide_coverage(frame, ["rep13", "rep14", "rep16"])
    assert observed.loc[0, "cells_rep16"] == 0
    assert observed.loc[0, "minimum_cells_across_replicates"] == 0


def test_inspect_h5_checks_barcodes_and_gene_signature(tmp_path: Path) -> None:
    import h5py
    import numpy as np

    path = tmp_path / "mini.h5"
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("shape", data=np.array([2, 2], dtype=np.int32))
        matrix.create_dataset("data", data=np.array([1], dtype=np.int32))
        matrix.create_dataset("barcodes", data=np.array([b"AA-1", b"BB-1"]))
        features = matrix.create_group("features")
        features.create_dataset("feature_type", data=np.array([b"Gene Expression", b"Peaks"]))
        features.create_dataset("name", data=np.array([b"GATA1", b"chr1:1-2"]))
        features.create_dataset("id", data=np.array([b"ENSG1", b"chr1:1-2"]))
    observed = inspect_h5(path, {"AA-1"})
    assert observed["metadata_missing_from_h5"] == 0
    assert observed["h5_cells_excluded_by_metadata_qc"] == 1
    assert observed["feature_types"] == {"Gene Expression": 1, "Peaks": 1}
    assert len(observed["gene_feature_signature"]) == 64
