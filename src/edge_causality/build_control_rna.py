"""Build leakage-free control RNA matrices from 10x multiome H5 files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def extract_rna_cells(
    h5_path: Path,
    requested_barcodes: set[str],
    replicate: str,
) -> tuple[sparse.csr_matrix, list[str], np.ndarray, np.ndarray]:
    """Extract requested QC cells, preserving H5 column order."""
    with h5py.File(h5_path) as handle:
        matrix = handle["matrix"]
        barcodes = np.char.decode(matrix["barcodes"][:])
        selected_columns = [
            i for i, barcode in enumerate(barcodes) if barcode in requested_barcodes
        ]
        selected_barcodes = [str(barcodes[i]) for i in selected_columns]
        missing = requested_barcodes - set(selected_barcodes)
        if missing:
            raise ValueError(f"{len(missing)} requested cells missing from {h5_path}")

        feature_types = np.char.decode(matrix["features"]["feature_type"][:])
        gene_rows = np.flatnonzero(feature_types == "Gene Expression")
        if not np.array_equal(gene_rows, np.arange(len(gene_rows))):
            raise ValueError("Gene Expression features must be the leading H5 rows")
        n_genes = len(gene_rows)
        gene_names = np.char.decode(matrix["features"]["name"][:n_genes])
        gene_ids = np.char.decode(matrix["features"]["id"][:n_genes])
        source_indptr = matrix["indptr"][:]

        data_parts: list[np.ndarray] = []
        index_parts: list[np.ndarray] = []
        output_indptr = [0]
        for column in selected_columns:
            start, stop = int(source_indptr[column]), int(source_indptr[column + 1])
            indices = matrix["indices"][start:stop]
            data = matrix["data"][start:stop]
            keep = indices < n_genes
            index_parts.append(indices[keep].astype(np.int32, copy=False))
            data_parts.append(data[keep])
            output_indptr.append(output_indptr[-1] + int(keep.sum()))

    if data_parts:
        data = np.concatenate(data_parts)
        indices = np.concatenate(index_parts)
    else:
        data = np.array([], dtype=np.int32)
        indices = np.array([], dtype=np.int32)
    gene_by_cell = sparse.csc_matrix(
        (data, indices, np.asarray(output_indptr, dtype=np.int64)),
        shape=(n_genes, len(selected_columns)),
    )
    cell_ids = [f"{replicate}_{barcode}" for barcode in selected_barcodes]
    return gene_by_cell.T.tocsr(), cell_ids, gene_names, gene_ids


def select_candidate_genes(
    counts: sparse.csr_matrix,
    gene_names: np.ndarray,
    minimum_detection_fraction: float,
    minimum_mean_cpm: float,
    excluded_prefixes: list[str],
) -> pd.DataFrame:
    n_cells = counts.shape[0]
    detection_fraction = np.asarray((counts > 0).sum(axis=0)).ravel() / n_cells
    library_sizes = np.asarray(counts.sum(axis=1)).ravel()
    inverse_sizes = np.divide(
        1_000_000.0,
        library_sizes,
        out=np.zeros_like(library_sizes, dtype=float),
        where=library_sizes > 0,
    )
    mean_cpm = np.asarray(counts.multiply(inverse_sizes[:, None]).mean(axis=0)).ravel()
    excluded = np.zeros(len(gene_names), dtype=bool)
    for prefix in excluded_prefixes:
        excluded |= np.char.startswith(gene_names.astype(str), prefix)
    eligible = (
        (detection_fraction >= minimum_detection_fraction)
        & (mean_cpm >= minimum_mean_cpm)
        & ~excluded
    )
    return pd.DataFrame(
        {
            "gene_name": gene_names,
            "detection_fraction": detection_fraction,
            "mean_cpm": mean_cpm,
            "excluded_prefix": excluded,
            "candidate_eligible": eligible,
        }
    )


def build(
    config_path: Path,
    output_dir: Path,
    all_timepoints: bool,
    all_guides: bool = False,
) -> dict:
    config = load_config(config_path)
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0).dropna(
        subset=["replicate"]
    )
    if not all_guides:
        metadata = metadata.loc[
            metadata.perturbation_name.isin(
                config["data"]["discovery_control_guides"]
            )
        ].copy()
    if not all_timepoints:
        metadata = metadata.loc[
            metadata.Timepoint.eq(config["mvp"]["timepoint"])
            & metadata.new_CellType.isin(config["mvp"]["cell_types"])
        ].copy()

    matrices = []
    cell_ids: list[str] = []
    reference_names = reference_ids = None
    replicate_order = sorted(
        metadata.replicate.unique(), key=lambda x: int(str(x).removeprefix("rep"))
    )
    for replicate in replicate_order:
        rep_metadata = metadata.loc[metadata.replicate.eq(replicate)]
        requested = {cell.split("_", 1)[1] for cell in rep_metadata.index.astype(str)}
        path = Path(config["data"]["h5_template"].format(replicate=replicate))
        counts, ids, names, gene_ids = extract_rna_cells(path, requested, replicate)
        if reference_names is None:
            reference_names, reference_ids = names, gene_ids
        elif not (
            np.array_equal(reference_names, names)
            and np.array_equal(reference_ids, gene_ids)
        ):
            raise ValueError("Gene features differ across libraries")
        matrices.append(counts)
        cell_ids.extend(ids)

    counts = sparse.vstack(matrices, format="csr")
    cell_metadata = metadata.loc[cell_ids].copy()
    discovery_rows = cell_metadata.perturbation_name.isin(
        config["data"]["discovery_control_guides"]
    ).to_numpy()
    features = select_candidate_genes(
        counts[discovery_rows],
        reference_names,
        float(config["edge_discovery"]["minimum_detection_fraction"]),
        float(config["edge_discovery"]["minimum_mean_cpm"]),
        list(config["edge_discovery"]["exclude_gene_prefixes"]),
    )
    features.insert(0, "gene_id", reference_ids)
    duplicate_symbol = features.gene_name.duplicated(keep=False)
    features.insert(
        0,
        "feature_key",
        np.where(
            duplicate_symbol,
            features.gene_name.astype(str) + "|" + features.gene_id.astype(str),
            features.gene_name.astype(str),
        ),
    )
    if features.feature_key.duplicated().any():
        raise ValueError("Feature keys remain non-unique after Ensembl disambiguation")

    output_dir.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(output_dir / "rna_counts_cells_by_genes.npz", counts)
    cell_metadata.to_csv(output_dir / "cell_metadata.csv.gz", compression="gzip")
    features.to_csv(output_dir / "gene_features.csv.gz", index=False, compression="gzip")
    summary = {
        "scope": (
            "all_timepoints" if all_timepoints else "day14_late_erythroid"
        ),
        "guide_scope": "all_guides" if all_guides else "discovery_controls_only",
        "controls_used_for_feature_selection": list(
            config["data"]["discovery_control_guides"]
        ),
        "cells": int(counts.shape[0]),
        "genes_total": int(counts.shape[1]),
        "genes_candidate_eligible": int(features.candidate_eligible.sum()),
        "duplicate_gene_symbol_rows_disambiguated": int(duplicate_symbol.sum()),
        "nonzero_counts": int(counts.nnz),
        "replicate_counts": {
            str(k): int(v) for k, v in cell_metadata.replicate.value_counts().items()
        },
        "timepoint_counts": {
            str(k): int(v) for k, v in cell_metadata.Timepoint.value_counts().items()
        },
        "cell_type_counts": {
            str(k): int(v) for k, v in cell_metadata.new_CellType.value_counts().items()
        },
    }
    with (output_dir / "build_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--all-timepoints", action="store_true")
    parser.add_argument("--all-guides", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            build(args.config, args.output, args.all_timepoints, args.all_guides),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
