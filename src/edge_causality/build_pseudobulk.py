"""Build RNA pseudobulks directly from filtered 10x multiome H5 files.

The author-provided metadata is the cell whitelist. Aggregation is stratified by
replicate, guide, and cell type so the initial guide screen does not mistake gross
cell-state composition shifts for a molecular perturbation signal.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def aggregate_replicate(
    h5_path: Path,
    metadata: pd.DataFrame,
    cell_types: list[str] | None,
    group_by_cell_type: bool = True,
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray, np.ndarray]:
    """Return group-by-gene counts and group metadata for one replicate."""
    if cell_types is None:
        selected = metadata.copy()
    else:
        selected = metadata.loc[metadata.new_CellType.isin(cell_types)].copy()
    selected["barcode"] = selected.index.str.split("_", n=1).str[1]
    if group_by_cell_type:
        selected["group"] = (
            selected["perturbation_name"].astype(str)
            + "||"
            + selected["new_CellType"].astype(str)
        )
    else:
        selected["group"] = selected["perturbation_name"].astype(str)

    with h5py.File(h5_path) as handle:
        matrix = handle["matrix"]
        barcodes = np.char.decode(matrix["barcodes"][:])
        barcode_to_col = {barcode: i for i, barcode in enumerate(barcodes)}
        missing = sorted(set(selected.barcode) - set(barcode_to_col))
        if missing:
            raise ValueError(f"{len(missing)} metadata barcodes missing from {h5_path}")

        group_names = sorted(selected.group.unique())
        group_to_row = {name: i for i, name in enumerate(group_names)}
        cell_group = np.full(len(barcodes), -1, dtype=np.int32)
        for row in selected.itertuples():
            cell_group[barcode_to_col[row.barcode]] = group_to_row[row.group]

        feature_types = np.char.decode(matrix["features"]["feature_type"][:])
        gene_rows = np.flatnonzero(feature_types == "Gene Expression")
        if not np.array_equal(gene_rows, np.arange(len(gene_rows))):
            raise ValueError("Gene Expression features are not the leading H5 rows")
        gene_names = np.char.decode(matrix["features"]["name"][: len(gene_rows)])
        gene_ids = np.char.decode(matrix["features"]["id"][: len(gene_rows)])

        indptr = matrix["indptr"][:]
        output = np.zeros((len(group_names), len(gene_rows)), dtype=np.int64)
        block_cells = 256
        for start_col in range(0, len(barcodes), block_cells):
            stop_col = min(start_col + block_cells, len(barcodes))
            groups = cell_group[start_col:stop_col]
            keep_columns = np.flatnonzero(groups >= 0)
            if not len(keep_columns):
                continue
            data_start = int(indptr[start_col])
            data_stop = int(indptr[stop_col])
            block_indices = matrix["indices"][data_start:data_stop]
            block_data = matrix["data"][data_start:data_stop]
            block_ptr = indptr[start_col : stop_col + 1] - data_start
            block = sparse.csc_matrix(
                (block_data, block_indices, block_ptr),
                shape=(int(matrix["shape"][0]), stop_col - start_col),
            )[: len(gene_rows), keep_columns]
            local_groups = groups[keep_columns]
            assignment = sparse.csr_matrix(
                (
                    np.ones(len(local_groups), dtype=np.int8),
                    (local_groups, np.arange(len(local_groups))),
                ),
                shape=(len(group_names), len(local_groups)),
            )
            output += (assignment @ block.T).toarray().astype(np.int64, copy=False)

    counts = selected.groupby("group", observed=True).size()
    if group_by_cell_type:
        guide, cell_type = zip(*(name.split("||", 1) for name in group_names))
    else:
        guide = group_names
        cell_type = ["ALL"] * len(group_names)
    group_metadata = pd.DataFrame(
        {
            "group": group_names,
            "guide": guide,
            "cell_type": cell_type,
            "n_cells": [int(counts[name]) for name in group_names],
        }
    )
    return output, group_metadata, gene_names, gene_ids


def build_all_timepoint_total(config_path: Path, output_dir: Path) -> None:
    """Aggregate all author-QC cells by guide and library at every timepoint.

    This is the primary time-resolved data layer. Because collection time is an
    exogenous context, it avoids conditioning the primary interaction test on a
    cell label that may itself have changed after perturbation.
    """
    config = load_config(config_path)
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0).dropna(
        subset=["replicate"]
    )
    ordered_timepoints = list(config["state_dependence"]["ordered_timepoints"])
    replicate_timepoint = {
        replicate: timepoint
        for timepoint, replicates in config["data"]["replicate_timepoints"].items()
        for replicate in replicates
    }

    matrices = []
    annotations = []
    reference_names = reference_ids = None
    for timepoint in ordered_timepoints:
        for replicate in config["data"]["replicate_timepoints"][timepoint]:
            rep_metadata = metadata.loc[metadata.replicate.eq(replicate)]
            h5_path = Path(
                config["data"]["h5_template"].format(replicate=replicate)
            )
            counts, groups, gene_names, gene_ids = aggregate_replicate(
                h5_path,
                rep_metadata,
                cell_types=None,
                group_by_cell_type=False,
            )
            if reference_names is None:
                reference_names, reference_ids = gene_names, gene_ids
            elif not (
                np.array_equal(reference_names, gene_names)
                and np.array_equal(reference_ids, gene_ids)
            ):
                raise ValueError("Gene features differ across replicate H5 files")
            groups.insert(0, "timepoint", replicate_timepoint[replicate])
            groups.insert(0, "replicate", replicate)
            guide_target = (
                rep_metadata[["perturbation_name", "target"]]
                .drop_duplicates()
                .set_index("perturbation_name")["target"]
            )
            groups["target"] = groups.guide.map(guide_target)
            matrices.append(counts)
            annotations.append(groups)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "all_timepoint_total_rna_pseudobulk_counts.npz",
        counts=np.vstack(matrices),
        gene_name=reference_names,
        gene_id=reference_ids,
    )
    pd.concat(annotations, ignore_index=True).to_csv(
        output_dir / "all_timepoint_total_rna_pseudobulk_groups.csv", index=False
    )


def build(config_path: Path, output_dir: Path) -> None:
    config = load_config(config_path)
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0).dropna(
        subset=["replicate"]
    )
    metadata = metadata.loc[metadata.Timepoint.eq(config["mvp"]["timepoint"])]
    cell_types = list(config["mvp"]["cell_types"])

    matrices = []
    annotations = []
    reference_names = reference_ids = None
    for replicate in config["data"]["day14_replicates"]:
        rep_metadata = metadata.loc[metadata.replicate.eq(replicate)]
        h5_path = Path(
            config["data"]["h5_template"].format(replicate=replicate)
        )
        counts, groups, gene_names, gene_ids = aggregate_replicate(
            h5_path, rep_metadata, cell_types
        )
        if reference_names is None:
            reference_names, reference_ids = gene_names, gene_ids
        elif not (
            np.array_equal(reference_names, gene_names)
            and np.array_equal(reference_ids, gene_ids)
        ):
            raise ValueError("Gene features differ across replicate H5 files")
        groups.insert(0, "replicate", replicate)
        matrices.append(counts)
        annotations.append(groups)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "day14_rna_pseudobulk_counts.npz",
        counts=np.vstack(matrices),
        gene_name=reference_names,
        gene_id=reference_ids,
    )
    pd.concat(annotations, ignore_index=True).to_csv(
        output_dir / "day14_rna_pseudobulk_groups.csv", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/pseudobulk"))
    parser.add_argument(
        "--all-timepoints-total",
        action="store_true",
        help="Aggregate all QC cells by guide and library across all timepoints",
    )
    args = parser.parse_args()
    if args.all_timepoints_total:
        build_all_timepoint_total(args.config, args.output)
    else:
        build(args.config, args.output)


if __name__ == "__main__":
    main()
