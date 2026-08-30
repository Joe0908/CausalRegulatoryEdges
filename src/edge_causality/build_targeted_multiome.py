"""Build a compact paired RNA/ATAC matrix around frozen candidate genes.

Peak calls differ between 10x libraries.  We therefore collect all local peak
intervals in a fixed TSS window and merge overlapping intervals across
libraries before extracting counts.  Only QC-passing cells needed for the
candidate perturbations and controls are retained.
"""

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


def parse_interval(value: str) -> tuple[str, int, int]:
    chrom, coordinates = str(value).split(":", 1)
    start, end = coordinates.split("-", 1)
    return chrom, int(start), int(end)


def merge_intervals(
    intervals: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Return the union of overlapping (but not merely adjacent) intervals."""
    merged: list[list[object]] = []
    for chrom, start, end in sorted(intervals, key=lambda x: (x[0], x[1], x[2])):
        if not merged or merged[-1][0] != chrom or start >= int(merged[-1][2]):
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(int(merged[-1][2]), end)
    return [(str(chrom), int(start), int(end)) for chrom, start, end in merged]


def local_peak_rows(
    feature_type: np.ndarray,
    intervals: np.ndarray,
    chromosome: str,
    window_start: int,
    window_end: int,
) -> list[tuple[int, str, int, int]]:
    output = []
    for row in np.flatnonzero(feature_type == "Peaks"):
        chrom, start, end = parse_interval(str(intervals[row]))
        if chrom == chromosome and start < window_end and end > window_start:
            output.append((int(row), chrom, start, end))
    return output


def build_consensus_peaks(
    h5_paths: dict[str, Path], candidates: pd.DataFrame, window_bp: int
) -> tuple[pd.DataFrame, dict[str, dict[int, int]]]:
    """Create candidate-specific consensus intervals and source-row mappings."""
    source: dict[str, dict[str, list[tuple[int, str, int, int]]]] = {}
    for replicate, path in h5_paths.items():
        with h5py.File(path) as handle:
            features = handle["matrix/features"]
            feature_type = np.char.decode(features["feature_type"][:])
            intervals = np.char.decode(features["interval"][:])
        source[replicate] = {}
        for row in candidates.itertuples(index=False):
            source[replicate][row.target] = local_peak_rows(
                feature_type,
                intervals,
                row.chromosome,
                max(0, int(row.tss) - window_bp),
                int(row.tss) + window_bp,
            )

    records: list[dict] = []
    lookup: dict[tuple[str, str, int, int], int] = {}
    for candidate in candidates.itertuples(index=False):
        intervals = [
            (chrom, start, end)
            for replicate in h5_paths
            for _, chrom, start, end in source[replicate][candidate.target]
        ]
        for chrom, start, end in merge_intervals(intervals):
            peak_index = len(records)
            lookup[(candidate.target, chrom, start, end)] = peak_index
            records.append(
                {
                    "peak_index": peak_index,
                    "peak_id": f"{chrom}:{start}-{end}",
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "TF": candidate.TF,
                    "target": candidate.target,
                    "gene_id": candidate.gene_id,
                    "tss": int(candidate.tss),
                    "distance_to_tss": int((start + end) // 2 - candidate.tss),
                    "candidate_role": candidate.role,
                }
            )
    peaks = pd.DataFrame(records)

    row_maps: dict[str, dict[int, int]] = {}
    presence: dict[int, set[str]] = {i: set() for i in peaks.peak_index}
    for replicate in h5_paths:
        row_maps[replicate] = {}
        for candidate in candidates.itertuples(index=False):
            candidate_peaks = peaks.loc[peaks.target.eq(candidate.target)]
            for source_row, chrom, start, end in source[replicate][candidate.target]:
                hit = candidate_peaks.loc[
                    candidate_peaks.chromosome.eq(chrom)
                    & candidate_peaks.start.lt(end)
                    & candidate_peaks.end.gt(start)
                ]
                if len(hit) != 1:
                    raise ValueError(
                        f"Source peak maps to {len(hit)} consensus peaks: "
                        f"{replicate} {chrom}:{start}-{end}"
                    )
                consensus = int(hit.iloc[0].peak_index)
                if source_row in row_maps[replicate]:
                    raise ValueError("A source peak was assigned to multiple genes")
                row_maps[replicate][source_row] = consensus
                presence[consensus].add(replicate)
    peaks["libraries_present"] = peaks.peak_index.map(
        lambda i: len(presence[int(i)])
    )
    peaks["replicates_present"] = peaks.peak_index.map(
        lambda i: ";".join(sorted(presence[int(i)], key=lambda x: int(x[3:])))
    )
    return peaks, row_maps


def extract_replicate(
    h5_path: Path,
    replicate: str,
    requested_barcodes: set[str],
    candidates: pd.DataFrame,
    row_to_consensus: dict[int, int],
    n_consensus: int,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, list[str]]:
    """Extract candidate gene and harmonized ATAC counts for one library."""
    with h5py.File(h5_path) as handle:
        matrix = handle["matrix"]
        barcodes = np.char.decode(matrix["barcodes"][:])
        selected_columns = np.flatnonzero(
            np.fromiter((b in requested_barcodes for b in barcodes), dtype=bool)
        )
        selected_barcodes = [str(barcodes[i]) for i in selected_columns]
        missing = requested_barcodes - set(selected_barcodes)
        if missing:
            raise ValueError(f"{len(missing)} QC cells absent from {h5_path}")

        names = np.char.decode(matrix["features/name"][:])
        ids = np.char.decode(matrix["features/id"][:])
        gene_row_map: dict[int, int] = {}
        for gene_index, candidate in enumerate(candidates.itertuples(index=False)):
            hit = np.flatnonzero(
                (names == candidate.target) & (ids == candidate.gene_id)
            )
            if len(hit) != 1:
                raise ValueError(
                    f"Expected one feature for {candidate.target}/{candidate.gene_id}"
                )
            gene_row_map[int(hit[0])] = gene_index

        row_to_gene = np.full(len(names), -1, dtype=np.int32)
        for source_row, gene_index in gene_row_map.items():
            row_to_gene[source_row] = gene_index
        row_to_peak = np.full(len(names), -1, dtype=np.int32)
        for source_row, peak_index in row_to_consensus.items():
            row_to_peak[source_row] = peak_index

        source_indptr = matrix["indptr"][:]
        gene_cell: list[int] = []
        gene_feature: list[int] = []
        gene_data: list[float] = []
        peak_cell: list[int] = []
        peak_feature: list[int] = []
        peak_data: list[float] = []
        for output_cell, column in enumerate(selected_columns):
            start, stop = int(source_indptr[column]), int(source_indptr[column + 1])
            source_rows = matrix["indices"][start:stop]
            values = matrix["data"][start:stop]
            mapped_gene = row_to_gene[source_rows]
            keep_gene = mapped_gene >= 0
            gene_cell.extend([output_cell] * int(keep_gene.sum()))
            gene_feature.extend(mapped_gene[keep_gene].tolist())
            gene_data.extend(values[keep_gene].tolist())
            mapped_peak = row_to_peak[source_rows]
            keep_peak = mapped_peak >= 0
            peak_cell.extend([output_cell] * int(keep_peak.sum()))
            peak_feature.extend(mapped_peak[keep_peak].tolist())
            peak_data.extend(values[keep_peak].tolist())

    genes = sparse.coo_matrix(
        (gene_data, (gene_cell, gene_feature)),
        shape=(len(selected_columns), len(candidates)),
    ).tocsr()
    atac = sparse.coo_matrix(
        (peak_data, (peak_cell, peak_feature)),
        shape=(len(selected_columns), n_consensus),
    ).tocsr()
    cell_ids = [f"{replicate}_{barcode}" for barcode in selected_barcodes]
    return genes, atac, cell_ids


def build(config_path: Path, output_dir: Path) -> dict:
    config = load_config(config_path)
    settings = config["chromatin_mechanism"]
    candidates = pd.DataFrame(settings["candidates"])
    replicate_order = [
        rep
        for reps in config["data"]["replicate_timepoints"].values()
        for rep in reps
    ]
    h5_paths = {
        rep: Path(config["data"]["h5_template"].format(replicate=rep))
        for rep in replicate_order
    }
    missing_h5 = [str(path) for path in h5_paths.values() if not path.exists()]
    if missing_h5:
        raise FileNotFoundError(f"Missing {len(missing_h5)} H5 files: {missing_h5[:2]}")

    peaks, row_maps = build_consensus_peaks(
        h5_paths, candidates, int(settings["window_bp"])
    )
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0).dropna(
        subset=["replicate"]
    )
    target_tfs = candidates.TF.unique().tolist()
    control_guides = list(
        dict.fromkeys(
            config["data"]["discovery_control_guides"]
            + config["data"]["intervention_reference_guides"]
        )
    )
    metadata = metadata.loc[
        metadata.target.isin(target_tfs)
        | metadata.perturbation_name.isin(control_guides)
    ].copy()

    gene_blocks = []
    atac_blocks = []
    cell_ids: list[str] = []
    for replicate in replicate_order:
        selected = metadata.loc[metadata.replicate.eq(replicate)]
        requested = {str(cell).split("_", 1)[1] for cell in selected.index}
        genes, atac, ids = extract_replicate(
            h5_paths[replicate],
            replicate,
            requested,
            candidates,
            row_maps[replicate],
            len(peaks),
        )
        gene_blocks.append(genes)
        atac_blocks.append(atac)
        cell_ids.extend(ids)

    genes = sparse.vstack(gene_blocks, format="csr")
    atac = sparse.vstack(atac_blocks, format="csr")
    cell_metadata = metadata.loc[cell_ids].copy()
    cell_metadata.insert(0, "cell_id", cell_metadata.index)
    cell_metadata.reset_index(drop=True, inplace=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(output_dir / "candidate_rna_counts.npz", genes)
    sparse.save_npz(output_dir / "candidate_atac_counts.npz", atac)
    cell_metadata.to_csv(
        output_dir / "candidate_cell_metadata.csv.gz", index=False, compression="gzip"
    )
    candidates.insert(0, "gene_index", np.arange(len(candidates)))
    candidates.to_csv(output_dir / "candidate_genes.csv", index=False)
    peaks.to_csv(output_dir / "candidate_consensus_peaks.csv.gz", index=False)
    summary = {
        "cells": int(len(cell_metadata)),
        "candidate_genes": int(len(candidates)),
        "consensus_peaks": int(len(peaks)),
        "rna_nonzero": int(genes.nnz),
        "atac_nonzero": int(atac.nnz),
        "peak_libraries_present": {
            "minimum": int(peaks.libraries_present.min()),
            "median": float(peaks.libraries_present.median()),
            "maximum": int(peaks.libraries_present.max()),
        },
        "cells_by_timepoint": {
            str(k): int(v) for k, v in cell_metadata.Timepoint.value_counts().items()
        },
        "cells_by_perturbation_class": {
            "GATA1": int(cell_metadata.target.eq("GATA1").sum()),
            "NFE2": int(cell_metadata.target.eq("NFE2").sum()),
            "controls": int(cell_metadata.perturbation_name.isin(control_guides).sum()),
        },
    }
    with (output_dir / "build_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--output", type=Path, default=Path("data/processed/targeted_multiome")
    )
    args = parser.parse_args()
    print(json.dumps(build(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
