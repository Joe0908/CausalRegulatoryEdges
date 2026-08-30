"""Audit GSE274113 metadata and day-14 10x multiome matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_metadata(path: Path) -> tuple[pd.DataFrame, int]:
    raw = pd.read_csv(path, index_col=0)
    blank_rows = int(raw["replicate"].isna().sum())
    metadata = raw.dropna(subset=["replicate"]).copy()
    if not metadata.index.is_unique:
        raise ValueError("Metadata cell identifiers are not unique")
    return metadata, blank_rows


def inspect_h5(path: Path, expected_barcodes: set[str]) -> dict:
    with h5py.File(path) as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        feature_types = np.char.decode(matrix["features"]["feature_type"][:])
        labels, counts = np.unique(feature_types, return_counts=True)
        gene_mask = feature_types == "Gene Expression"
        gene_names = np.char.decode(matrix["features"]["name"][:])[gene_mask]
        gene_ids = np.char.decode(matrix["features"]["id"][:])[gene_mask]
        gene_signature = hashlib.sha256(
            "\n".join(f"{gene_id}\t{name}" for gene_id, name in zip(gene_ids, gene_names)).encode()
        ).hexdigest()
        observed = {value.decode() for value in matrix["barcodes"][:]}
        nnz = int(matrix["data"].shape[0])

    overlap = observed & expected_barcodes
    return {
        "path": str(path),
        "shape_features_by_cells": list(shape),
        "nonzero_entries": nnz,
        "feature_types": {str(k): int(v) for k, v in zip(labels, counts)},
        "gene_feature_signature": gene_signature,
        "h5_barcodes": len(observed),
        "metadata_barcodes": len(expected_barcodes),
        "intersection": len(overlap),
        "metadata_missing_from_h5": len(expected_barcodes - observed),
        "h5_cells_excluded_by_metadata_qc": len(observed - expected_barcodes),
    }


def audit_all_timepoints(config_path: Path, output_dir: Path) -> dict:
    """Verify every public library against the author-QC metadata whitelist."""
    config = load_config(config_path)
    metadata, blank_rows = load_metadata(Path(config["data"]["metadata"]))
    mapping = config["data"]["replicate_timepoints"]
    expected_map = {
        replicate: timepoint
        for timepoint, replicates in mapping.items()
        for replicate in replicates
    }
    observed_replicates = set(metadata.replicate.astype(str))
    if observed_replicates != set(expected_map):
        raise ValueError(
            "Replicate mapping mismatch: "
            f"metadata-only={sorted(observed_replicates - set(expected_map))}, "
            f"config-only={sorted(set(expected_map) - observed_replicates)}"
        )

    audits = []
    for replicate in sorted(expected_map, key=lambda x: int(x.removeprefix("rep"))):
        rep_metadata = metadata.loc[metadata.replicate.eq(replicate)]
        observed_timepoints = set(rep_metadata.Timepoint.astype(str))
        if observed_timepoints != {expected_map[replicate]}:
            raise ValueError(
                f"{replicate} timepoint mismatch: {sorted(observed_timepoints)}"
            )
        expected_barcodes = {
            cell.split("_", 1)[1] for cell in rep_metadata.index.astype(str)
        }
        path = Path(config["data"]["h5_template"].format(replicate=replicate))
        result = inspect_h5(path, expected_barcodes)
        result.update(
            {
                "replicate": replicate,
                "timepoint": expected_map[replicate],
                "qc_metadata_cells": int(len(rep_metadata)),
            }
        )
        audits.append(result)

    gene_signatures = {row["gene_feature_signature"] for row in audits}
    missing_cells = sum(row["metadata_missing_from_h5"] for row in audits)
    result = {
        "metadata_rows_raw": int(len(metadata) + blank_rows),
        "blank_rows_removed": blank_rows,
        "metadata_cells_qc_pass": int(len(metadata)),
        "libraries_audited": len(audits),
        "all_metadata_barcodes_recovered": missing_cells == 0,
        "metadata_barcodes_missing_total": int(missing_cells),
        "gene_features_identical_across_libraries": len(gene_signatures) == 1,
        "gene_feature_signatures": sorted(gene_signatures),
        "timepoint_cell_counts": {
            str(k): int(v) for k, v in metadata.Timepoint.value_counts().items()
        },
        "libraries": audits,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audits).to_csv(output_dir / "all_library_audit.csv", index=False)
    pd.crosstab(metadata.Timepoint, metadata.new_CellType).to_csv(
        output_dir / "timepoint_by_cell_type.csv"
    )
    pd.crosstab(metadata.Timepoint, metadata.target).to_csv(
        output_dir / "timepoint_by_target.csv"
    )
    with (output_dir / "all_timepoint_audit.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    return result


def guide_coverage(day14: pd.DataFrame, replicates: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for (target, guide), group in day14.groupby(
        ["target", "perturbation_name"], observed=True
    ):
        counts = group.groupby("replicate", observed=True).size().to_dict()
        rows.append(
            {
                "target": target,
                "guide": guide,
                "cells_total": int(len(group)),
                **{f"cells_{rep}": int(counts.get(rep, 0)) for rep in replicates},
                "minimum_cells_across_replicates": int(
                    min(counts.get(rep, 0) for rep in replicates)
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["target", "cells_total"], ascending=[True, False]
    )


def target_summary(coverage: pd.DataFrame, config: dict) -> pd.DataFrame:
    minimum_total = int(config["mvp"]["minimum_cells_per_guide_total"])
    minimum_per_rep = int(config["mvp"]["minimum_cells_per_guide_per_replicate"])
    eligible = coverage.assign(
        representation_eligible=lambda x: (x.cells_total >= minimum_total)
        & (x.minimum_cells_across_replicates >= minimum_per_rep)
    )
    rows = []
    for target, group in eligible.groupby("target", observed=True):
        rows.append(
            {
                "target": target,
                "cells_total": int(group.cells_total.sum()),
                "guides_total": int(len(group)),
                "guides_representation_eligible": int(
                    group.representation_eligible.sum()
                ),
                "smallest_guide_total": int(group.cells_total.min()),
                "smallest_guide_replicate_count": int(
                    group.minimum_cells_across_replicates.min()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["guides_representation_eligible", "cells_total"], ascending=False
    )


def audit(config_path: Path, output_dir: Path) -> dict:
    config = load_config(config_path)
    metadata, blank_rows = load_metadata(Path(config["data"]["metadata"]))
    timepoint = config["mvp"]["timepoint"]
    replicates = list(config["data"]["day14_replicates"])
    subset = metadata.loc[metadata["Timepoint"].eq(timepoint)].copy()
    unexpected_replicates = sorted(set(subset["replicate"]) - set(replicates))
    if unexpected_replicates:
        raise ValueError(f"Unexpected {timepoint} replicates: {unexpected_replicates}")

    coverage = guide_coverage(subset, replicates)
    targets = target_summary(coverage, config)
    minimum_guides = int(config["mvp"]["minimum_guides_per_target"])
    shortlist_size = int(config["mvp"]["shortlist_size"])
    shortlist = targets.loc[
        (targets.target != config["data"]["non_targeting_label"])
        & (targets.guides_representation_eligible >= minimum_guides)
    ].head(shortlist_size)

    h5_audits = []
    for replicate in replicates:
        path = Path(config["data"]["h5_template"].format(replicate=replicate))
        prefixed = metadata.index[metadata.replicate.eq(replicate)]
        expected = {cell.split("_", 1)[1] for cell in prefixed}
        h5_audits.append(inspect_h5(path, expected))

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output_dir / "day14_guide_coverage.csv", index=False)
    targets.to_csv(output_dir / "day14_target_summary.csv", index=False)
    shortlist.to_csv(output_dir / "day14_representation_shortlist.csv", index=False)

    result = {
        "metadata_rows_raw": int(len(metadata) + blank_rows),
        "blank_rows_removed": blank_rows,
        "metadata_cells_qc_pass": int(len(metadata)),
        "timepoint_cells": int(len(subset)),
        "replicate_cell_counts": {
            str(k): int(v) for k, v in subset.replicate.value_counts().items()
        },
        "cell_type_counts": {
            str(k): int(v) for k, v in subset.new_CellType.value_counts().items()
        },
        "representation_shortlist": shortlist.target.tolist(),
        "shortlist_is_not_efficacy_validated": True,
        "h5_audits": h5_audits,
    }
    with (output_dir / "audit_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument("--output", type=Path, default=Path("reports/metadata_audit"))
    parser.add_argument(
        "--all-timepoints",
        action="store_true",
        help="Audit all configured libraries instead of only the day-14 MVP",
    )
    args = parser.parse_args()
    result = (
        audit_all_timepoints(args.config, args.output)
        if args.all_timepoints
        else audit(args.config, args.output)
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
