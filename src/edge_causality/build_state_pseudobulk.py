"""Build candidate-gene pseudobulks for the secondary erythroid-state model."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from edge_causality.build_pseudobulk import aggregate_replicate


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def build(
    config_path: Path,
    edge_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    states = list(config["state_dependence"]["secondary_erythroid_states"])
    metadata = pd.read_csv(config["data"]["metadata"], index_col=0).dropna(
        subset=["replicate"]
    )
    edges = pd.read_csv(edge_path)
    candidate_keys = set(edges.target.astype(str))

    matrices = []
    annotations = []
    selected_indices = None
    selected_keys = selected_names = selected_ids = None
    for timepoint in config["state_dependence"]["ordered_timepoints"]:
        for replicate in config["data"]["replicate_timepoints"][timepoint]:
            rep_metadata = metadata.loc[metadata.replicate.eq(replicate)]
            path = Path(config["data"]["h5_template"].format(replicate=replicate))
            counts, groups, names, ids = aggregate_replicate(
                path, rep_metadata, states, group_by_cell_type=True
            )
            duplicate = pd.Series(names).duplicated(keep=False).to_numpy()
            keys = np.where(duplicate, names + "|" + ids, names)
            if selected_indices is None:
                selected_indices = np.flatnonzero(np.isin(keys, list(candidate_keys)))
                selected_keys = keys[selected_indices]
                selected_names = names[selected_indices]
                selected_ids = ids[selected_indices]
                missing = candidate_keys - set(selected_keys)
                if missing:
                    raise ValueError(f"Candidate targets missing from H5: {sorted(missing)[:5]}")
            elif not np.array_equal(keys[selected_indices], selected_keys):
                raise ValueError("Candidate gene features differ across libraries")
            groups.insert(0, "timepoint", timepoint)
            groups.insert(0, "replicate", replicate)
            guide_target = (
                rep_metadata[["perturbation_name", "target"]]
                .drop_duplicates()
                .set_index("perturbation_name")["target"]
            )
            groups["guide_target"] = groups.guide.map(guide_target)
            matrices.append(counts[:, selected_indices])
            annotations.append(groups)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "erythroid_state_candidate_counts.npz",
        counts=np.vstack(matrices),
        feature_key=selected_keys,
        gene_name=selected_names,
        gene_id=selected_ids,
    )
    groups = pd.concat(annotations, ignore_index=True)
    groups.to_csv(output_dir / "erythroid_state_candidate_groups.csv", index=False)
    return {
        "profiles": int(len(groups)),
        "cells": int(groups.n_cells.sum()),
        "candidate_genes": int(len(selected_keys)),
        "states": states,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--edges",
        type=Path,
        default=Path("reports/validation/E0_to_E1_edge_matrix.csv.gz"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/processed/state_pseudobulk")
    )
    args = parser.parse_args()
    print(build(args.config, args.edges, args.output))


if __name__ == "__main__":
    main()
