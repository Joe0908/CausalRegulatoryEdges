"""External temporal validation of targeted chromatin evidence.

This module intentionally does not call observational data causal evidence.  It
asks a narrower transport question: do perturbation-sensitive GATA1 peaks become
accessible during the progenitor transition before their target genes undergo
their largest erythroid induction in an independent human atlas?
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact


INTERVAL_PATTERN = re.compile(r"^(chr[A-Za-z0-9]+)[:-](\d+)-(\d+)$")


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def parse_interval(value: str) -> tuple[str, int, int]:
    match = INTERVAL_PATTERN.match(str(value))
    if match is None:
        raise ValueError(f"Cannot parse genomic interval: {value}")
    chromosome, start, end = match.groups()
    return chromosome, int(start), int(end)


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if interval_overlap(start_a, end_a, start_b, end_b) > 0:
        return 0
    if end_a <= start_b:
        return start_b - end_a
    return start_a - end_b


def frozen_peak_sets(evidence: pd.DataFrame, primary_tf: str) -> pd.DataFrame:
    """Return the positive and linked comparison peaks frozen in targeted analysis."""
    tf = evidence.loc[evidence.TF.eq(primary_tf)].copy()
    positive = tf.loc[tf.E2_peak.astype(bool)].copy()
    positive["external_set"] = "targeted_E2"
    comparison = tf.loc[
        tf.link_pass.astype(bool) & ~tf.E2_peak.astype(bool)
    ].copy()
    comparison["external_set"] = "linked_non_E2"
    selected = pd.concat([positive, comparison], ignore_index=True)
    columns = [
        "TF",
        "target",
        "peak_id",
        "chromosome",
        "start",
        "end",
        "distance_to_tss",
        "candidate_role",
        "external_set",
    ]
    return selected[columns].drop_duplicates("peak_id").reset_index(drop=True)


def _decode_csv_name(value: bytes) -> str:
    return value.strip().strip(b'"').decode("utf-8")


def extract_overlapping_peak_matrix(
    matrix_path: Path,
    candidates: pd.DataFrame,
    nearest_max_distance: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream the wide ATAC CSV and retain only intervals near frozen peaks.

    The upstream file is more than two GB.  Parsing numeric values only for
    matched rows keeps this operation memory-bounded.
    """
    by_chromosome: dict[str, list[dict]] = {}
    for row in candidates.to_dict("records"):
        by_chromosome.setdefault(str(row["chromosome"]), []).append(row)
    exact_records: list[dict] = []
    nearest: dict[str, tuple[int, str, np.ndarray, int, int]] = {}
    values: dict[str, np.ndarray] = {}
    with matrix_path.open("rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b",")
        cell_names = [_decode_csv_name(value) for value in header[1:]]
        for line in handle:
            separator = line.find(b",")
            if separator < 0:
                continue
            external_peak_id = _decode_csv_name(line[:separator])
            try:
                chromosome, start, end = parse_interval(external_peak_id)
            except ValueError:
                continue
            chromosome_candidates = by_chromosome.get(chromosome, [])
            if not chromosome_candidates:
                continue
            parsed_values: np.ndarray | None = None
            for candidate in chromosome_candidates:
                overlap = interval_overlap(
                    int(candidate["start"]), int(candidate["end"]), start, end
                )
                distance = interval_distance(
                    int(candidate["start"]), int(candidate["end"]), start, end
                )
                if overlap > 0:
                    if parsed_values is None:
                        parsed_values = np.fromstring(
                            line[separator + 1 :].decode("utf-8"), sep=","
                        )
                    values[external_peak_id] = parsed_values
                    exact_records.append(
                        {
                            "peak_id": candidate["peak_id"],
                            "external_peak_id": external_peak_id,
                            "external_start": start,
                            "external_end": end,
                            "overlap_bp": overlap,
                            "distance_bp": 0,
                            "mapping_type": "overlap",
                        }
                    )
                current = nearest.get(str(candidate["peak_id"]))
                if distance <= nearest_max_distance and (
                    current is None or distance < current[0]
                ):
                    if parsed_values is None:
                        parsed_values = np.fromstring(
                            line[separator + 1 :].decode("utf-8"), sep=","
                        )
                    nearest[str(candidate["peak_id"])] = (
                        distance,
                        external_peak_id,
                        parsed_values.copy(),
                        start,
                        end,
                    )
    mapped_candidates = {record["peak_id"] for record in exact_records}
    sensitivity_records = []
    for candidate_peak, (distance, external_peak, row_values, start, end) in nearest.items():
        if candidate_peak in mapped_candidates:
            continue
        values[external_peak] = row_values
        sensitivity_records.append(
            {
                "peak_id": candidate_peak,
                "external_peak_id": external_peak,
                "external_start": start,
                "external_end": end,
                "overlap_bp": 0,
                "distance_bp": distance,
                "mapping_type": "nearest_sensitivity",
            }
        )
    mapping = pd.DataFrame(exact_records + sensitivity_records)
    if values:
        matrix = pd.DataFrame(values, index=cell_names).T
        matrix.index.name = "external_peak_id"
    else:
        matrix = pd.DataFrame(columns=cell_names)
        matrix.index.name = "external_peak_id"
    return mapping, matrix


def standardized_state_effect(
    values: np.ndarray, states: np.ndarray, early: str, committed: str
) -> float:
    selected = np.isin(states, [early, committed]) & np.isfinite(values)
    standard_deviation = np.std(values[selected], ddof=1)
    if selected.sum() < 2 or standard_deviation == 0:
        return float("nan")
    early_mean = np.mean(values[selected & (states == early)])
    committed_mean = np.mean(values[selected & (states == committed)])
    return float((committed_mean - early_mean) / standard_deviation)


def stratified_bootstrap_state_effect(
    values: np.ndarray,
    states: np.ndarray,
    strata: np.ndarray,
    early: str,
    committed: str,
    iterations: int,
    confidence: float,
    rng: np.random.Generator,
) -> tuple[float, float]:
    selected = np.isin(states, [early, committed]) & np.isfinite(values)
    groups = [
        np.flatnonzero(selected & (states == state) & (strata == stratum))
        for state in [early, committed]
        for stratum in np.unique(strata[selected])
    ]
    groups = [group for group in groups if len(group)]
    effects = []
    for _ in range(iterations):
        sampled = np.concatenate(
            [rng.choice(group, size=len(group), replace=True) for group in groups]
        )
        effect = standardized_state_effect(
            values[sampled], states[sampled], early, committed
        )
        if np.isfinite(effect):
            effects.append(effect)
    if not effects:
        return float("nan"), float("nan")
    alpha = 1 - confidence
    return tuple(np.quantile(effects, [alpha / 2, 1 - alpha / 2]).tolist())


def atac_peak_evidence(
    candidates: pd.DataFrame,
    mapping: pd.DataFrame,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    settings: dict,
    random_seed: int,
) -> pd.DataFrame:
    early = str(settings["atac_states"]["early"])
    committed = str(settings["atac_states"]["committed"])
    primary_cells = metadata["predicted.id"].isin([early, committed])
    metadata = metadata.loc[primary_cells].copy()
    matrix = matrix.loc[:, metadata.index]
    rng = np.random.default_rng(random_seed)
    records = []
    for candidate in candidates.itertuples(index=False):
        candidate_mapping = mapping.loc[
            mapping.peak_id.eq(candidate.peak_id)
            & mapping.mapping_type.eq("overlap")
        ]
        if candidate_mapping.empty:
            records.append(
                {
                    **candidate._asdict(),
                    "mapped_external_peaks": 0,
                    "external_peak_ids": "",
                    "atac_early_cells": int(
                        metadata["predicted.id"].eq(early).sum()
                    ),
                    "atac_committed_cells": int(
                        metadata["predicted.id"].eq(committed).sum()
                    ),
                    "atac_standardized_effect": float("nan"),
                    "atac_bootstrap_ci_low": float("nan"),
                    "atac_bootstrap_ci_high": float("nan"),
                    "atac_establishment_pass": False,
                }
            )
            continue
        external_ids = candidate_mapping.external_peak_id.unique().tolist()
        external_values = matrix.loc[external_ids].to_numpy(dtype=float)
        selected_values = external_values[:, primary_cells.loc[metadata.index].to_numpy()]
        means = selected_values.mean(axis=1, keepdims=True)
        standard_deviations = selected_values.std(axis=1, ddof=1, keepdims=True)
        standard_deviations[standard_deviations == 0] = 1
        profile = ((selected_values - means) / standard_deviations).mean(axis=0)
        states = metadata["predicted.id"].astype(str).to_numpy()
        batches = metadata["batch"].astype(str).to_numpy()
        effect = standardized_state_effect(profile, states, early, committed)
        ci_low, ci_high = stratified_bootstrap_state_effect(
            profile,
            states,
            batches,
            early,
            committed,
            int(settings["bootstrap_iterations"]),
            float(settings["bootstrap_confidence"]),
            rng,
        )
        passes = np.isfinite(effect) and effect >= float(
            settings["minimum_standardized_atac_establishment_effect"]
        )
        if settings.get("require_atac_bootstrap_ci_above_zero", False):
            passes = passes and ci_low > 0
        records.append(
            {
                **candidate._asdict(),
                "mapped_external_peaks": len(external_ids),
                "external_peak_ids": ";".join(external_ids),
                "atac_early_cells": int(np.sum(states == early)),
                "atac_committed_cells": int(np.sum(states == committed)),
                "atac_early_mean_z": float(np.mean(profile[states == early])),
                "atac_committed_mean_z": float(
                    np.mean(profile[states == committed])
                ),
                "atac_standardized_effect": effect,
                "atac_bootstrap_ci_low": ci_low,
                "atac_bootstrap_ci_high": ci_high,
                "atac_establishment_pass": bool(passes),
            }
        )
    return pd.DataFrame(records)


def commitment_metrics(early: float, committed: float, late: float) -> dict:
    dynamic_range = late - early
    fraction = (
        (committed - early) / dynamic_range
        if dynamic_range > 0
        else float("nan")
    )
    late_log2_fold_change = np.log2(late + 0.5) - np.log2(committed + 0.5)
    return {
        "rna_dynamic_range_cpm": float(dynamic_range),
        "rna_commitment_fraction": float(fraction),
        "rna_late_log2_fold_change": float(late_log2_fold_change),
    }


def bootstrap_rna_late_fold_change(
    pseudobulk: pd.DataFrame,
    committed: str,
    late: str,
    iterations: int,
    confidence: float,
    rng: np.random.Generator,
) -> tuple[float, float]:
    committed_values = pseudobulk.loc[
        pseudobulk.state.eq(committed), "cpm"
    ].to_numpy(dtype=float)
    late_values = pseudobulk.loc[pseudobulk.state.eq(late), "cpm"].to_numpy(
        dtype=float
    )
    if len(committed_values) == 0 or len(late_values) == 0:
        return float("nan"), float("nan")
    effects = np.empty(iterations, dtype=float)
    for iteration in range(iterations):
        committed_mean = rng.choice(
            committed_values, len(committed_values), replace=True
        ).mean()
        late_mean = rng.choice(late_values, len(late_values), replace=True).mean()
        effects[iteration] = np.log2(late_mean + 0.5) - np.log2(
            committed_mean + 0.5
        )
    alpha = 1 - confidence
    return tuple(np.quantile(effects, [alpha / 2, 1 - alpha / 2]).tolist())


def extract_rna_evidence(
    h5ad_path: Path,
    targets: list[str],
    settings: dict,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        import anndata as ad
    except ImportError as error:  # pragma: no cover - exercised by CLI environments
        raise RuntimeError("External validation requires the 'anndata' package") from error

    adata = ad.read_h5ad(h5ad_path)
    states = settings["rna_states"]
    ordered_states = [states["early"], states["committed"], states["late"]]
    selected_cells = adata.obs["Cluster"].isin(ordered_states).to_numpy()
    selected_obs = adata.obs.loc[selected_cells].copy()
    selected_targets = [target for target in targets if target in adata.var_names]
    missing_targets = sorted(set(targets) - set(selected_targets))
    counts = adata.layers["counts"][selected_cells, :][
        :, adata.var_names.get_indexer(selected_targets)
    ]
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=float)
    totals = selected_obs["n_counts"].to_numpy(dtype=float)
    rows = []
    for state in ordered_states:
        for sample in selected_obs.loc[selected_obs.Cluster.eq(state), "sample"].unique():
            cell_mask = (
                selected_obs.Cluster.eq(state) & selected_obs["sample"].eq(sample)
            ).to_numpy()
            if not cell_mask.any():
                continue
            sample_total = float(totals[cell_mask].sum())
            for gene_index, target in enumerate(selected_targets):
                rows.append(
                    {
                        "target": target,
                        "state": state,
                        "sample": str(sample),
                        "n_cells": int(cell_mask.sum()),
                        "count": float(counts[cell_mask, gene_index].sum()),
                        "library_total": sample_total,
                        "cpm": float(
                            counts[cell_mask, gene_index].sum()
                            / max(sample_total, 1)
                            * 1_000_000
                        ),
                    }
                )
    pseudobulk = pd.DataFrame(rows)
    rng = np.random.default_rng(random_seed)
    records = []
    for target in targets:
        if target in missing_targets:
            records.append(
                {
                    "target": target,
                    "rna_gene_present": False,
                    "rna_delay_pass": False,
                }
            )
            continue
        gene = pseudobulk.loc[pseudobulk.target.eq(target)]
        means = gene.groupby("state", observed=True).cpm.mean()
        sample_counts = gene.groupby("state", observed=True)["sample"].nunique()
        early, committed, late = ordered_states
        metrics = commitment_metrics(means[early], means[committed], means[late])
        ci_low, ci_high = bootstrap_rna_late_fold_change(
            gene,
            committed,
            late,
            int(settings["bootstrap_iterations"]),
            float(settings["bootstrap_confidence"]),
            rng,
        )
        enough_samples = all(
            sample_counts.get(state, 0) >= int(settings["minimum_samples_per_rna_state"])
            for state in ordered_states
        )
        passes = (
            enough_samples
            and metrics["rna_dynamic_range_cpm"] > 0
            and metrics["rna_late_log2_fold_change"]
            >= float(settings["minimum_late_rna_log2_fold_change"])
            and metrics["rna_commitment_fraction"]
            <= float(settings["maximum_rna_commitment_fraction"])
        )
        if settings.get("require_late_rna_bootstrap_ci_above_zero", False):
            passes = passes and ci_low > 0
        records.append(
            {
                "target": target,
                "rna_gene_present": True,
                "rna_early_mean_cpm": float(means[early]),
                "rna_committed_mean_cpm": float(means[committed]),
                "rna_late_mean_cpm": float(means[late]),
                "rna_early_samples": int(sample_counts[early]),
                "rna_committed_samples": int(sample_counts[committed]),
                "rna_late_samples": int(sample_counts[late]),
                **metrics,
                "rna_late_bootstrap_ci_low": ci_low,
                "rna_late_bootstrap_ci_high": ci_high,
                "rna_delay_pass": bool(passes),
            }
        )
    return pd.DataFrame(records), pseudobulk


def chromvar_context(
    chromvar_path: Path,
    metadata: pd.DataFrame,
    settings: dict,
) -> pd.DataFrame:
    chromvar = pd.read_csv(chromvar_path, index_col=0)
    motif_rows = [index for index in chromvar.index if str(index).startswith("GATA1-")]
    shared_cells = metadata.index.intersection(chromvar.columns)
    values = chromvar.loc[motif_rows, shared_cells].mean(axis=0)
    frame = metadata.loc[shared_cells, ["predicted.id", "batch", "origin"]].copy()
    frame["GATA1_chromvar_mean"] = values.loc[shared_cells].to_numpy(dtype=float)
    return frame.reset_index(names="cell_id")


def summarize_results(
    peak_evidence: pd.DataFrame, rna_evidence: pd.DataFrame
) -> tuple[pd.DataFrame, dict]:
    combined = peak_evidence.merge(rna_evidence, on="target", how="left")
    combined["external_peak_pass"] = (
        combined["atac_establishment_pass"].fillna(False)
        & combined["rna_delay_pass"].fillna(False)
    )
    edge = (
        combined.groupby(["TF", "target", "external_set"], as_index=False)
        .agg(
            frozen_peaks=("peak_id", "nunique"),
            mapped_peaks=("mapped_external_peaks", lambda value: int((value > 0).sum())),
            establishing_peaks=("atac_establishment_pass", "sum"),
            external_passing_peaks=("external_peak_pass", "sum"),
            rna_delay_pass=("rna_delay_pass", "first"),
        )
    )
    edge["external_edge_pass"] = edge.external_passing_peaks.gt(0)
    positive = combined.external_set.eq("targeted_E2")
    comparison = combined.external_set.eq("linked_non_E2")
    table = np.array(
        [
            [
                int((positive & combined.atac_establishment_pass).sum()),
                int((positive & ~combined.atac_establishment_pass).sum()),
            ],
            [
                int((comparison & combined.atac_establishment_pass).sum()),
                int((comparison & ~combined.atac_establishment_pass).sum()),
            ],
        ]
    )
    odds_ratio, p_value = fisher_exact(table) if table.sum() else (np.nan, np.nan)
    summary = {
        "frozen_E2_peaks": int(positive.sum()),
        "frozen_linked_non_E2_peaks": int(comparison.sum()),
        "mapped_E2_peaks": int((positive & combined.mapped_external_peaks.gt(0)).sum()),
        "mapped_linked_non_E2_peaks": int(
            (comparison & combined.mapped_external_peaks.gt(0)).sum()
        ),
        "establishing_E2_peaks": int(
            (positive & combined.atac_establishment_pass).sum()
        ),
        "establishing_linked_non_E2_peaks": int(
            (comparison & combined.atac_establishment_pass).sum()
        ),
        "external_edges_passed": int(
            edge.loc[edge.external_set.eq("targeted_E2"), "external_edge_pass"].sum()
        ),
        "E2_vs_non_E2_establishment_odds_ratio": float(odds_ratio),
        "E2_vs_non_E2_establishment_fisher_p_value": float(p_value),
        "interpretation": "external_temporal_validation_not_external_causal_validation",
    }
    return edge, summary


def run(config_path: Path, evidence_path: Path, output_dir: Path) -> None:
    config = load_config(config_path)
    settings = config["external_validation"]
    random_seed = int(config["project"]["random_seed"])
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(evidence_path)
    candidates = frozen_peak_sets(evidence, str(settings["primary_tf"]))
    metadata = pd.read_csv(Path(settings["atac_metadata"]), index_col=0)
    mapping, matrix = extract_overlapping_peak_matrix(
        Path(settings["atac_peak_matrix"]),
        candidates,
        int(settings["coordinate_mapping"]["nearest_peak_max_distance_bp"]),
    )
    peak_results = atac_peak_evidence(
        candidates, mapping, matrix, metadata, settings, random_seed
    )
    rna_results, rna_pseudobulk = extract_rna_evidence(
        Path(settings["rna_object"]),
        sorted(candidates.target.unique().tolist()),
        settings,
        random_seed,
    )
    chromvar = chromvar_context(Path(settings["atac_chromvar"]), metadata, settings)
    edges, summary = summarize_results(peak_results, rna_results)
    chromvar_states = chromvar["predicted.id"].astype(str).to_numpy()
    chromvar_values = chromvar.GATA1_chromvar_mean.to_numpy(dtype=float)
    chromvar_batches = chromvar.batch.astype(str).to_numpy()
    early = str(settings["atac_states"]["early"])
    committed = str(settings["atac_states"]["committed"])
    chromvar_effect = standardized_state_effect(
        chromvar_values, chromvar_states, early, committed
    )
    chromvar_ci = stratified_bootstrap_state_effect(
        chromvar_values,
        chromvar_states,
        chromvar_batches,
        early,
        committed,
        int(settings["bootstrap_iterations"]),
        float(settings["bootstrap_confidence"]),
        np.random.default_rng(random_seed),
    )
    summary.update(
        {
            "GATA1_chromvar_early_mean": float(
                np.mean(chromvar_values[chromvar_states == early])
            ),
            "GATA1_chromvar_committed_mean": float(
                np.mean(chromvar_values[chromvar_states == committed])
            ),
            "GATA1_chromvar_standardized_effect": chromvar_effect,
            "GATA1_chromvar_bootstrap_ci_low": chromvar_ci[0],
            "GATA1_chromvar_bootstrap_ci_high": chromvar_ci[1],
        }
    )
    mapping.to_csv(output_dir / "external_peak_mapping.csv", index=False)
    matrix.to_csv(output_dir / "targeted_external_atac_matrix.csv.gz")
    peak_results.to_csv(output_dir / "external_peak_evidence.csv", index=False)
    rna_results.to_csv(output_dir / "external_rna_evidence.csv", index=False)
    rna_pseudobulk.to_csv(output_dir / "external_rna_pseudobulk.csv", index=False)
    chromvar.to_csv(output_dir / "external_GATA1_chromvar.csv.gz", index=False)
    edges.to_csv(output_dir / "external_edge_summary.csv", index=False)
    with (output_dir / "external_validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--evidence",
        type=Path,
        default=Path(
            "reports/chromatin_mechanism/candidate_peak_final_E2_evidence.csv.gz"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("reports/external_validation"),
    )
    args = parser.parse_args()
    run(args.config, args.evidence, args.output_dir)


if __name__ == "__main__":
    main()
