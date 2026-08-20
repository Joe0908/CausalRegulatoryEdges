"""Stage-resolved external validation in adult human erythropoiesis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact

from edge_causality.external_validation import frozen_peak_sets, interval_overlap


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def liftover_candidates(candidates: pd.DataFrame, chain_path: Path) -> pd.DataFrame:
    try:
        from pyliftover import LiftOver
    except ImportError as error:  # pragma: no cover
        raise RuntimeError("Terminal validation requires the 'pyliftover' package") from error
    converter = LiftOver(str(chain_path))
    records = []
    for candidate in candidates.to_dict("records"):
        start_hits = converter.convert_coordinate(
            str(candidate["chromosome"]), int(candidate["start"])
        )
        end_hits = converter.convert_coordinate(
            str(candidate["chromosome"]), int(candidate["end"]) - 1
        )
        compatible = [
            (start, end)
            for start in start_hits
            for end in end_hits
            if start[0] == end[0] and start[2] == end[2]
        ]
        if not compatible:
            records.append(
                {
                    **candidate,
                    "hg19_chromosome": None,
                    "hg19_start": np.nan,
                    "hg19_end": np.nan,
                    "liftover_pass": False,
                }
            )
            continue
        start_hit, end_hit = compatible[0]
        coordinates = [int(start_hit[1]), int(end_hit[1])]
        records.append(
            {
                **candidate,
                "hg19_chromosome": start_hit[0],
                "hg19_start": min(coordinates),
                "hg19_end": max(coordinates) + 1,
                "liftover_pass": True,
            }
        )
    return pd.DataFrame(records)


def map_ludwig_peaks(
    lifted: pd.DataFrame, peaks: pd.DataFrame
) -> pd.DataFrame:
    records = []
    for candidate in lifted.loc[lifted.liftover_pass].itertuples(index=False):
        chromosome_peaks = peaks.loc[peaks.chromosome.eq(candidate.hg19_chromosome)]
        for peak in chromosome_peaks.itertuples(index=True):
            overlap = interval_overlap(
                int(candidate.hg19_start),
                int(candidate.hg19_end),
                int(peak.start),
                int(peak.end),
            )
            if overlap:
                records.append(
                    {
                        "peak_id": candidate.peak_id,
                        "ludwig_peak_index": int(peak.Index),
                        "ludwig_peak_id": f"{peak.chromosome}:{peak.start}-{peak.end}",
                        "overlap_bp": overlap,
                    }
                )
    return pd.DataFrame(records)


def parse_library_name(name: str) -> dict:
    fields = str(name).split("_")
    population = next(field for field in fields if field.startswith("P"))
    donor = next(field for field in fields if field.startswith("Donor"))
    replicate = next(field for field in fields if field.startswith("Rep"))
    return {
        "library": str(name),
        "donor": donor,
        "population": population,
        "replicate": replicate,
    }


def counts_to_log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    totals = counts.sum(axis=0).to_numpy(dtype=float)
    return np.log2(counts.div(np.maximum(totals, 1), axis=1) * 1_000_000 + 0.5)


def activation_stage(
    stage_means: np.ndarray,
    activation_fraction: float,
    minimum_dynamic_range: float,
) -> tuple[float, float]:
    baseline = float(stage_means[0])
    maximum = float(np.max(stage_means))
    dynamic_range = maximum - baseline
    maximum_stage = int(np.argmax(stage_means))
    if dynamic_range < minimum_dynamic_range or maximum_stage == 0:
        return float("nan"), dynamic_range
    threshold = baseline + activation_fraction * dynamic_range
    crossing = np.flatnonzero(stage_means >= threshold)
    crossing = crossing[crossing > 0]
    if len(crossing) == 0:
        return float("nan"), dynamic_range
    return float(crossing[0]), dynamic_range


def stage_means(
    values: np.ndarray, populations: np.ndarray, ordered: list[str]
) -> np.ndarray:
    return np.asarray(
        [np.mean(values[populations == population]) for population in ordered],
        dtype=float,
    )


def bootstrap_lead_fraction(
    atac_values: np.ndarray,
    atac_populations: np.ndarray,
    rna_values: np.ndarray,
    rna_populations: np.ndarray,
    ordered: list[str],
    activation_fraction: float,
    minimum_atac_dynamic_range: float,
    minimum_rna_dynamic_range: float,
    minimum_lead_stages: int,
    iterations: int,
    rng: np.random.Generator,
) -> float:
    successes = 0
    atac_groups = [
        np.flatnonzero(atac_populations == population) for population in ordered
    ]
    rna_groups = [
        np.flatnonzero(rna_populations == population) for population in ordered
    ]
    for _ in range(iterations):
        atac_sample = np.concatenate(
            [rng.choice(group, len(group), replace=True) for group in atac_groups]
        )
        rna_sample = np.concatenate(
            [rng.choice(group, len(group), replace=True) for group in rna_groups]
        )
        atac_onset, _ = activation_stage(
            stage_means(
                atac_values[atac_sample], atac_populations[atac_sample], ordered
            ),
            activation_fraction,
            minimum_atac_dynamic_range,
        )
        rna_onset, _ = activation_stage(
            stage_means(rna_values[rna_sample], rna_populations[rna_sample], ordered),
            activation_fraction,
            minimum_rna_dynamic_range,
        )
        if np.isfinite(atac_onset) and np.isfinite(rna_onset):
            successes += int(rna_onset - atac_onset >= minimum_lead_stages)
    # Draws that fail either pre-registered dynamic-range requirement are
    # failures, not missing observations.  Using all iterations as the
    # denominator prevents unstable low-signal features from acquiring inflated
    # conditional support.
    return float(successes / iterations) if iterations else float("nan")


def build_trajectories(
    candidates: pd.DataFrame,
    mapping: pd.DataFrame,
    atac_counts: pd.DataFrame,
    rna_counts: pd.DataFrame,
    ordered_populations: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    atac_metadata = pd.DataFrame(
        [parse_library_name(column) for column in atac_counts.columns]
    ).set_index("library")
    rna_numeric = rna_counts.drop(columns="genes").copy()
    rna_metadata = pd.DataFrame(
        [parse_library_name(column) for column in rna_numeric.columns]
    ).set_index("library")
    atac_log_cpm = counts_to_log2_cpm(atac_counts)
    rna_numeric.index = rna_counts.genes.astype(str)
    rna_numeric = rna_numeric.groupby(level=0).sum()
    rna_log_cpm = counts_to_log2_cpm(rna_numeric)

    atac_rows = []
    for candidate in candidates.itertuples(index=False):
        indices = mapping.loc[
            mapping.peak_id.eq(candidate.peak_id), "ludwig_peak_index"
        ].to_numpy(dtype=int)
        if len(indices) == 0:
            continue
        raw_profile = atac_counts.iloc[indices].sum(axis=0).to_frame().T
        profile = counts_to_log2_cpm_with_totals(raw_profile, atac_counts.sum(axis=0))
        for library, value in profile.iloc[0].items():
            atac_rows.append(
                {
                    "peak_id": candidate.peak_id,
                    "target": candidate.target,
                    "external_set": candidate.external_set,
                    "library": library,
                    "population": atac_metadata.loc[library, "population"],
                    "donor": atac_metadata.loc[library, "donor"],
                    "log2_cpm": float(value),
                }
            )
    atac_trajectory = pd.DataFrame(atac_rows)
    rna_rows = []
    for target in sorted(candidates.target.unique()):
        if target not in rna_log_cpm.index:
            continue
        for library, value in rna_log_cpm.loc[target].items():
            rna_rows.append(
                {
                    "target": target,
                    "library": library,
                    "population": rna_metadata.loc[library, "population"],
                    "donor": rna_metadata.loc[library, "donor"],
                    "log2_cpm": float(value),
                }
            )
    return atac_trajectory, pd.DataFrame(rna_rows), atac_metadata, rna_metadata


def counts_to_log2_cpm_with_totals(
    counts: pd.DataFrame, totals: pd.Series
) -> pd.DataFrame:
    return np.log2(counts.div(np.maximum(totals.to_numpy(dtype=float), 1), axis=1) * 1_000_000 + 0.5)


def evaluate_peaks(
    candidates: pd.DataFrame,
    mapping: pd.DataFrame,
    atac_trajectory: pd.DataFrame,
    rna_trajectory: pd.DataFrame,
    settings: dict,
    random_seed: int,
) -> pd.DataFrame:
    ordered = list(settings["ordered_populations"])
    rng = np.random.default_rng(random_seed)
    records = []
    for candidate in candidates.itertuples(index=False):
        atac = atac_trajectory.loc[
            atac_trajectory.peak_id.eq(candidate.peak_id)
        ]
        rna = rna_trajectory.loc[rna_trajectory.target.eq(candidate.target)]
        mapped = mapping.loc[mapping.peak_id.eq(candidate.peak_id)]
        base = {
            **candidate._asdict(),
            "mapped_ludwig_peaks": int(mapped.ludwig_peak_index.nunique())
            if not mapped.empty
            else 0,
            "ludwig_peak_ids": ";".join(mapped.ludwig_peak_id.astype(str))
            if not mapped.empty
            else "",
        }
        if atac.empty or rna.empty:
            records.append(
                {
                    **base,
                    "atac_activation_stage_index": np.nan,
                    "rna_activation_stage_index": np.nan,
                    "bootstrap_lead_fraction": np.nan,
                    "terminal_temporal_pass": False,
                }
            )
            continue
        atac_values = atac.log2_cpm.to_numpy(dtype=float)
        atac_populations = atac.population.astype(str).to_numpy()
        rna_values = rna.log2_cpm.to_numpy(dtype=float)
        rna_populations = rna.population.astype(str).to_numpy()
        atac_curve = stage_means(atac_values, atac_populations, ordered)
        rna_curve = stage_means(rna_values, rna_populations, ordered)
        atac_onset, atac_dynamic = activation_stage(
            atac_curve,
            float(settings["activation_fraction"]),
            float(settings["minimum_atac_dynamic_range_log2_cpm"]),
        )
        rna_onset, rna_dynamic = activation_stage(
            rna_curve,
            float(settings["activation_fraction"]),
            float(settings["minimum_rna_dynamic_range_log2_cpm"]),
        )
        support = bootstrap_lead_fraction(
            atac_values,
            atac_populations,
            rna_values,
            rna_populations,
            ordered,
            float(settings["activation_fraction"]),
            float(settings["minimum_atac_dynamic_range_log2_cpm"]),
            float(settings["minimum_rna_dynamic_range_log2_cpm"]),
            int(settings["minimum_atac_lead_stages"]),
            int(settings["bootstrap_iterations"]),
            rng,
        )
        lead = rna_onset - atac_onset
        passes = (
            np.isfinite(lead)
            and lead >= int(settings["minimum_atac_lead_stages"])
            and support >= float(settings["minimum_bootstrap_lead_fraction"])
        )
        record = {
            **base,
            "atac_activation_stage_index": atac_onset,
            "rna_activation_stage_index": rna_onset,
            "atac_activation_population": ordered[int(atac_onset)]
            if np.isfinite(atac_onset)
            else "",
            "rna_activation_population": ordered[int(rna_onset)]
            if np.isfinite(rna_onset)
            else "",
            "atac_dynamic_range_log2_cpm": atac_dynamic,
            "rna_dynamic_range_log2_cpm": rna_dynamic,
            "atac_lead_stages": lead,
            "bootstrap_lead_fraction": support,
            "terminal_temporal_pass": bool(passes),
        }
        for index, population in enumerate(ordered):
            record[f"atac_mean_{population}"] = float(atac_curve[index])
            record[f"rna_mean_{population}"] = float(rna_curve[index])
        records.append(record)
    return pd.DataFrame(records)


def summarize(evidence: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    edges = (
        evidence.groupby(["TF", "target", "external_set"], as_index=False)
        .agg(
            frozen_peaks=("peak_id", "nunique"),
            mapped_peaks=("mapped_ludwig_peaks", lambda value: int((value > 0).sum())),
            temporal_passing_peaks=("terminal_temporal_pass", "sum"),
        )
    )
    edges["terminal_external_edge_pass"] = edges.temporal_passing_peaks.gt(0)
    positive = evidence.external_set.eq("chapter3_E2")
    comparison = evidence.external_set.eq("linked_non_E2")
    table = np.array(
        [
            [
                int((positive & evidence.terminal_temporal_pass).sum()),
                int((positive & ~evidence.terminal_temporal_pass).sum()),
            ],
            [
                int((comparison & evidence.terminal_temporal_pass).sum()),
                int((comparison & ~evidence.terminal_temporal_pass).sum()),
            ],
        ]
    )
    odds, p_value = fisher_exact(table)
    positive_edges = edges.external_set.eq("chapter3_E2")
    comparison_edges = edges.external_set.eq("linked_non_E2")
    edge_table = np.array(
        [
            [
                int(
                    (positive_edges & edges.terminal_external_edge_pass).sum()
                ),
                int(
                    (positive_edges & ~edges.terminal_external_edge_pass).sum()
                ),
            ],
            [
                int(
                    (comparison_edges & edges.terminal_external_edge_pass).sum()
                ),
                int(
                    (comparison_edges & ~edges.terminal_external_edge_pass).sum()
                ),
            ],
        ]
    )
    edge_odds, edge_p_value = fisher_exact(edge_table)
    summary = {
        "frozen_E2_peaks": int(positive.sum()),
        "mapped_E2_peaks": int((positive & evidence.mapped_ludwig_peaks.gt(0)).sum()),
        "temporal_passing_E2_peaks": int(
            (positive & evidence.terminal_temporal_pass).sum()
        ),
        "frozen_linked_non_E2_peaks": int(comparison.sum()),
        "mapped_linked_non_E2_peaks": int(
            (comparison & evidence.mapped_ludwig_peaks.gt(0)).sum()
        ),
        "temporal_passing_linked_non_E2_peaks": int(
            (comparison & evidence.terminal_temporal_pass).sum()
        ),
        "terminal_external_edges_passed": int(
            edges.loc[
                edges.external_set.eq("chapter3_E2"),
                "terminal_external_edge_pass",
            ].sum()
        ),
        "E2_vs_non_E2_temporal_odds_ratio": float(odds),
        "E2_vs_non_E2_temporal_fisher_p_value": float(p_value),
        "E2_edges_tested": int(positive_edges.sum()),
        "linked_non_E2_edges_tested": int(comparison_edges.sum()),
        "linked_non_E2_edges_passed": int(
            (comparison_edges & edges.terminal_external_edge_pass).sum()
        ),
        "E2_vs_non_E2_edge_odds_ratio": float(edge_odds),
        "E2_vs_non_E2_edge_fisher_p_value": float(edge_p_value),
        "interpretation": "stage_resolved_external_temporal_validation_not_causal_retest",
    }
    return edges, summary


def run(config_path: Path, chapter3_path: Path, output_dir: Path) -> None:
    config = load_config(config_path)
    parent = config["external_validation"]
    settings = parent["terminal_erythroid_fallback"]
    candidates = frozen_peak_sets(
        pd.read_csv(chapter3_path), str(parent["primary_tf"])
    )
    lifted = liftover_candidates(candidates, Path(settings["liftover_chain"]))
    peaks = pd.read_csv(
        settings["atac_peaks"],
        sep="\t",
        header=None,
        names=["chromosome", "start", "end"],
    )
    mapping = map_ludwig_peaks(lifted, peaks)
    atac_counts = pd.read_csv(settings["atac_counts"], sep="\t")
    rna_counts = pd.read_csv(settings["rna_counts"], sep="\t")
    atac_trajectory, rna_trajectory, _, _ = build_trajectories(
        candidates,
        mapping,
        atac_counts,
        rna_counts,
        list(settings["ordered_populations"]),
    )
    evidence = evaluate_peaks(
        candidates,
        mapping,
        atac_trajectory,
        rna_trajectory,
        settings,
        int(config["project"]["random_seed"]),
    )
    edges, summary = summarize(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    lifted.to_csv(output_dir / "ludwig_lifted_candidates.csv", index=False)
    mapping.to_csv(output_dir / "ludwig_peak_mapping.csv", index=False)
    atac_trajectory.to_csv(output_dir / "ludwig_atac_trajectory.csv", index=False)
    rna_trajectory.to_csv(output_dir / "ludwig_rna_trajectory.csv", index=False)
    evidence.to_csv(output_dir / "ludwig_peak_temporal_evidence.csv", index=False)
    edges.to_csv(output_dir / "ludwig_edge_summary.csv", index=False)
    with (output_dir / "ludwig_validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--chapter3",
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
    run(args.config, args.chapter3, args.output_dir)


if __name__ == "__main__":
    main()
