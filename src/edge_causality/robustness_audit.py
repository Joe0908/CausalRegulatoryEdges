"""Post-freeze threshold and inferential-unit audit.

This module only reclassifies already-computed evidence. It never refits the
frozen E0, E1 or E2 models and therefore cannot silently change the primary
analysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu
import yaml


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def edge_key(frame: pd.DataFrame) -> pd.Series:
    return frame["TF"].astype(str) + "->" + frame["target"].astype(str)


def e0_stability_grid(edges: pd.DataFrame, grid: dict) -> pd.DataFrame:
    frozen = edges["bootstrap_top5_frequency"].ge(0.70)
    frozen_keys = set(edge_key(edges.loc[frozen]))
    rows = []
    for percent in grid["selection_percent"]:
        column = f"bootstrap_top{int(percent)}_frequency"
        for frequency in grid["minimum_bootstrap_frequency"]:
            selected = edges[column].ge(float(frequency))
            keys = set(edge_key(edges.loc[selected]))
            union = frozen_keys | keys
            row = {
                "selection_percent": int(percent),
                "minimum_bootstrap_frequency": float(frequency),
                "selected_edges": int(selected.sum()),
                "selected_TFs": int(edges.loc[selected, "TF"].nunique()),
                "overlap_with_frozen_edges": int(len(keys & frozen_keys)),
                "jaccard_with_frozen": float(len(keys & frozen_keys) / len(union))
                if union
                else np.nan,
                "median_absolute_association": float(
                    edges.loc[selected, "absolute_association"].median()
                )
                if selected.any()
                else np.nan,
            }
            for truth in ("author_effect_025", "author_supported_concordant"):
                if truth in edges:
                    row[f"{truth}_rate"] = float(edges.loc[selected, truth].mean())
                    row[f"{truth}_edges"] = int(edges.loc[selected, truth].sum())
            rows.append(row)
    return pd.DataFrame(rows)


def call_e1(
    edges: pd.DataFrame,
    fdr_max: float = 0.05,
    minimum_absolute_log2_fold_change: float = 0.25,
    require_guide_direction_consistency: bool = True,
    require_leave_one_guide_out: bool = True,
) -> pd.Series:
    selected = (
        edges["effective_guides_used"].ge(2)
        & edges["perturbation_fdr"].lt(fdr_max)
        & edges["perturbation_log2fc"].abs().ge(minimum_absolute_log2_fold_change)
    )
    if require_guide_direction_consistency:
        selected &= edges["guide_direction_consistent"].fillna(False)
    if require_leave_one_guide_out:
        selected &= edges["leave_one_guide_out_direction_consistent"].fillna(False)
    return selected


def e1_sensitivity(edges: pd.DataFrame, grid: dict) -> pd.DataFrame:
    rows = []

    def append(parameter: str, value: float | str, selected: pd.Series) -> None:
        rows.append(
            {
                "parameter": parameter,
                "value": value,
                "supported_edges": int(selected.sum()),
                "direction_concordant_edges": int(
                    (selected & edges["direction_concordant_with_knockout"]).sum()
                ),
                "TFs_with_support": int(edges.loc[selected, "TF"].nunique()),
            }
        )

    for value in grid["fdr_max"]:
        append("fdr_max", float(value), call_e1(edges, fdr_max=float(value)))
    for value in grid["minimum_absolute_log2_fold_change"]:
        append(
            "minimum_absolute_log2_fold_change",
            float(value),
            call_e1(edges, minimum_absolute_log2_fold_change=float(value)),
        )
    append("guide_robustness", "pooled_effect_only", call_e1(edges, require_guide_direction_consistency=False, require_leave_one_guide_out=False))
    append("guide_robustness", "guide_consistency", call_e1(edges, require_leave_one_guide_out=False))
    append("guide_robustness", "guide_and_leave_one_out", call_e1(edges))
    return pd.DataFrame(rows)


FROZEN_E2 = {
    "link_fdr_max": 0.05,
    "minimum_absolute_link_correlation": 0.03,
    "minimum_bootstrap_sign_fraction": 0.80,
    "minimum_peak_libraries": 4,
    "atac_fdr_max": 0.05,
    "minimum_absolute_atac_effect": 0.20,
    "minimum_consistent_guides": 2,
    "minimum_absolute_rna_effect": 0.25,
    "motif_relative_score_threshold": 0.85,
}


def call_e2(evidence: pd.DataFrame, **overrides: float) -> pd.Series:
    settings = FROZEN_E2 | overrides
    linked = (
        evidence["link_fdr"].lt(settings["link_fdr_max"])
        & evidence["link_correlation"].abs().ge(
            settings["minimum_absolute_link_correlation"]
        )
        & evidence["link_bootstrap_sign_fraction"].ge(
            settings["minimum_bootstrap_sign_fraction"]
        )
        & evidence["libraries_present"].ge(settings["minimum_peak_libraries"])
    )
    atac = (
        evidence["strongest_effect_fdr"].lt(settings["atac_fdr_max"])
        & evidence["strongest_effect"].abs().ge(
            settings["minimum_absolute_atac_effect"]
        )
        & evidence["strongest_consistent_guides"].ge(
            settings["minimum_consistent_guides"]
        )
    )
    predicted_rna_sign = np.sign(
        evidence["link_correlation"] * evidence["strongest_effect"]
    )
    rna = (
        np.sign(evidence["rna_effect_at_strongest_atac_timepoint"])
        == predicted_rna_sign
    ) & evidence["rna_effect_at_strongest_atac_timepoint"].abs().ge(
        settings["minimum_absolute_rna_effect"]
    )
    motif = evidence["motif_best_relative_score"].ge(
        settings["motif_relative_score_threshold"]
    )
    return linked & atac & rna & motif


def e2_sensitivity(evidence: pd.DataFrame, grid: dict) -> pd.DataFrame:
    frozen = call_e2(evidence)
    frozen_keys = set(evidence.loc[frozen, "peak_id"].astype(str))
    rows = []
    for parameter, values in grid.items():
        for value in values:
            selected = call_e2(evidence, **{parameter: float(value)})
            keys = set(evidence.loc[selected, "peak_id"].astype(str))
            union = keys | frozen_keys
            rows.append(
                {
                    "parameter": parameter,
                    "value": float(value),
                    "E2_peaks": int(selected.sum()),
                    "E2_edges": int(edge_key(evidence.loc[selected]).nunique()),
                    "overlap_with_frozen_peaks": int(len(keys & frozen_keys)),
                    "jaccard_with_frozen": float(len(keys & frozen_keys) / len(union))
                    if union
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _direction_stats(values: pd.Series) -> dict:
    values = values.dropna().astype(float)
    negative = int(values.lt(0).sum())
    return {
        "n": int(len(values)),
        "negative": negative,
        "negative_fraction": float(negative / len(values)) if len(values) else np.nan,
        "median_correlation": float(values.median()) if len(values) else np.nan,
        "one_sided_sign_p": float(
            binomtest(negative, len(values), 0.5, alternative="greater").pvalue
        )
        if len(values)
        else np.nan,
    }


def chromatin_inferential_unit_audit(evidence: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    frozen = call_e2(evidence)
    linked = evidence["link_pass"].astype(bool)
    peak_groups = {
        "E2": evidence.loc[frozen],
        "linked_non_E2": evidence.loc[linked & ~frozen],
    }
    rows = []
    for group, frame in peak_groups.items():
        for outcome, column in (
            ("ATAC", "baseline_access_vs_absolute_atac_effect_correlation"),
            ("RNA", "baseline_access_vs_absolute_rna_effect_correlation"),
        ):
            rows.append({"unit": "peak", "group": group, "outcome": outcome, **_direction_stats(frame[column])})

    edge_rows = []
    for _, frame in evidence.loc[linked].groupby(["TF", "target"], observed=True):
        selected = frame.loc[frozen.loc[frame.index]] if frozen.loc[frame.index].any() else frame
        edge_rows.append(
            {
                "TF": frame["TF"].iloc[0],
                "target": frame["target"].iloc[0],
                "group": "E2" if frozen.loc[frame.index].any() else "linked_non_E2",
                "ATAC": selected["baseline_access_vs_absolute_atac_effect_correlation"].median(),
                "RNA": selected["baseline_access_vs_absolute_rna_effect_correlation"].median(),
            }
        )
    edge_frame = pd.DataFrame(edge_rows)
    for group, frame in edge_frame.groupby("group", observed=True):
        for outcome in ("ATAC", "RNA"):
            rows.append({"unit": "edge", "group": group, "outcome": outcome, **_direction_stats(frame[outcome])})

    comparisons = {}
    for unit, frames in (("peak", peak_groups), ("edge", {g: x for g, x in edge_frame.groupby("group", observed=True)})):
        for outcome, column in (
            ("ATAC", "baseline_access_vs_absolute_atac_effect_correlation"),
            ("RNA", "baseline_access_vs_absolute_rna_effect_correlation"),
        ):
            if unit == "peak":
                positive = frames["E2"][column].dropna()
                comparison = frames["linked_non_E2"][column].dropna()
            else:
                positive = frames["E2"][outcome].dropna()
                comparison = frames["linked_non_E2"][outcome].dropna()
            comparisons[f"{unit}_{outcome}_E2_less_than_linked_non_E2_mannwhitney_p"] = float(
                mannwhitneyu(positive, comparison, alternative="less").pvalue
            )
    return pd.DataFrame(rows), comparisons


def run(
    grid_path: Path,
    observational_path: Path,
    e1_path: Path,
    e2_path: Path,
    output_dir: Path,
) -> dict:
    grid = load_yaml(grid_path)
    observational = pd.read_csv(observational_path)
    e1_edges = pd.read_csv(e1_path)
    evidence = pd.read_csv(e2_path)
    e0 = e0_stability_grid(observational, grid["e0"])
    e1 = e1_sensitivity(e1_edges, grid["e1"])
    e2 = e2_sensitivity(evidence, grid["e2"])
    units, comparisons = chromatin_inferential_unit_audit(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    e0.to_csv(output_dir / "e0_stability_grid.csv", index=False)
    e1.to_csv(output_dir / "e1_threshold_sensitivity.csv", index=False)
    e2.to_csv(output_dir / "e2_threshold_sensitivity.csv", index=False)
    units.to_csv(output_dir / "chromatin_inferential_unit_audit.csv", index=False)
    summary = {
        "frozen_E0_edges": int(
            e0.loc[e0.selection_percent.eq(5) & e0.minimum_bootstrap_frequency.eq(0.70), "selected_edges"].iloc[0]
        ),
        "frozen_E1_edges": int(call_e1(e1_edges).sum()),
        "frozen_E2_peaks": int(call_e2(evidence).sum()),
        "frozen_E2_edges": int(edge_key(evidence.loc[call_e2(evidence)]).nunique()),
        "chromatin_comparisons": comparisons,
        "interpretation": "post-freeze sensitivity; primary estimates unchanged",
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=Path, default=Path("config/robustness.yaml"))
    parser.add_argument(
        "--observational",
        type=Path,
        default=Path("reports/author_truth/observational_edges_with_author_truth.csv.gz"),
    )
    parser.add_argument(
        "--e1", type=Path, default=Path("reports/validation/E0_to_E1_edge_matrix.csv.gz")
    )
    parser.add_argument(
        "--e2",
        type=Path,
        default=Path("reports/chromatin_mechanism/candidate_peak_final_E2_evidence.csv.gz"),
    )
    parser.add_argument("--output", type=Path, default=Path("reports/robustness"))
    args = parser.parse_args()
    print(json.dumps(run(args.grid, args.observational, args.e1, args.e2, args.output), indent=2))


if __name__ == "__main__":
    main()
