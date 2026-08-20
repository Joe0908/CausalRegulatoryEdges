"""Targeted chromatin-mechanism tests for frozen TF-to-gene candidates.

The analysis deliberately separates three questions:

1. Is a local ACR associated with its gene in non-targeting control cells?
2. Does perturbing the candidate TF change accessibility of that ACR?
3. Is the accessibility change directionally compatible with the RNA effect?

Motif evidence is added in a separate, sequence-based step.  Results from this
module are therefore labelled ``provisional_E2`` rather than final E2.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr
import yaml

from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import (
    fit_interaction_model,
    interaction_design,
    per_guide_context_effects,
)


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def residualize(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    return y - design @ (np.linalg.pinv(design) @ y)


def safe_correlation(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def covariate_design(metadata: pd.DataFrame) -> np.ndarray:
    """Depth plus library/state effects; all are pre-perturbation or controls-only."""
    log_rna = np.log1p(metadata.nCount_RNA.to_numpy(dtype=float))
    log_atac = np.log1p(metadata.nCount_atac.to_numpy(dtype=float))
    continuous = np.column_stack(
        [
            (log_rna - log_rna.mean()) / max(log_rna.std(), 1e-8),
            (log_atac - log_atac.mean()) / max(log_atac.std(), 1e-8),
        ]
    )
    categories = pd.get_dummies(
        metadata[["replicate", "new_CellType"]].astype(str), drop_first=True
    ).to_numpy(dtype=float)
    return np.column_stack([np.ones(len(metadata)), continuous, categories])


def bootstrap_sign_fraction(
    x: np.ndarray,
    y: np.ndarray,
    strata: np.ndarray,
    observed: float,
    iterations: int,
    rng: np.random.Generator,
) -> float:
    signs = []
    stratum_indices = [np.flatnonzero(strata == value) for value in np.unique(strata)]
    for _ in range(iterations):
        sampled = np.concatenate(
            [rng.choice(index, size=len(index), replace=True) for index in stratum_indices]
        )
        if np.std(x[sampled]) == 0 or np.std(y[sampled]) == 0:
            continue
        signs.append(np.sign(np.corrcoef(x[sampled], y[sampled])[0, 1]))
    if not signs or observed == 0:
        return 0.0
    return float(np.mean(np.asarray(signs) == np.sign(observed)))


def test_control_links(
    genes: sparse.csr_matrix,
    atac: sparse.csr_matrix,
    metadata: pd.DataFrame,
    candidates: pd.DataFrame,
    peaks: pd.DataFrame,
    settings: dict,
    random_seed: int,
) -> pd.DataFrame:
    controls = metadata.perturbation_name.isin(settings["link_control_guides"])
    rng = np.random.default_rng(random_seed)
    records = []
    for peak in peaks.itertuples(index=False):
        present_replicates = set(str(peak.replicates_present).split(";"))
        selected = controls & metadata.replicate.isin(present_replicates)
        cell_index = np.flatnonzero(selected.to_numpy())
        candidate = candidates.loc[candidates.target.eq(peak.target)].iloc[0]
        gene_index = int(candidate.gene_index)
        gene_count = genes[cell_index, gene_index].toarray().ravel()
        peak_count = atac[cell_index, int(peak.peak_index)].toarray().ravel()
        selected_metadata = metadata.iloc[cell_index]
        gene_value = np.log2(
            gene_count
            / np.maximum(selected_metadata.nCount_RNA.to_numpy(dtype=float), 1)
            * 1_000_000
            + 0.5
        )
        peak_value = np.log1p(
            peak_count
            / np.maximum(selected_metadata.nCount_atac.to_numpy(dtype=float), 1)
            * 10_000
        )
        design = covariate_design(selected_metadata)
        gene_residual = residualize(gene_value, design)
        peak_residual = residualize(peak_value, design)
        if np.std(gene_residual) == 0 or np.std(peak_residual) == 0:
            correlation, p_value = 0.0, 1.0
        else:
            correlation, p_value = pearsonr(peak_residual, gene_residual)
        sign_fraction = bootstrap_sign_fraction(
            peak_residual,
            gene_residual,
            selected_metadata.replicate.astype(str).to_numpy(),
            float(correlation),
            int(settings["bootstrap_iterations"]),
            rng,
        )
        records.append(
            {
                "peak_index": int(peak.peak_index),
                "link_cells": int(len(cell_index)),
                "gene_detection_fraction": float(np.mean(gene_count > 0)),
                "peak_detection_fraction": float(np.mean(peak_count > 0)),
                "link_correlation": float(correlation),
                "link_p_value": float(p_value),
                "link_bootstrap_sign_fraction": sign_fraction,
            }
        )
    output = pd.DataFrame(records)
    output["link_fdr"] = bh_adjust(output.link_p_value.to_numpy())
    library_eligible = (
        peaks.set_index("peak_index")
        .loc[output.peak_index, "libraries_present"]
        .to_numpy()
        .astype(int)
        >= int(settings["minimum_peak_libraries"])
    )
    output["link_pass"] = (
        output.link_fdr.lt(float(settings["link_fdr_max"]))
        & output.link_correlation.abs().ge(
            float(settings["minimum_absolute_link_correlation"])
        )
        & output.link_bootstrap_sign_fraction.ge(
            float(settings["minimum_bootstrap_sign_fraction"])
        )
        & library_eligible
    )
    return output


def make_pseudobulk(
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    selected: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    block_metadata = metadata.iloc[selected].copy()
    keys = pd.MultiIndex.from_frame(
        block_metadata[["replicate", "Timepoint", "perturbation_name", "target"]]
        .astype(str)
        .rename(columns={"perturbation_name": "guide"})
    )
    codes, unique = pd.factorize(keys, sort=True)
    assignment = sparse.csr_matrix(
        (np.ones(len(selected)), (codes, np.arange(len(selected)))),
        shape=(len(unique), len(selected)),
    )
    aggregated = (assignment @ counts[selected]).toarray().astype(float)
    groups = unique.to_frame(index=False)
    groups.columns = ["replicate", "timepoint", "guide", "target"]
    groups["n_cells"] = np.bincount(codes)
    groups["library_total"] = np.bincount(
        codes, weights=block_metadata.nCount_atac.to_numpy(dtype=float)
    )
    return aggregated, groups


def fit_feature_effects(
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    feature_indices: np.ndarray,
    feature_table: pd.DataFrame,
    tf: str,
    reference_guides: list[str],
    timepoints: list[str],
    minimum_cells: int,
    denominator_column: str,
) -> pd.DataFrame:
    target_guides = sorted(
        metadata.loc[metadata.target.eq(tf), "perturbation_name"].unique().tolist()
    )
    selected = np.flatnonzero(
        (
            metadata.target.eq(tf)
            | metadata.perturbation_name.isin(reference_guides)
        ).to_numpy()
    )
    pseudobulk, groups = make_pseudobulk(counts, metadata, selected)
    if denominator_column == "nCount_RNA":
        block_metadata = metadata.iloc[selected]
        keys = pd.MultiIndex.from_frame(
            block_metadata[
                ["replicate", "Timepoint", "perturbation_name", "target"]
            ]
            .astype(str)
            .rename(columns={"perturbation_name": "guide"})
        )
        codes, _ = pd.factorize(keys, sort=True)
        groups["library_total"] = np.bincount(
            codes, weights=block_metadata.nCount_RNA.to_numpy(dtype=float)
        )
    keep = groups.n_cells.ge(minimum_cells).to_numpy()
    groups = groups.loc[keep].reset_index(drop=True)
    log_cpm = np.log2(
        pseudobulk[keep]
        / np.maximum(groups.library_total.to_numpy()[:, None], 1)
        * 1_000_000
        + 0.5
    )
    reduced, full, contrasts, _ = interaction_design(
        groups, target_guides, timepoints
    )
    weights = groups.n_cells.to_numpy(dtype=float)
    weights /= np.median(weights)
    model = fit_interaction_model(
        log_cpm[:, feature_indices], reduced, full, weights, contrasts
    )
    output = feature_table.copy().reset_index(drop=True)
    output.insert(0, "TF_tested", tf)
    output["targeting_guides"] = ";".join(target_guides)
    for name, values in model.items():
        output[name] = values
    guide_effects = per_guide_context_effects(
        log_cpm,
        groups,
        feature_indices,
        target_guides,
        reference_guides,
        timepoints,
    )
    for timepoint in timepoints:
        pooled = output[f"effect_{timepoint}"].to_numpy()
        matrix = guide_effects[timepoint]
        output[f"consistent_guides_{timepoint}"] = np.sum(
            np.sign(matrix) == np.sign(pooled)[None, :], axis=0
        )
    return output


def add_fdr_and_sensitivity(
    effects: pd.DataFrame, settings: dict, timepoints: list[str]
) -> pd.DataFrame:
    effects = effects.copy()
    effects["interaction_fdr"] = bh_adjust(effects.interaction_p_value.to_numpy())
    for timepoint in timepoints:
        effects[f"effect_fdr_{timepoint}"] = bh_adjust(
            effects[f"effect_p_value_{timepoint}"].to_numpy()
        )
    effect_matrix = effects[[f"effect_{t}" for t in timepoints]].to_numpy()
    strongest = np.argmax(np.abs(effect_matrix), axis=1)
    effects["strongest_timepoint"] = [timepoints[i] for i in strongest]
    effects["strongest_effect"] = effect_matrix[np.arange(len(effects)), strongest]
    effects["strongest_effect_fdr"] = [
        effects.iloc[i][f"effect_fdr_{timepoints[j]}"]
        for i, j in enumerate(strongest)
    ]
    effects["strongest_consistent_guides"] = [
        effects.iloc[i][f"consistent_guides_{timepoints[j]}"]
        for i, j in enumerate(strongest)
    ]
    effects["perturbation_sensitive"] = (
        effects.strongest_effect_fdr.lt(0.05)
        & effects.strongest_effect.abs().ge(
            float(settings["minimum_absolute_atac_effect"])
        )
        & effects.strongest_consistent_guides.ge(
            int(settings["minimum_consistent_guides"])
        )
    )
    effects["state_dependent"] = effects.interaction_fdr.lt(
        float(settings["atac_interaction_fdr_max"])
    )
    return effects


def add_baseline_gating_diagnostics(
    evidence: pd.DataFrame,
    genes: sparse.csr_matrix,
    atac: sparse.csr_matrix,
    metadata: pd.DataFrame,
    candidates: pd.DataFrame,
    control_guides: list[str],
    timepoints: list[str],
) -> pd.DataFrame:
    """Compare control accessibility with effect magnitude across collection days.

    A simple permissive-gating model predicts a positive relationship: the TF
    perturbation should matter most where the linked peak is already accessible.
    A negative relationship instead suggests an establishment/priming phase in
    which perturbation prevents a regulatory element from becoming fully open.
    With only four collection days these correlations are descriptive, not a
    separate significance test.
    """
    output = evidence.copy()
    controls = metadata.perturbation_name.isin(control_guides).to_numpy()
    gene_lookup = candidates.set_index("target").gene_index.astype(int).to_dict()
    access_profiles = []
    expression_profiles = []
    access_vs_atac = []
    access_vs_rna = []
    expression_vs_rna = []
    for _, row in output.iterrows():
        access = []
        expression = []
        gene_index = gene_lookup[row["target"]]
        for timepoint in timepoints:
            selected = np.flatnonzero(
                controls & metadata.Timepoint.eq(timepoint).to_numpy()
            )
            peak_count = atac[selected, int(row["peak_index"])]
            gene_count = genes[selected, gene_index].toarray().ravel()
            access.append(float((peak_count > 0).mean()))
            expression.append(
                float(
                    np.mean(
                        np.log2(
                            gene_count
                            / np.maximum(
                                metadata.iloc[selected].nCount_RNA.to_numpy(dtype=float),
                                1,
                            )
                            * 1_000_000
                            + 0.5
                        )
                    )
                )
            )
        atac_effect = np.abs(
            np.array([row[f"effect_{timepoint}"] for timepoint in timepoints])
        )
        rna_effect = np.abs(
            np.array([row[f"rna_effect_{timepoint}"] for timepoint in timepoints])
        )
        access_array = np.asarray(access)
        expression_array = np.asarray(expression)
        access_profiles.append(access)
        expression_profiles.append(expression)
        access_vs_atac.append(safe_correlation(access_array, atac_effect))
        access_vs_rna.append(safe_correlation(access_array, rna_effect))
        expression_vs_rna.append(safe_correlation(expression_array, rna_effect))
    for index, timepoint in enumerate(timepoints):
        output[f"control_accessible_fraction_{timepoint}"] = [
            profile[index] for profile in access_profiles
        ]
        output[f"control_gene_expression_{timepoint}"] = [
            profile[index] for profile in expression_profiles
        ]
    output["baseline_access_vs_absolute_atac_effect_correlation"] = access_vs_atac
    output["baseline_access_vs_absolute_rna_effect_correlation"] = access_vs_rna
    output["baseline_expression_vs_absolute_rna_effect_correlation"] = (
        expression_vs_rna
    )
    return output


def run(
    config_path: Path,
    input_dir: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["chromatin_mechanism"]
    timepoints = list(config["state_dependence"]["ordered_timepoints"])
    references = list(config["data"]["intervention_reference_guides"])
    genes = sparse.load_npz(input_dir / "candidate_rna_counts.npz").tocsr()
    atac = sparse.load_npz(input_dir / "candidate_atac_counts.npz").tocsr()
    metadata = pd.read_csv(input_dir / "candidate_cell_metadata.csv.gz")
    candidates = pd.read_csv(input_dir / "candidate_genes.csv")
    peaks = pd.read_csv(input_dir / "candidate_consensus_peaks.csv.gz")

    links = test_control_links(
        genes,
        atac,
        metadata,
        candidates,
        peaks,
        settings,
        int(config["project"]["random_seed"]),
    )

    atac_blocks = []
    rna_blocks = []
    for tf in candidates.TF.unique():
        tf_peaks = peaks.loc[peaks.TF.eq(tf)].copy()
        atac_blocks.append(
            fit_feature_effects(
                atac,
                metadata,
                tf_peaks.peak_index.to_numpy(dtype=int),
                tf_peaks,
                tf,
                references,
                timepoints,
                int(settings["minimum_cells_per_pseudobulk"]),
                "nCount_atac",
            )
        )
        tf_genes = candidates.loc[candidates.TF.eq(tf)].copy()
        rna_blocks.append(
            fit_feature_effects(
                genes,
                metadata,
                tf_genes.gene_index.to_numpy(dtype=int),
                tf_genes,
                tf,
                references,
                timepoints,
                int(settings["minimum_cells_per_pseudobulk"]),
                "nCount_RNA",
            )
        )
    atac_effects = add_fdr_and_sensitivity(
        pd.concat(atac_blocks, ignore_index=True), settings, timepoints
    )
    rna_effects = pd.concat(rna_blocks, ignore_index=True)
    rna_effects["interaction_fdr"] = bh_adjust(
        rna_effects.interaction_p_value.to_numpy()
    )
    for timepoint in timepoints:
        rna_effects[f"effect_fdr_{timepoint}"] = bh_adjust(
            rna_effects[f"effect_p_value_{timepoint}"].to_numpy()
        )

    evidence = atac_effects.merge(links, on="peak_index", validate="one_to_one")
    rna_columns = ["TF_tested", "target", "interaction_fdr"] + [
        item
        for timepoint in timepoints
        for item in (f"effect_{timepoint}", f"effect_fdr_{timepoint}")
    ]
    renamed = {
        column: f"rna_{column}"
        for column in rna_columns
        if column not in {"TF_tested", "target"}
    }
    evidence = evidence.merge(
        rna_effects[rna_columns].rename(columns=renamed),
        on=["TF_tested", "target"],
        validate="many_to_one",
    )
    rna_at_strongest = []
    directional = []
    pattern_correlations = []
    for _, row in evidence.iterrows():
        strongest = row["strongest_timepoint"]
        rna_effect = float(row[f"rna_effect_{strongest}"])
        rna_at_strongest.append(rna_effect)
        atac_pattern = np.array(
            [row[f"effect_day {t}"] for t in [7, 9, 11, 14]], dtype=float
        )
        rna_pattern = np.array(
            [row[f"rna_effect_day {t}"] for t in [7, 9, 11, 14]], dtype=float
        )
        pattern_correlations.append(
            float(np.corrcoef(atac_pattern, rna_pattern)[0, 1])
            if np.std(atac_pattern) > 0 and np.std(rna_pattern) > 0
            else np.nan
        )
        predicted = np.sign(row["link_correlation"] * row["strongest_effect"])
        directional.append(
            predicted != 0
            and np.sign(rna_effect) == predicted
            and abs(rna_effect) >= float(config["state_dependence"]["minimum_on_effect"])
        )
    evidence["rna_effect_at_strongest_atac_timepoint"] = rna_at_strongest
    evidence["atac_rna_effect_pattern_correlation"] = pattern_correlations
    evidence["mechanistic_direction_concordant"] = directional
    evidence["provisional_E2_peak"] = (
        evidence.link_pass
        & evidence.perturbation_sensitive
        & evidence.mechanistic_direction_concordant
    )
    evidence["chromatin_gating_peak"] = (
        evidence.provisional_E2_peak & evidence.state_dependent
    )
    evidence = add_baseline_gating_diagnostics(
        evidence,
        genes,
        atac,
        metadata,
        candidates,
        list(settings["link_control_guides"]),
        timepoints,
    )

    ranked = evidence.sort_values(
        [
            "target",
            "chromatin_gating_peak",
            "provisional_E2_peak",
            "link_pass",
            "perturbation_sensitive",
            "link_correlation",
        ],
        ascending=[True, False, False, False, False, False],
    )
    best = ranked.groupby(["TF", "target", "candidate_role"], as_index=False).first()
    counts = evidence.groupby(["TF", "target"]).agg(
        local_peaks=("peak_index", "size"),
        linked_peaks=("link_pass", "sum"),
        perturbation_sensitive_peaks=("perturbation_sensitive", "sum"),
        provisional_E2_peaks=("provisional_E2_peak", "sum"),
        chromatin_gating_peaks=("chromatin_gating_peak", "sum"),
    ).reset_index()
    best_columns = [
        "TF",
        "target",
        "peak_id",
        "distance_to_tss",
        "link_correlation",
        "link_fdr",
        "strongest_timepoint",
        "strongest_effect",
        "strongest_effect_fdr",
        "rna_effect_at_strongest_atac_timepoint",
        "mechanistic_direction_concordant",
        "atac_rna_effect_pattern_correlation",
        "provisional_E2_peak",
        "chromatin_gating_peak",
    ]
    edge_summary = counts.merge(best[best_columns], on=["TF", "target"])
    edge_summary = edge_summary.merge(
        candidates[["TF", "target", "role"]], on=["TF", "target"]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(
        output_dir / "candidate_peak_mechanism_evidence.csv.gz",
        index=False,
        compression="gzip",
    )
    edge_summary.to_csv(output_dir / "candidate_edge_mechanism_summary.csv", index=False)
    rna_effects.to_csv(output_dir / "candidate_rna_timepoint_effects.csv", index=False)
    summary = {
        "candidate_edges": int(len(edge_summary)),
        "candidate_peaks_tested": int(len(evidence)),
        "control_linked_peaks": int(evidence.link_pass.sum()),
        "atac_perturbation_sensitive_peaks": int(
            evidence.perturbation_sensitive.sum()
        ),
        "provisional_E2_peaks": int(evidence.provisional_E2_peak.sum()),
        "chromatin_gating_peaks": int(evidence.chromatin_gating_peak.sum()),
        "edges_with_provisional_E2": int(
            edge_summary.provisional_E2_peaks.gt(0).sum()
        ),
        "motif_status": "not_yet_tested; provisional_E2 is not final E2",
        "multiple_testing_scope": "193 frozen candidate-local consensus peaks",
    }
    with (output_dir / "chromatin_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--input", type=Path, default=Path("data/processed/targeted_multiome")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/chromatin_mechanism")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
