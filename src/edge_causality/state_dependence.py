"""Test whether perturbational TF-to-gene effects vary across collection time.

The primary context is the exogenous collection time (day 7/9/11/14), not an
author cell-type label observed after perturbation.  For each TF, a weighted
pseudobulk model compares all three targeting guides with the three AAVS1
reference guides while absorbing library effects.  State dependence is tested
with a joint perturbation-by-timepoint interaction F test.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution
from scipy.stats import fisher_exact
from scipy.stats import t as student_t
import yaml

from edge_causality.score_perturbations import bh_adjust


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def feature_keys(gene_name: np.ndarray, gene_id: np.ndarray) -> np.ndarray:
    duplicate = pd.Series(gene_name).duplicated(keep=False).to_numpy()
    return np.where(duplicate, gene_name + "|" + gene_id, gene_name)


def interaction_design(
    groups: pd.DataFrame,
    target_guides: list[str],
    ordered_timepoints: list[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], list[int]]:
    """Create reduced/full designs and timepoint-specific effect contrasts."""
    condition = groups.guide.isin(target_guides).astype(float).to_numpy()
    replicate = pd.get_dummies(
        groups.replicate.astype(str), drop_first=True
    ).to_numpy(dtype=float)
    reduced = np.column_stack([np.ones(len(groups)), condition, replicate])

    interaction_columns = []
    contrasts: dict[str, np.ndarray] = {}
    for timepoint in ordered_timepoints[1:]:
        indicator = groups.timepoint.eq(timepoint).astype(float).to_numpy()
        interaction_columns.append(condition * indicator)
    full = np.column_stack([reduced, *interaction_columns])
    interaction_indices = list(range(reduced.shape[1], full.shape[1]))

    for timepoint in ordered_timepoints:
        contrast = np.zeros(full.shape[1], dtype=float)
        contrast[1] = 1.0
        if timepoint != ordered_timepoints[0]:
            contrast[interaction_indices[ordered_timepoints[1:].index(timepoint)]] = 1.0
        contrasts[timepoint] = contrast
    return reduced, full, contrasts, interaction_indices


def fit_interaction_model(
    response: np.ndarray,
    reduced: np.ndarray,
    full: np.ndarray,
    weights: np.ndarray,
    contrasts: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Fit vectorized WLS and return context effects plus the joint F test."""
    root_weight = np.sqrt(weights).reshape(-1, 1)
    x_full = full * root_weight
    x_reduced = reduced * root_weight
    y = response * root_weight

    inverse_full = np.linalg.pinv(x_full.T @ x_full)
    beta_full = inverse_full @ x_full.T @ y
    beta_reduced = np.linalg.pinv(x_reduced.T @ x_reduced) @ x_reduced.T @ y
    residual_full = response - full @ beta_full
    residual_reduced = response - reduced @ beta_reduced
    sse_full = (weights[:, None] * residual_full**2).sum(axis=0)
    sse_reduced = (weights[:, None] * residual_reduced**2).sum(axis=0)

    rank_full = np.linalg.matrix_rank(full)
    rank_reduced = np.linalg.matrix_rank(reduced)
    df_denominator = max(1, len(response) - rank_full)
    df_numerator = max(1, rank_full - rank_reduced)
    numerator = np.maximum(sse_reduced - sse_full, 0) / df_numerator
    denominator = sse_full / df_denominator
    f_stat = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )
    interaction_p = f_distribution.sf(f_stat, df_numerator, df_denominator)

    sigma2 = sse_full / df_denominator
    output: dict[str, np.ndarray] = {
        "interaction_f_stat": f_stat,
        "interaction_p_value": interaction_p,
        "model_df_numerator": np.repeat(df_numerator, response.shape[1]),
        "model_df_denominator": np.repeat(df_denominator, response.shape[1]),
    }
    for timepoint, contrast in contrasts.items():
        effect = contrast @ beta_full
        variance_multiplier = float(contrast @ inverse_full @ contrast)
        standard_error = np.sqrt(np.maximum(sigma2 * variance_multiplier, 0))
        t_stat = np.divide(
            effect,
            standard_error,
            out=np.zeros_like(effect),
            where=standard_error > 0,
        )
        output[f"effect_{timepoint}"] = effect
        output[f"effect_se_{timepoint}"] = standard_error
        output[f"effect_p_value_{timepoint}"] = 2 * student_t.sf(
            np.abs(t_stat), df_denominator
        )
    return output


def per_guide_context_effects(
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    gene_indices: np.ndarray,
    target_guides: list[str],
    reference_guides: list[str],
    ordered_timepoints: list[str],
) -> dict[str, np.ndarray]:
    """Estimate each targeting guide against within-library AAVS1 controls."""
    output: dict[str, np.ndarray] = {}
    for timepoint in ordered_timepoints:
        effects = []
        for guide in target_guides:
            differences = []
            guide_rows = groups.index[
                groups.guide.eq(guide) & groups.timepoint.eq(timepoint)
            ]
            for row in guide_rows:
                reference_rows = groups.index[
                    groups.guide.isin(reference_guides)
                    & groups.replicate.eq(groups.loc[row, "replicate"])
                ]
                if len(reference_rows):
                    differences.append(
                        log_cpm[row, gene_indices]
                        - log_cpm[reference_rows][:, gene_indices].mean(axis=0)
                    )
            if differences:
                effects.append(np.mean(differences, axis=0))
            else:
                effects.append(np.full(len(gene_indices), np.nan))
        output[timepoint] = np.vstack(effects)
    return output


def classify_edge(row: pd.Series, settings: dict, timepoints: list[str]) -> str:
    effects = np.array([row[f"effect_{timepoint}"] for timepoint in timepoints])
    consistent = np.array(
        [row[f"consistent_guides_{timepoint}"] for timepoint in timepoints]
    )
    on = float(settings["minimum_on_effect"])
    off = float(settings["maximum_off_effect"])
    minimum_range = float(settings["minimum_effect_range"])
    minimum_guides = int(settings["minimum_consistent_guides"])
    significant = row.interaction_fdr < float(settings["interaction_fdr_max"])
    strongest = int(np.argmax(np.abs(effects)))

    if consistent[strongest] < minimum_guides:
        return "unstable"
    if significant and effects.max() >= on and effects.min() <= -on:
        return "reversed"
    if significant and np.max(np.abs(effects)) >= on and np.min(np.abs(effects)) <= off:
        return "gated"
    if (
        significant
        and np.max(np.abs(effects)) >= on
        and np.ptp(np.abs(effects)) >= minimum_range
    ):
        return "amplified"
    if (
        not significant
        and np.all(np.abs(effects) >= on)
        and (np.all(effects > 0) or np.all(effects < 0))
        and np.all(consistent >= minimum_guides)
    ):
        return "constitutive"
    if significant:
        return "state_dependent_other"
    return "no_detected_interaction"


def fit_tf(
    tf: str,
    genes: pd.DataFrame,
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    reference_guides: list[str],
    ordered_timepoints: list[str],
    minimum_cells: int,
) -> pd.DataFrame:
    target_guides = sorted(
        groups.loc[groups.target.eq(tf), "guide"].dropna().unique().tolist()
    )
    selected = groups.index[
        groups.guide.isin(target_guides + reference_guides)
        & (groups.n_cells >= minimum_cells)
    ].to_numpy()
    selected_groups = groups.loc[selected].reset_index(drop=True)
    if len(target_guides) < 2 or selected_groups.timepoint.nunique() < len(ordered_timepoints):
        raise ValueError(f"Insufficient guide/context coverage for {tf}")

    reduced, full, contrasts, _ = interaction_design(
        selected_groups, target_guides, ordered_timepoints
    )
    weights = selected_groups.n_cells.to_numpy(dtype=float)
    weights /= np.median(weights)
    indices = genes.feature_index.to_numpy(dtype=int)
    model = fit_interaction_model(
        log_cpm[selected][:, indices], reduced, full, weights, contrasts
    )
    output = genes.copy().reset_index(drop=True)
    output.insert(0, "TF", tf)
    output["targeting_guides"] = ";".join(target_guides)
    for name, values in model.items():
        output[name] = values

    per_guide = per_guide_context_effects(
        log_cpm,
        groups,
        indices,
        target_guides,
        reference_guides,
        ordered_timepoints,
    )
    for timepoint in ordered_timepoints:
        pooled = output[f"effect_{timepoint}"].to_numpy()
        guide_matrix = per_guide[timepoint]
        output[f"consistent_guides_{timepoint}"] = np.sum(
            np.sign(guide_matrix) == np.sign(pooled)[None, :], axis=0
        )
        output[f"guide_effect_range_{timepoint}"] = (
            np.nanmax(guide_matrix, axis=0) - np.nanmin(guide_matrix, axis=0)
        )
    return output


def run(
    config_path: Path,
    pseudobulk_dir: Path,
    edge_matrix_path: Path,
    author_truth_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["state_dependence"]
    timepoints = list(settings["ordered_timepoints"])
    reference_guides = list(config["data"]["intervention_reference_guides"])

    z = np.load(pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_counts.npz")
    counts = z["counts"].astype(np.float64)
    gene_name = z["gene_name"].astype(str)
    gene_id = z["gene_id"].astype(str)
    keys = feature_keys(gene_name, gene_id)
    key_to_index = {key: i for i, key in enumerate(keys)}
    groups = pd.read_csv(
        pseudobulk_dir / "all_timepoint_total_rna_pseudobulk_groups.csv"
    )
    library_total = counts.sum(axis=1, keepdims=True)
    cpm = counts / np.maximum(library_total, 1) * 1_000_000
    log_cpm = np.log2(cpm + 0.5)

    e0 = pd.read_csv(edge_matrix_path)
    e0["feature_index"] = e0.target.map(key_to_index)
    if e0.feature_index.isna().any():
        missing = e0.loc[e0.feature_index.isna(), "target"].tolist()
        raise ValueError(f"E0 targets absent from pseudobulk features: {missing[:5]}")

    blocks = []
    for tf, edges in e0.groupby("TF", observed=True):
        gene_table = edges[["target", "feature_index"]].copy()
        block = fit_tf(
            str(tf),
            gene_table,
            log_cpm,
            groups,
            reference_guides,
            timepoints,
            int(settings["minimum_cells_per_pseudobulk"]),
        )
        blocks.append(block)
    results = pd.concat(blocks, ignore_index=True)
    results["interaction_fdr"] = bh_adjust(results.interaction_p_value.to_numpy())
    for timepoint in timepoints:
        results[f"effect_fdr_{timepoint}"] = bh_adjust(
            results[f"effect_p_value_{timepoint}"].to_numpy()
        )
    results = results.merge(
        e0.drop(columns=["feature_index"]), on=["TF", "target"], validate="one_to_one"
    )
    results["edge_class"] = results.apply(
        classify_edge, axis=1, settings=settings, timepoints=timepoints
    )
    for timepoint in timepoints:
        results[f"knockout_direction_concordant_{timepoint}"] = (
            np.sign(results.signed_association)
            == -np.sign(results[f"effect_{timepoint}"])
        )

    author = pd.read_csv(author_truth_path)[
        [
            "TF",
            "target",
            "author_TF_sensitive",
            "author_effect_025",
            "author_direction_concordant",
        ]
    ]
    results = results.merge(author, on=["TF", "target"], how="left")
    results["atlas_supported_day14_unsupported"] = (
        results.author_effect_025.fillna(False) & ~results.E1_supported.fillna(False)
    )
    results["state_dependent"] = results.interaction_fdr < float(
        settings["interaction_fdr_max"]
    )

    # GATA2 and SPI1 were frozen as positive controls. They are evaluated
    # transcriptome-wide but kept separate from the controls-only E0 edge set.
    mean_cpm = cpm.mean(axis=0)
    eligible = (mean_cpm >= 1.0) & ~np.char.startswith(gene_name, "MT-")
    positive_gene_table = pd.DataFrame(
        {
            "target": keys[eligible],
            "feature_index": np.flatnonzero(eligible),
        }
    )
    positive_blocks = []
    for tf in settings["positive_control_tfs"]:
        positive_blocks.append(
            fit_tf(
                tf,
                positive_gene_table,
                log_cpm,
                groups,
                reference_guides,
                timepoints,
                int(settings["minimum_cells_per_pseudobulk"]),
            )
        )
    positive = pd.concat(positive_blocks, ignore_index=True)
    positive["interaction_fdr"] = bh_adjust(positive.interaction_p_value.to_numpy())
    for timepoint in timepoints:
        positive[f"effect_fdr_{timepoint}"] = bh_adjust(
            positive[f"effect_p_value_{timepoint}"].to_numpy()
        )
    positive["edge_class"] = positive.apply(
        classify_edge, axis=1, settings=settings, timepoints=timepoints
    )

    target_group = results.atlas_supported_day14_unsupported
    table = pd.crosstab(target_group, results.state_dependent).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    odds_ratio, enrichment_p = fisher_exact(table.to_numpy())
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir / "timepoint_interaction_edges.csv.gz",
        index=False,
        compression="gzip",
    )
    positive.to_csv(
        output_dir / "positive_control_transcriptome.csv.gz",
        index=False,
        compression="gzip",
    )
    class_counts = {
        str(key): int(value) for key, value in results.edge_class.value_counts().items()
    }
    positive_summary = {}
    for tf, block in positive.groupby("TF", observed=True):
        positive_summary[str(tf)] = {
            "genes_tested": int(len(block)),
            "interaction_fdr_lt_0_05": int(
                (block.interaction_fdr < float(settings["interaction_fdr_max"])).sum()
            ),
            "edge_classes": {
                str(key): int(value)
                for key, value in block.edge_class.value_counts().items()
            },
        }
    summary = {
        "primary_context": "collection_timepoint",
        "timepoints": timepoints,
        "E0_edges_tested": int(len(results)),
        "state_dependent_interaction_fdr_lt_0_05": int(results.state_dependent.sum()),
        "edge_classes": class_counts,
        "atlas_supported_day14_unsupported": int(target_group.sum()),
        "atlas_supported_day14_unsupported_state_dependent": int(
            (target_group & results.state_dependent).sum()
        ),
        "heterogeneity_enrichment_odds_ratio": float(odds_ratio),
        "heterogeneity_enrichment_fisher_p_value": float(enrichment_p),
        "positive_controls": positive_summary,
        "limitations": [
            "replicate labels are library/batch strata, not independent donors",
            "total effects may include perturbation-induced lineage-composition shifts",
            "cell-type-conditioned and control-trajectory analyses remain secondary",
            "lack of a significant interaction is not proof of context invariance",
        ],
    }
    with (output_dir / "state_dependence_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--pseudobulk", type=Path, default=Path("data/processed/pseudobulk")
    )
    parser.add_argument(
        "--edge-matrix",
        type=Path,
        default=Path("reports/validation/E0_to_E1_edge_matrix.csv.gz"),
    )
    parser.add_argument(
        "--author-truth",
        type=Path,
        default=Path("reports/author_truth/observational_edges_with_author_truth.csv.gz"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("reports/state_dependence")
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.config,
                args.pseudobulk,
                args.edge_matrix,
                args.author_truth,
                args.output,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
