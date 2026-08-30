"""Within-screen pseudobulk validation of controls-only E0 candidate edges."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
import yaml

from edge_causality.score_perturbations import bh_adjust


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def weighted_ols(
    response: np.ndarray, design: np.ndarray, weights: np.ndarray, coefficient: int
) -> tuple[np.ndarray, np.ndarray]:
    root_weight = np.sqrt(weights).reshape(-1, 1)
    x = design * root_weight
    y = response * root_weight
    inverse = np.linalg.pinv(x.T @ x)
    beta = inverse @ x.T @ y
    residual = response - design @ beta
    degrees_freedom = max(1, len(response) - np.linalg.matrix_rank(design))
    sigma2 = (weights[:, None] * residual**2).sum(axis=0) / degrees_freedom
    standard_error = np.sqrt(np.maximum(sigma2 * inverse[coefficient, coefficient], 0))
    t_stat = np.divide(
        beta[coefficient],
        standard_error,
        out=np.zeros_like(beta[coefficient]),
        where=standard_error > 0,
    )
    p_value = 2 * student_t.sf(np.abs(t_stat), degrees_freedom)
    return beta[coefficient], p_value


def make_design(groups: pd.DataFrame, target_guides: list[str]) -> np.ndarray:
    condition = groups.guide.isin(target_guides).astype(float).to_numpy()[:, None]
    categorical = pd.get_dummies(
        groups[["replicate", "cell_type"]].astype(str), drop_first=True
    ).to_numpy(dtype=float)
    return np.column_stack([np.ones(len(groups)), condition, categorical])


def guide_effects(
    log_cpm: np.ndarray,
    groups: pd.DataFrame,
    target_guides: list[str],
    reference_guides: list[str],
) -> dict[str, np.ndarray]:
    output: dict[str, np.ndarray] = {}
    for guide in target_guides:
        differences = []
        guide_rows = groups.index[groups.guide.eq(guide)]
        for row in guide_rows:
            stratum = groups.loc[row, ["replicate", "cell_type"]]
            control_rows = groups.index[
                groups.guide.isin(reference_guides)
                & groups.replicate.eq(stratum.replicate)
                & groups.cell_type.eq(stratum.cell_type)
            ]
            if len(control_rows):
                differences.append(log_cpm[row] - log_cpm[control_rows].mean(axis=0))
        if differences:
            output[guide] = np.mean(differences, axis=0)
    return output


def validate(
    config_path: Path,
    pseudobulk_dir: Path,
    residual_edges_path: Path,
    grnboost_edges_path: Path,
    guide_score_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["causal_validation"]
    z = np.load(pseudobulk_dir / "day14_rna_pseudobulk_counts.npz")
    counts = z["counts"].astype(np.float64)
    gene_name = z["gene_name"].astype(str)
    gene_id = z["gene_id"].astype(str)
    groups = pd.read_csv(pseudobulk_dir / "day14_rna_pseudobulk_groups.csv")
    duplicate = pd.Series(gene_name).duplicated(keep=False).to_numpy()
    feature_key = np.where(duplicate, gene_name + "|" + gene_id, gene_name)
    gene_lookup = {key: i for i, key in enumerate(feature_key)}
    library_total = counts.sum(axis=1, keepdims=True)
    log_cpm = np.log2(counts / np.maximum(library_total, 1) * 1_000_000 + 0.5)

    residual = pd.read_csv(residual_edges_path)
    grn = pd.read_csv(grnboost_edges_path)
    e0 = residual.merge(
        grn[["TF", "target", "importance", "rank_within_TF", "percentile_within_TF"]],
        on=["TF", "target"],
        how="left",
        validate="one_to_one",
    )
    for cutoff in config["edge_discovery"]["evaluation_cutoffs_percent"]:
        e0[f"grnboost_top{int(cutoff)}"] = (
            e0.percentile_within_TF <= float(cutoff) / 100
        )

    guide_scores = pd.read_csv(guide_score_path)
    effective_by_tf = (
        guide_scores.loc[guide_scores.effective_guide]
        .groupby("target", observed=True)
        .guide.apply(list)
        .to_dict()
    )
    reference_guides = list(config["data"]["intervention_reference_guides"])
    all_results = []
    for tf, candidates in e0.groupby("TF", observed=True):
        target_guides = effective_by_tf.get(tf, [])
        if not target_guides:
            block = candidates.copy()
            block["effective_guides_used"] = 0
            block["perturbation_log2fc"] = np.nan
            block["perturbation_p_value"] = np.nan
            block["guide_direction_consistent"] = False
            block["leave_one_guide_out_direction_consistent"] = False
            all_results.append(block)
            continue
        selected_rows = groups.index[
            groups.guide.isin(target_guides + reference_guides)
        ].to_numpy()
        selected_groups = groups.loc[selected_rows].reset_index(drop=True)
        response = log_cpm[selected_rows]
        weights = selected_groups.n_cells.to_numpy(dtype=float)
        weights /= np.median(weights)
        design = make_design(selected_groups, target_guides)
        coefficient, p_value = weighted_ols(response, design, weights, coefficient=1)
        per_guide = guide_effects(
            log_cpm, groups, target_guides, reference_guides
        )
        loo_coefficients = []
        if len(target_guides) >= 2:
            for omitted in target_guides:
                retained = [guide for guide in target_guides if guide != omitted]
                loo_rows = groups.index[
                    groups.guide.isin(retained + reference_guides)
                ].to_numpy()
                loo_groups = groups.loc[loo_rows].reset_index(drop=True)
                loo_weights = loo_groups.n_cells.to_numpy(dtype=float)
                loo_weights /= np.median(loo_weights)
                loo_design = make_design(loo_groups, retained)
                loo_beta, _ = weighted_ols(
                    log_cpm[loo_rows], loo_design, loo_weights, coefficient=1
                )
                loo_coefficients.append(loo_beta)

        block = candidates.copy()
        indices = np.array([gene_lookup[target] for target in block.target])
        effect = coefficient[indices]
        block["effective_guides_used"] = len(target_guides)
        block["effective_guide_names"] = ";".join(target_guides)
        block["perturbation_log2fc"] = effect
        block["perturbation_p_value"] = p_value[indices]
        if per_guide:
            guide_matrix = np.vstack([per_guide[guide][indices] for guide in target_guides])
            block["guide_direction_consistent"] = np.all(
                np.sign(guide_matrix) == np.sign(effect)[None, :], axis=0
            )
            block["minimum_absolute_guide_log2fc"] = np.min(
                np.abs(guide_matrix), axis=0
            )
        else:
            block["guide_direction_consistent"] = False
            block["minimum_absolute_guide_log2fc"] = np.nan
        if loo_coefficients:
            loo = np.vstack([beta[indices] for beta in loo_coefficients])
            block["leave_one_guide_out_direction_consistent"] = np.all(
                np.sign(loo) == np.sign(effect)[None, :], axis=0
            )
        else:
            block["leave_one_guide_out_direction_consistent"] = False
        all_results.append(block)

    results = pd.concat(all_results, ignore_index=True)
    tested = results.perturbation_p_value.notna()
    results["perturbation_fdr"] = np.nan
    results.loc[tested, "perturbation_fdr"] = bh_adjust(
        results.loc[tested, "perturbation_p_value"].to_numpy()
    )
    results["E1_supported"] = (
        (results.effective_guides_used >= int(settings["minimum_effective_guides"]))
        & (results.perturbation_fdr < float(settings["fdr_max"]))
        & (
            results.perturbation_log2fc.abs()
            >= float(settings["minimum_absolute_log2_fold_change"])
        )
        & results.guide_direction_consistent.fillna(False)
        & results.leave_one_guide_out_direction_consistent.fillna(False)
    )
    observational_sign = np.sign(results.signed_association)
    perturbation_sign = np.sign(results.perturbation_log2fc)
    # Knockout support predicts an opposite perturbation sign for an activating
    # observational edge and the same sign for a repressive edge.
    results["direction_concordant_with_knockout"] = (
        observational_sign == -perturbation_sign
    )
    results["E1_direction_concordant"] = (
        results.E1_supported & results.direction_concordant_with_knockout
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "E0_to_E1_edge_matrix.csv.gz", index=False, compression="gzip")
    results.loc[results.E1_supported].to_csv(
        output_dir / "E1_supported_edges.csv", index=False
    )
    e0.to_csv(output_dir / "E0_candidate_edges.csv", index=False)
    summary = {
        "E0_stable_residual_edges": int(len(e0)),
        "E0_GRNBoost2_consensus": {
            f"top_{int(cutoff)}_percent": int(e0[f"grnboost_top{int(cutoff)}"].sum())
            for cutoff in config["edge_discovery"]["evaluation_cutoffs_percent"]
        },
        "effective_guides_by_TF": {
            tf: guides for tf, guides in effective_by_tf.items()
        },
        "E1_supported_edges": int(results.E1_supported.sum()),
        "E1_direction_concordant_edges": int(results.E1_direction_concordant.sum()),
        "E1_by_TF": {
            str(k): int(v)
            for k, v in results.loc[results.E1_supported].TF.value_counts().items()
        },
        "interpretation": (
            "Within-screen library/state pseudobulk evidence; library replicates "
            "are not independent donors."
        ),
    }
    with (output_dir / "validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--pseudobulk", type=Path, default=Path("data/processed/pseudobulk")
    )
    parser.add_argument(
        "--residual-edges",
        type=Path,
        default=Path("reports/residualized_grn/stable_residualized_edges.csv"),
    )
    parser.add_argument(
        "--grnboost-edges",
        type=Path,
        default=Path("reports/grnboost2/grnboost2_edges.csv.gz"),
    )
    parser.add_argument(
        "--guide-scores",
        type=Path,
        default=Path("reports/perturbation_score/guide_crossfit_scores.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("reports/validation"))
    args = parser.parse_args()
    print(
        json.dumps(
            validate(
                args.config,
                args.pseudobulk,
                args.residual_edges,
                args.grnboost_edges,
                args.guide_scores,
                args.output,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
