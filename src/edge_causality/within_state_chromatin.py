"""Secondary within-state sensitivity analysis for ATAC perturbation effects.

Author cell states are measured after perturbation, so this is not the primary
causal estimand: conditioning can remove a genuine fate-mediated pathway or
introduce collider bias.  The analysis is nevertheless useful for asking
whether a total ATAC effect remains visible among similarly labelled erythroid
cells, a higher bar for a direct chromatin-mechanism claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as student_t
import yaml

from edge_causality.score_perturbations import bh_adjust


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def fit_wls_condition(
    response: np.ndarray,
    groups: pd.DataFrame,
    target_guides: list[str],
    adjust_state: bool,
) -> tuple[np.ndarray, np.ndarray]:
    condition = groups.guide.isin(target_guides).astype(float).to_numpy()
    categorical = ["replicate"] + (["cell_state"] if adjust_state else [])
    dummy = pd.get_dummies(groups[categorical].astype(str), drop_first=True).to_numpy(
        dtype=float
    )
    design = np.column_stack([np.ones(len(groups)), condition, dummy])
    weights = groups.n_cells.to_numpy(dtype=float)
    weights /= np.median(weights)
    root = np.sqrt(weights)[:, None]
    weighted = design * root
    inverse = np.linalg.pinv(weighted.T @ weighted)
    beta = inverse @ weighted.T @ (response * root)
    residual = response - design @ beta
    df = max(1, len(groups) - np.linalg.matrix_rank(design))
    sigma2 = (weights[:, None] * residual**2).sum(axis=0) / df
    standard_error = np.sqrt(np.maximum(sigma2 * inverse[1, 1], 0))
    statistic = np.divide(
        beta[1],
        standard_error,
        out=np.where(beta[1] == 0, 0.0, np.sign(beta[1]) * np.inf),
        where=standard_error > 0,
    )
    p_value = 2 * student_t.sf(np.abs(statistic), df)
    return beta[1], p_value


def build_state_pseudobulk(
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    tf: str,
    reference_guides: list[str],
    states: list[str],
    minimum_cells: int,
) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    guides = sorted(
        metadata.loc[metadata.target.eq(tf), "perturbation_name"].unique().tolist()
    )
    selected = np.flatnonzero(
        metadata.new_CellType.isin(states).to_numpy()
        & (
            metadata.target.eq(tf)
            | metadata.perturbation_name.isin(reference_guides)
        ).to_numpy()
    )
    block = metadata.iloc[selected].copy()
    keys = pd.MultiIndex.from_frame(
        block[
            ["replicate", "Timepoint", "perturbation_name", "target", "new_CellType"]
        ]
        .astype(str)
        .rename(
            columns={
                "Timepoint": "timepoint",
                "perturbation_name": "guide",
                "new_CellType": "cell_state",
            }
        )
    )
    codes, unique = pd.factorize(keys, sort=True)
    assignment = sparse.csr_matrix(
        (np.ones(len(selected)), (codes, np.arange(len(selected)))),
        shape=(len(unique), len(selected)),
    )
    aggregate = (assignment @ counts[selected]).toarray().astype(float)
    groups = unique.to_frame(index=False)
    groups.columns = ["replicate", "timepoint", "guide", "target", "cell_state"]
    groups["n_cells"] = np.bincount(codes)
    groups["library_total"] = np.bincount(
        codes, weights=block.nCount_atac.to_numpy(dtype=float)
    )
    keep = groups.n_cells.ge(minimum_cells).to_numpy()
    groups = groups.loc[keep].reset_index(drop=True)
    log_cpm = np.log2(
        aggregate[keep]
        / np.maximum(groups.library_total.to_numpy()[:, None], 1)
        * 1_000_000
        + 0.5
    )
    return log_cpm, groups, guides


def run(config_path: Path, input_dir: Path, mechanism_dir: Path) -> dict:
    config = load_config(config_path)
    settings = config["chromatin_mechanism"]
    states = list(config["state_dependence"]["secondary_erythroid_states"])
    references = list(config["data"]["intervention_reference_guides"])
    metadata = pd.read_csv(input_dir / "candidate_cell_metadata.csv.gz")
    atac = sparse.load_npz(input_dir / "candidate_atac_counts.npz").tocsr()
    peaks = pd.read_csv(input_dir / "candidate_consensus_peaks.csv.gz")
    final = pd.read_csv(mechanism_dir / "candidate_peak_final_E2_evidence.csv.gz")

    blocks = []
    for tf in peaks.TF.unique():
        log_cpm, groups, guides = build_state_pseudobulk(
            atac,
            metadata,
            tf,
            references,
            states,
            int(settings["minimum_cells_per_pseudobulk"]),
        )
        feature = peaks.loc[peaks.TF.eq(tf), ["peak_index", "TF", "target"]].copy()
        indices = feature.peak_index.to_numpy(dtype=int)
        unadjusted_effect, unadjusted_p = fit_wls_condition(
            log_cpm[:, indices], groups, guides, adjust_state=False
        )
        adjusted_effect, adjusted_p = fit_wls_condition(
            log_cpm[:, indices], groups, guides, adjust_state=True
        )
        guide_effects = []
        for guide in guides:
            selected_groups = groups.loc[
                groups.guide.isin([guide] + references)
            ].reset_index(drop=True)
            selected_response = log_cpm[
                groups.guide.isin([guide] + references).to_numpy()
            ][:, indices]
            effect, _ = fit_wls_condition(
                selected_response, selected_groups, [guide], adjust_state=True
            )
            guide_effects.append(effect)
        guide_matrix = np.vstack(guide_effects)
        feature["erythroid_unadjusted_atac_effect"] = unadjusted_effect
        feature["erythroid_unadjusted_atac_p_value"] = unadjusted_p
        feature["within_state_atac_effect"] = adjusted_effect
        feature["within_state_atac_p_value"] = adjusted_p
        feature["within_state_consistent_guides"] = np.sum(
            np.sign(guide_matrix) == np.sign(adjusted_effect)[None, :], axis=0
        )
        feature["within_state_groups"] = len(groups)
        blocks.append(feature)
    robustness = pd.concat(blocks, ignore_index=True)
    robustness["within_state_atac_fdr"] = bh_adjust(
        robustness.within_state_atac_p_value.to_numpy()
    )
    robustness["within_state_atac_supported"] = (
        robustness.within_state_atac_fdr.lt(0.05)
        & robustness.within_state_atac_effect.abs().ge(
            float(settings["minimum_absolute_atac_effect"])
        )
        & robustness.within_state_consistent_guides.ge(
            int(settings["minimum_consistent_guides"])
        )
    )
    robustness["absolute_effect_attenuation_after_state_adjustment"] = 1 - np.divide(
        robustness.within_state_atac_effect.abs(),
        robustness.erythroid_unadjusted_atac_effect.abs(),
        out=np.full(len(robustness), np.nan),
        where=robustness.erythroid_unadjusted_atac_effect.abs().to_numpy() > 0,
    )
    combined = final.merge(
        robustness.drop(columns=["TF", "target"]),
        on="peak_index",
        validate="one_to_one",
    )
    combined["E2_total_effect"] = combined.E2_peak
    combined["E2_state_robust"] = (
        combined.E2_total_effect & combined.within_state_atac_supported
    )
    combined.to_csv(
        mechanism_dir / "candidate_peak_E2_state_robustness.csv.gz",
        index=False,
        compression="gzip",
    )
    edge_summary = combined.groupby(["TF", "target", "candidate_role"]).agg(
        total_effect_E2_peaks=("E2_total_effect", "sum"),
        state_robust_E2_peaks=("E2_state_robust", "sum"),
        median_state_attenuation_among_E2=(
            "absolute_effect_attenuation_after_state_adjustment",
            lambda values: float(
                np.nanmedian(
                    values[
                        combined.loc[values.index, "E2_total_effect"].to_numpy()
                    ]
                )
            )
            if combined.loc[values.index, "E2_total_effect"].any()
            else np.nan,
        ),
    ).reset_index()
    edge_summary.to_csv(
        mechanism_dir / "candidate_edge_E2_state_robustness.csv", index=False
    )
    e2 = combined.loc[combined.E2_total_effect]
    summary = {
        "total_effect_E2_peaks": int(len(e2)),
        "total_effect_E2_edges": int(e2.target.nunique()),
        "state_robust_E2_peaks": int(e2.E2_state_robust.sum()),
        "state_robust_E2_edges": int(
            e2.loc[e2.E2_state_robust, "target"].nunique()
        ),
        "median_absolute_effect_attenuation_among_E2": float(
            e2.absolute_effect_attenuation_after_state_adjustment.median()
        ),
        "interpretation": (
            "secondary sensitivity only: author cell state is post-perturbation; "
            "attenuation is compatible with fate mediation but does not prove it"
        ),
    }
    with (mechanism_dir / "within_state_chromatin_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--input", type=Path, default=Path("data/processed/targeted_multiome")
    )
    parser.add_argument(
        "--mechanism-dir", type=Path, default=Path("reports/chromatin_mechanism")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.input, args.mechanism_dir), indent=2))


if __name__ == "__main__":
    main()
