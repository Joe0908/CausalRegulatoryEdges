"""Library-aware robustness audit for frozen candidate peak-gene links.

The original E2 definition is not changed.  This post-freeze audit asks whether
its controls-only peak-gene links survive random-effects pooling, omission of
each source library, and a descriptive within-library sign check.  Libraries
are batch/replicate strata, not independent biological donors.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, pearsonr

from edge_causality.chromatin_mechanism import covariate_design, residualize
from edge_causality.score_perturbations import bh_adjust
from edge_causality.state_dependence import load_config


def _values(
    genes: sparse.csr_matrix,
    atac: sparse.csr_matrix,
    metadata: pd.DataFrame,
    rows: np.ndarray,
    gene_index: int,
    peak_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    block = metadata.iloc[rows]
    gene_count = genes[rows, gene_index].toarray().ravel()
    peak_count = atac[rows, peak_index].toarray().ravel()
    gene = np.log2(
        gene_count / np.maximum(block.nCount_RNA.to_numpy(float), 1) * 1_000_000
        + 0.5
    )
    peak = np.log1p(
        peak_count / np.maximum(block.nCount_atac.to_numpy(float), 1) * 10_000
    )
    design = covariate_design(block)
    return residualize(peak, design), residualize(gene, design)


def _correlation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), float("nan")
    result = pearsonr(x, y)
    return float(result.statistic), float(result.pvalue)


def _random_effects(r: np.ndarray, n: np.ndarray) -> dict[str, float]:
    keep = np.isfinite(r) & (n > 3)
    r = np.clip(r[keep], -0.999999, 0.999999)
    n = n[keep]
    if len(r) == 0:
        return {"meta_link_correlation": np.nan, "meta_link_p_value": np.nan,
                "meta_tau2": np.nan, "meta_Q": np.nan}
    z = np.arctanh(r)
    variance = 1.0 / (n - 3.0)
    fixed_weight = 1.0 / variance
    fixed = float(np.sum(fixed_weight * z) / np.sum(fixed_weight))
    q = float(np.sum(fixed_weight * (z - fixed) ** 2))
    c = float(np.sum(fixed_weight) - np.sum(fixed_weight**2) / np.sum(fixed_weight))
    tau2 = max(0.0, (q - (len(z) - 1)) / c) if c > 0 else 0.0
    weight = 1.0 / (variance + tau2)
    pooled_z = float(np.sum(weight * z) / np.sum(weight))
    se = float(np.sqrt(1.0 / np.sum(weight)))
    p_value = float(2 * norm.sf(abs(pooled_z / se))) if se > 0 else 1.0
    return {
        "meta_link_correlation": float(np.tanh(pooled_z)),
        "meta_link_p_value": p_value,
        "meta_tau2": tau2,
        "meta_Q": q,
    }


def audit(
    config_path: Path,
    input_dir: Path,
    evidence_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    settings = config["chromatin_mechanism"]
    min_cells = int(settings.get("library_robustness_minimum_control_cells", 30))
    sign_min = float(settings.get("library_robustness_minimum_sign_fraction", 0.80))
    meta_fdr_max = float(settings.get("library_robustness_meta_fdr_max", 0.05))
    min_abs = float(settings["minimum_absolute_link_correlation"])

    genes = sparse.load_npz(input_dir / "candidate_rna_counts.npz").tocsr()
    atac = sparse.load_npz(input_dir / "candidate_atac_counts.npz").tocsr()
    metadata = pd.read_csv(input_dir / "candidate_cell_metadata.csv.gz")
    candidates = pd.read_csv(input_dir / "candidate_genes.csv")
    peaks = pd.read_csv(input_dir / "candidate_consensus_peaks.csv.gz")
    evidence = pd.read_csv(evidence_path)
    gene_lookup = candidates.set_index("target").gene_index.astype(int).to_dict()
    control = metadata.perturbation_name.isin(settings["link_control_guides"]).to_numpy()
    libraries = sorted(metadata.replicate.astype(str).unique(), key=lambda x: int(x[3:]))

    library_records: list[dict] = []
    peak_records: list[dict] = []
    for peak in peaks.itertuples(index=False):
        present = set(str(peak.replicates_present).split(";"))
        effects = []
        for library in libraries:
            rows = np.flatnonzero(control & metadata.replicate.astype(str).eq(library).to_numpy())
            if library not in present or len(rows) < min_cells:
                r = p_value = np.nan
            else:
                x, y = _values(genes, atac, metadata, rows, gene_lookup[peak.target], int(peak.peak_index))
                r, p_value = _correlation(x, y)
            effects.append((library, len(rows), r, p_value))
            library_records.append({"peak_index": int(peak.peak_index), "peak_id": peak.peak_id,
                                    "TF": peak.TF, "target": peak.target, "library": library,
                                    "control_cells": int(len(rows)), "library_link_correlation": r,
                                    "library_link_p_value": p_value})
        r_values = np.array([x[2] for x in effects], float)
        n_values = np.array([x[1] for x in effects], float)
        meta = _random_effects(r_values, n_values)
        pooled = float(evidence.loc[evidence.peak_index.eq(peak.peak_index), "link_correlation"].iloc[0])
        evaluable = np.isfinite(r_values)
        sign_fraction = float(np.mean(np.sign(r_values[evaluable]) == np.sign(pooled))) if evaluable.any() else np.nan
        peak_records.append({"peak_index": int(peak.peak_index), "peak_id": peak.peak_id,
                             "TF": peak.TF, "target": peak.target,
                             "pooled_link_correlation": pooled,
                             "evaluable_libraries": int(evaluable.sum()),
                             "within_library_same_sign_fraction": sign_fraction,
                             **meta})

    peak_table = pd.DataFrame(peak_records)
    finite = peak_table.meta_link_p_value.notna()
    peak_table["meta_link_fdr"] = np.nan
    peak_table.loc[finite, "meta_link_fdr"] = bh_adjust(peak_table.loc[finite, "meta_link_p_value"].to_numpy())
    peak_table["random_effects_link_pass"] = (
        peak_table.meta_link_fdr.lt(meta_fdr_max)
        & peak_table.meta_link_correlation.abs().ge(min_abs)
        & (np.sign(peak_table.meta_link_correlation) == np.sign(peak_table.pooled_link_correlation))
    )
    peak_table["within_library_sign_pass"] = peak_table.within_library_same_sign_fraction.ge(sign_min)

    lolo_records: list[dict] = []
    for omitted in libraries:
        block = []
        for peak in peaks.itertuples(index=False):
            present = set(str(peak.replicates_present).split(";"))
            rows = np.flatnonzero(
                control & metadata.replicate.astype(str).ne(omitted).to_numpy()
                & metadata.replicate.astype(str).isin(present).to_numpy()
            )
            x, y = _values(genes, atac, metadata, rows, gene_lookup[peak.target], int(peak.peak_index))
            r, p_value = _correlation(x, y)
            block.append({"omitted_library": omitted, "peak_index": int(peak.peak_index),
                          "peak_id": peak.peak_id, "TF": peak.TF, "target": peak.target,
                          "control_cells": int(len(rows)), "lolo_link_correlation": r,
                          "lolo_link_p_value": p_value})
        block = pd.DataFrame(block)
        finite = block.lolo_link_p_value.notna()
        block["lolo_link_fdr"] = np.nan
        block.loc[finite, "lolo_link_fdr"] = bh_adjust(block.loc[finite, "lolo_link_p_value"].to_numpy())
        pooled_map = peak_table.set_index("peak_index").pooled_link_correlation
        block["lolo_link_pass"] = (
            block.lolo_link_fdr.lt(float(settings["link_fdr_max"]))
            & block.lolo_link_correlation.abs().ge(min_abs)
            & (np.sign(block.lolo_link_correlation) == np.sign(block.peak_index.map(pooled_map)))
        )
        lolo_records.extend(block.to_dict("records"))
    lolo = pd.DataFrame(lolo_records)
    lolo_summary = lolo.groupby("peak_index", observed=True).agg(
        lolo_refits=("lolo_link_pass", "size"), lolo_passes=("lolo_link_pass", "sum"),
        minimum_lolo_abs_correlation=("lolo_link_correlation", lambda x: float(np.nanmin(np.abs(x)))),
        maximum_lolo_fdr=("lolo_link_fdr", "max"),
    ).reset_index()
    peak_table = peak_table.merge(lolo_summary, on="peak_index", how="left")
    peak_table["all_lolo_link_pass"] = peak_table.lolo_passes.eq(peak_table.lolo_refits)
    e2 = evidence[["peak_index", "E2_peak"]].copy()
    peak_table = peak_table.merge(e2, on="peak_index", how="left")

    edge_table = peak_table.loc[peak_table.E2_peak.fillna(False)].groupby(["TF", "target"], observed=True).agg(
        E2_peaks=("peak_index", "size"),
        random_effects_robust_peaks=("random_effects_link_pass", "sum"),
        all_lolo_robust_peaks=("all_lolo_link_pass", "sum"),
        within_library_sign_robust_peaks=("within_library_sign_pass", "sum"),
    ).reset_index()
    edge_table["retains_random_effects_peak"] = edge_table.random_effects_robust_peaks.gt(0)
    edge_table["retains_all_lolo_peak"] = edge_table.all_lolo_robust_peaks.gt(0)
    edge_table["retains_sign_robust_peak"] = edge_table.within_library_sign_robust_peaks.gt(0)

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(library_records).to_csv(output_dir / "candidate_peak_library_effects.csv.gz", index=False, compression="gzip")
    lolo.to_csv(output_dir / "candidate_peak_leave_one_library_out.csv.gz", index=False, compression="gzip")
    peak_table.to_csv(output_dir / "candidate_peak_E2_library_robustness.csv.gz", index=False, compression="gzip")
    edge_table.to_csv(output_dir / "candidate_edge_E2_library_robustness.csv", index=False)

    e2_peaks = peak_table.loc[peak_table.E2_peak.fillna(False)]
    summary = {
        "frozen_E2_peaks": int(len(e2_peaks)),
        "frozen_E2_edges": int(len(edge_table)),
        "random_effects_robust_E2_peaks": int(e2_peaks.random_effects_link_pass.sum()),
        "all_lolo_robust_E2_peaks": int(e2_peaks.all_lolo_link_pass.sum()),
        "within_library_sign_robust_E2_peaks": int(e2_peaks.within_library_sign_pass.sum()),
        "edges_retaining_random_effects_peak": int(edge_table.retains_random_effects_peak.sum()),
        "edges_retaining_all_lolo_peak": int(edge_table.retains_all_lolo_peak.sum()),
        "edges_retaining_sign_robust_peak": int(edge_table.retains_sign_robust_peak.sum()),
        "libraries": libraries,
        "interpretation": "library-composition robustness, not independent-donor replication",
    }
    with (output_dir / "e2_library_robustness_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument("--input", type=Path, default=Path("data/processed/targeted_multiome"))
    parser.add_argument("--evidence", type=Path, default=Path("reports/chromatin_mechanism/candidate_peak_final_E2_evidence.csv.gz"))
    parser.add_argument("--output", type=Path, default=Path("reports/chromatin_mechanism"))
    args = parser.parse_args()
    print(json.dumps(audit(args.config, args.input, args.evidence, args.output), indent=2))


if __name__ == "__main__":
    main()
