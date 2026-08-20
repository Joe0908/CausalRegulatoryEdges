"""Evaluate observational edge rankings against published intervention truth."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
import yaml


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def matched_null_rates(
    edges: pd.DataFrame,
    selected: np.ndarray,
    truth_column: str,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    selected_rows = edges.loc[selected]
    selected_keys = set(zip(selected_rows.TF, selected_rows.target))
    pools: dict[tuple[str, int, int], np.ndarray] = {}
    fallback: dict[str, np.ndarray] = {}
    for tf, group in edges.groupby("TF", observed=True):
        allowed = ~group.apply(lambda row: (row.TF, row.target) in selected_keys, axis=1)
        fallback[str(tf)] = group.index[allowed].to_numpy()
        for (expression_bin, detection_bin), subgroup in group.loc[allowed].groupby(
            ["expression_bin", "detection_bin"], observed=True
        ):
            pools[(str(tf), int(expression_bin), int(detection_bin))] = subgroup.index.to_numpy()

    truth = edges[truth_column].to_numpy(dtype=float)
    truth_sums = np.zeros(iterations, dtype=float)
    selected_groups = selected_rows.groupby(
        ["TF", "expression_bin", "detection_bin"], observed=True
    ).size()
    for (tf, expression_bin, detection_bin), group_size in selected_groups.items():
        candidates = pools.get((str(tf), int(expression_bin), int(detection_bin)))
        if candidates is None or not len(candidates):
            candidates = fallback[str(tf)]
        sampled = rng.choice(
            candidates, size=(iterations, int(group_size)), replace=True
        )
        truth_sums += truth[sampled].sum(axis=1)
    return truth_sums / len(selected_rows)


def evaluate(
    config_path: Path,
    residual_path: Path,
    grnboost_path: Path,
    features_path: Path,
    author_table_path: Path,
    output_dir: Path,
) -> dict:
    config = load_config(config_path)
    rng = np.random.default_rng(int(config["project"]["random_seed"]))
    residual = pd.read_csv(residual_path)
    grn = pd.read_csv(grnboost_path)
    features = pd.read_csv(features_path)
    edges = residual.merge(
        grn[["TF", "target", "importance", "percentile_within_TF"]],
        on=["TF", "target"],
        validate="one_to_one",
    ).merge(
        features[["feature_key", "detection_fraction", "mean_cpm"]],
        left_on="target",
        right_on="feature_key",
        validate="many_to_one",
    )
    edges["residual_rank_fraction"] = edges.groupby(
        "TF", observed=True
    ).absolute_association.rank(method="first", ascending=False, pct=True)
    edges["residual_rank_score"] = 1 - edges.residual_rank_fraction
    edges["grnboost_rank_score"] = 1 - edges.percentile_within_TF
    edges["consensus_rank_score"] = (
        edges.residual_rank_score + edges.grnboost_rank_score
    ) / 2
    edges["consensus_rank_fraction"] = edges.groupby(
        "TF", observed=True
    ).consensus_rank_score.rank(method="first", ascending=False, pct=True)

    author = pd.read_excel(author_table_path, sheet_name="TF_sensitive_genes", header=2)
    author = author.rename(
        columns={
            "gene_ID": "target_symbol",
            "perturbation_name": "TF",
            "log2FC": "author_log2fc",
            "beta_weight": "author_beta_weight",
            "p_weight": "author_p_weight",
        }
    )
    author = author.loc[author.TF.isin(config["mvp"]["primary_tf_panel"])]
    # Two spreadsheet-autoconverted gene labels collapse to the same Excel date;
    # keep the strongest published row so joins remain deterministic.
    author = author.sort_values("author_p_weight").drop_duplicates(
        ["TF", "target_symbol"], keep="first"
    )
    edges = edges.merge(
        author,
        on=["TF", "target_symbol"],
        how="left",
        validate="many_to_one",
    )
    edges["author_TF_sensitive"] = edges.author_p_weight.notna()
    edges["author_effect_025"] = (
        edges.author_TF_sensitive
        & (edges.author_log2fc.abs() >= config["causal_validation"]["minimum_absolute_log2_fold_change"])
    )
    edges["author_direction_concordant"] = (
        np.sign(edges.signed_association) == -np.sign(edges.author_beta_weight)
    )
    edges["author_supported_concordant"] = (
        edges.author_effect_025 & edges.author_direction_concordant
    )

    edges["expression_bin"] = pd.qcut(
        edges.mean_cpm.rank(method="first"), 10, labels=False
    )
    edges["detection_bin"] = pd.qcut(
        edges.detection_fraction.rank(method="first"), 10, labels=False
    )
    metrics = []
    labels = ["author_TF_sensitive", "author_effect_025"]
    scores = {
        "residualized_association": "absolute_association",
        "GRNBoost2": "importance",
        "consensus": "consensus_rank_score",
    }
    for label in labels:
        y = edges[label].to_numpy(dtype=int)
        for method, score_column in scores.items():
            score = edges[score_column].to_numpy(dtype=float)
            metrics.append(
                {
                    "truth": label,
                    "method": method,
                    "AUROC": float(roc_auc_score(y, score)),
                    "AUPRC": float(average_precision_score(y, score)),
                    "positive_rate": float(y.mean()),
                    "positives": int(y.sum()),
                }
            )

    enrichment = []
    iterations = int(config["edge_discovery"]["matched_null_iterations"])
    rank_columns = {
        "residualized_association": "residual_rank_fraction",
        "GRNBoost2": "percentile_within_TF",
        "consensus": "consensus_rank_fraction",
    }
    for method, rank_column in rank_columns.items():
        for cutoff in config["edge_discovery"]["evaluation_cutoffs_percent"]:
            selected = edges[rank_column].to_numpy() <= float(cutoff) / 100
            for truth in ["author_effect_025", "author_supported_concordant"]:
                observed = float(edges.loc[selected, truth].mean())
                null = matched_null_rates(
                    edges, selected, truth, iterations, rng
                )
                null_mean = float(null.mean())
                enrichment.append(
                    {
                        "method": method,
                        "top_percent": int(cutoff),
                        "truth": truth,
                        "selected_edges": int(selected.sum()),
                        "observed_support_rate": observed,
                        "matched_null_mean": null_mean,
                        "enrichment_ratio": (
                            observed / null_mean if null_mean > 0 else np.nan
                        ),
                        "empirical_p_value": float(
                            (1 + np.sum(null >= observed)) / (iterations + 1)
                        ),
                        "null_q025": float(np.quantile(null, 0.025)),
                        "null_q975": float(np.quantile(null, 0.975)),
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(output_dir / "observational_edges_with_author_truth.csv.gz", index=False, compression="gzip")
    pd.DataFrame(metrics).to_csv(output_dir / "ranking_metrics.csv", index=False)
    pd.DataFrame(enrichment).to_csv(output_dir / "topk_matched_null_enrichment.csv", index=False)
    summary = {
        "candidate_edges": int(len(edges)),
        "author_TF_sensitive_edges": int(edges.author_TF_sensitive.sum()),
        "author_effect_025_edges": int(edges.author_effect_025.sum()),
        "author_supported_concordant_edges": int(
            edges.author_supported_concordant.sum()
        ),
        "metrics": metrics,
        "topk_enrichment": enrichment,
        "direction_note": (
            "Concordance compares observational sign with the opposite of the "
            "published perturbation-score regression beta; stronger perturbation "
            "represents stronger loss of TF function."
        ),
    }
    with (output_dir / "author_truth_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--residual-edges",
        type=Path,
        default=Path("reports/residualized_grn/residualized_signed_edges.csv.gz"),
    )
    parser.add_argument(
        "--grnboost-edges",
        type=Path,
        default=Path("reports/grnboost2/grnboost2_edges.csv.gz"),
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/processed/controls_day14/gene_features.csv.gz"),
    )
    parser.add_argument(
        "--author-table",
        type=Path,
        default=Path("external/paper/NIHMS2076308-supplement-Table_S3.xlsx"),
    )
    parser.add_argument("--output", type=Path, default=Path("reports/author_truth"))
    args = parser.parse_args()
    print(
        json.dumps(
            evaluate(
                args.config,
                args.residual_edges,
                args.grnboost_edges,
                args.features,
                args.author_table,
                args.output,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
