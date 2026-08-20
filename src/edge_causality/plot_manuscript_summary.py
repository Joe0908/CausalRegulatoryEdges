"""Generate compact summary figures for the initial manuscript."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


COLORS = {
    "observational": "#4C78A8",
    "intervention": "#E45756",
    "state": "#F2A541",
    "chromatin": "#7A5195",
    "external": "#2A9D8F",
    "neutral": "#B8B8B8",
}


def study_design(output: Path) -> None:
    figure, axis = plt.subplots(figsize=(14, 4.8))
    axis.set_xlim(0, 14)
    axis.set_ylim(0, 5)
    axis.axis("off")
    blocks = [
        (0.3, "1. Observational\nedge discovery", "43,740 edges\n619 stable E0", COLORS["observational"]),
        (3.05, "2. Perturbation\nvalidation", "13 strict E1\n8 sign-concordant", COLORS["intervention"]),
        (5.8, "3. Cell-state\ndependence", "53 time-dependent\n25 fate-shift-associated", COLORS["state"]),
        (8.55, "4. Chromatin\nmechanism", "11 E2 peaks\n4 GATA1 edges", COLORS["chromatin"]),
        (11.3, "5. External\ntransport", "8/11 E2 peaks\n0/10 non-E2", COLORS["external"]),
    ]
    for index, (x, title, result, color) in enumerate(blocks):
        box = FancyBboxPatch(
            (x, 1.35),
            2.35,
            2.3,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.4,
            edgecolor=color,
            facecolor="white",
        )
        axis.add_patch(box)
        axis.text(x + 1.175, 3.05, title, ha="center", va="center", fontsize=12, weight="bold", color=color)
        axis.text(x + 1.175, 2.0, result, ha="center", va="center", fontsize=11, color="#333333")
        if index < len(blocks) - 1:
            axis.annotate(
                "",
                xy=(x + 2.72, 2.5),
                xytext=(x + 2.4, 2.5),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=1.6),
            )
    axis.text(
        7,
        4.55,
        "From observational confidence to intervention, context, mechanism, and transport",
        ha="center",
        va="center",
        fontsize=15,
        weight="bold",
    )
    axis.text(
        7,
        0.55,
        "Primary question: why do some regulatory relationships survive intervention while others do not?",
        ha="center",
        va="center",
        fontsize=12,
        color="#444444",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(figure)


def network_state_summary(base: Path, output: Path) -> None:
    enrichment = pd.read_csv(base / "author_truth/topk_matched_null_enrichment.csv")
    interactions = pd.read_csv(base / "state_dependence/timepoint_interaction_edges.csv.gz")
    e1 = pd.read_csv(base / "validation/E1_supported_edges.csv")
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

    block = enrichment.loc[
        enrichment.truth.eq("author_effect_025") & enrichment.top_percent.eq(1)
    ]
    order = ["residualized_association", "GRNBoost2", "consensus"]
    labels = ["Residualized", "GRNBoost2", "Consensus"]
    block = block.set_index("method").loc[order]
    axes[0, 0].bar(labels, block.enrichment_ratio, color=[COLORS["observational"], COLORS["neutral"], "#8FB1CF"])
    axes[0, 0].axhline(1, color="#555555", linestyle="--", linewidth=1)
    axes[0, 0].set_ylabel("Matched-null enrichment")
    axes[0, 0].set_title("A  Intervention support at top 1%")
    for index, value in enumerate(block.enrichment_ratio):
        axes[0, 0].text(index, value + 0.08, f"{value:.2f}×", ha="center")

    tf_summary = interactions.groupby("TF", observed=True).state_dependent.agg(["sum", "count"])
    tf_summary["fraction"] = tf_summary["sum"] / tf_summary["count"]
    tf_summary = tf_summary.sort_values("fraction", ascending=False)
    axes[0, 1].bar(tf_summary.index, tf_summary.fraction * 100, color=COLORS["state"])
    axes[0, 1].set_ylabel("Time-dependent E0 edges (%)")
    axes[0, 1].set_title("B  Perturbation effects vary over time")
    axes[0, 1].tick_params(axis="x", rotation=35)

    axes[1, 0].bar(["Strict E1", "Sign-concordant E1"], [len(e1), int(e1.E1_direction_concordant.sum())], color=[COLORS["intervention"], "#B94140"])
    axes[1, 0].set_ylabel("Edges")
    axes[1, 0].set_title("C  Late-erythroid causal support is sparse")
    for index, value in enumerate([len(e1), int(e1.E1_direction_concordant.sum())]):
        axes[1, 0].text(index, value + 0.35, str(value), ha="center")

    values = [25, 28]
    axes[1, 1].bar(["Fate-shift associated", "Not fate-shift associated"], values, color=["#D4841B", "#E9C58D"])
    axes[1, 1].set_ylabel("Time-dependent E0 edges")
    axes[1, 1].set_title("D  Developmental redistribution explains many interactions")
    axes[1, 1].tick_params(axis="x", rotation=18)
    for index, value in enumerate(values):
        axes[1, 1].text(index, value + 0.6, str(value), ha="center")

    for axis in axes.ravel():
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.15)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reports", type=Path, default=Path("reports"))
    parser.add_argument("--output", type=Path, default=Path("manuscript/figures"))
    args = parser.parse_args()
    study_design(args.output / "figure1_study_design.png")
    network_state_summary(args.reports, args.output / "figure2_network_state_summary.png")


if __name__ == "__main__":
    main()
