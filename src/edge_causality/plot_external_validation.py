"""Plot representative temporal transport results for Chapter 4."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPRESENTATIVES = [
    ("ALAS2", "chrX:55027806-55028737"),
    ("SLC25A37", "chr8:23508478-23509397"),
    ("CPEB4", "chr5:173860317-173861249"),
]


def normalized_summary(
    data: pd.DataFrame, populations: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    grouped = data.groupby("population").log2_cpm.agg(["mean", "sem"]).reindex(populations)
    dynamic = grouped["mean"].max() - grouped["mean"].min()
    if dynamic == 0:
        return np.zeros(len(grouped)), np.zeros(len(grouped))
    mean = (grouped["mean"] - grouped["mean"].min()) / dynamic
    sem = grouped["sem"] / dynamic
    return mean.to_numpy(dtype=float), sem.to_numpy(dtype=float)


def plot(
    atac_path: Path,
    rna_path: Path,
    evidence_path: Path,
    summary_path: Path,
    output_path: Path,
) -> None:
    atac = pd.read_csv(atac_path)
    rna = pd.read_csv(rna_path)
    evidence = pd.read_csv(evidence_path)
    with summary_path.open(encoding="utf-8") as handle:
        summary = json.load(handle)
    populations = [f"P{index}" for index in range(1, 9)]
    labels = ["MyP", "CFU-E", "ProE1", "ProE2", "BasoE", "PolyE", "OrthoE", "Orth/Ret"]
    figure, axes = plt.subplots(2, 2, figsize=(12.6, 8.2))
    x = np.arange(len(populations))
    for axis, (target, peak_id) in zip(axes.ravel()[:3], REPRESENTATIVES, strict=True):
        atac_block = atac.loc[atac.peak_id.eq(peak_id)]
        rna_block = rna.loc[rna.target.eq(target)]
        atac_mean, atac_sem = normalized_summary(atac_block, populations)
        rna_mean, rna_sem = normalized_summary(rna_block, populations)
        axis.errorbar(
            x,
            atac_mean,
            yerr=atac_sem,
            marker="o",
            linewidth=2,
            color="#d95f02",
            label="ATAC",
        )
        axis.errorbar(
            x,
            rna_mean,
            yerr=rna_sem,
            marker="s",
            linewidth=2,
            color="#1b9e77",
            label="RNA",
        )
        axis.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
        row = evidence.loc[evidence.peak_id.eq(peak_id)].iloc[0]
        axis.set_title(
            f"{target}: ATAC {row.atac_activation_population} → RNA {row.rna_activation_population}\n"
            f"bootstrap lead support = {row.bootstrap_lead_fraction:.3f}"
        )
        axis.set_ylim(-0.08, 1.12)
        axis.set_xticks(x, labels, rotation=40, ha="right")
        axis.set_ylabel("Within-feature relative signal")
        axis.grid(axis="x", alpha=0.12)
    axes[0, 0].legend(frameon=False, loc="lower right")

    axis = axes[1, 1]
    positive_total = int(summary["frozen_E2_peaks"])
    comparison_total = int(summary["frozen_linked_non_E2_peaks"])
    positive_pass = int(summary["temporal_passing_E2_peaks"])
    comparison_pass = int(summary["temporal_passing_linked_non_E2_peaks"])
    rates = [positive_pass / positive_total, comparison_pass / comparison_total]
    bars = axis.bar(
        [0, 1], rates, color=["#6a3d9a", "#bdbdbd"], width=0.62
    )
    axis.set_xticks([0, 1], ["Chapter-3 E2", "Linked non-E2"])
    axis.set_ylabel("Fraction with ATAC-before-RNA")
    axis.set_ylim(0, 1)
    for bar, passed, total in zip(
        bars, [positive_pass, comparison_pass], [positive_total, comparison_total], strict=True
    ):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.04,
            f"{passed}/{total}",
            ha="center",
            fontsize=11,
        )
    axis.set_title(
        "Frozen peak-set discrimination\n"
        f"peak p = {summary['E2_vs_non_E2_temporal_fisher_p_value']:.4g}; "
        f"edge p = {summary['E2_vs_non_E2_edge_fisher_p_value']:.4g}"
    )
    axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(
        "Perturbation-sensitive GATA1 peaks precede target induction\n"
        "in independent adult human erythropoiesis",
        fontsize=14,
        y=0.995,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.94])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    base = Path("reports/external_validation")
    parser.add_argument("--atac", type=Path, default=base / "ludwig_atac_trajectory.csv")
    parser.add_argument("--rna", type=Path, default=base / "ludwig_rna_trajectory.csv")
    parser.add_argument(
        "--evidence", type=Path, default=base / "ludwig_peak_temporal_evidence.csv"
    )
    parser.add_argument(
        "--summary", type=Path, default=base / "ludwig_validation_summary.json"
    )
    parser.add_argument(
        "--output", type=Path, default=base / "terminal_temporal_validation.png"
    )
    args = parser.parse_args()
    plot(args.atac, args.rna, args.evidence, args.summary, args.output)


if __name__ == "__main__":
    main()
