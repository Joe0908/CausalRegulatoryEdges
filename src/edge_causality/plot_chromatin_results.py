"""Plot the key baseline-accessibility versus perturbation-effect diagnostic."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPRESENTATIVES = [
    ("CPEB4", "chr5:173860317-173861249", "positive control; total-effect E2"),
    ("SLC25A37", "chr8:23508478-23509397", "fate-shift comparison; total-effect E2"),
    ("OSBP2", "chr22:30649129-30650044", "primary candidate; linked but not E2"),
]


def plot(input_path: Path, output_path: Path) -> None:
    data = pd.read_csv(input_path)
    timepoints = ["day 7", "day 9", "day 11", "day 14"]
    x = range(len(timepoints))
    figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), sharex=True)
    legend_handles = legend_labels = None
    for axis, (target, peak_id, subtitle) in zip(axes, REPRESENTATIVES, strict=True):
        row = data.loc[data.target.eq(target) & data.peak_id.eq(peak_id)].iloc[0]
        access = [row[f"control_accessible_fraction_{timepoint}"] for timepoint in timepoints]
        atac = [abs(row[f"effect_{timepoint}"]) for timepoint in timepoints]
        rna = [abs(row[f"rna_effect_{timepoint}"]) for timepoint in timepoints]
        axis.plot(x, access, marker="o", color="#222222", label="control accessibility")
        axis.set_ylim(0, 0.65)
        axis.set_ylabel("Accessible-cell fraction", color="#222222")
        axis.tick_params(axis="y", labelcolor="#222222")
        second = axis.twinx()
        second.plot(x, atac, marker="s", color="#d95f02", label="|ATAC effect|")
        second.plot(x, rna, marker="^", color="#1b9e77", label="|RNA effect|")
        second.set_ylim(0, 1.1)
        second.set_ylabel("Absolute log2 effect", color="#555555")
        axis.set_xticks(list(x), ["D7", "D9", "D11", "D14"])
        axis.set_title(f"{target}\n{subtitle}", fontsize=10)
        axis.grid(axis="x", alpha=0.15)
        if axis is axes[0]:
            handles1, labels1 = axis.get_legend_handles_labels()
            handles2, labels2 = second.get_legend_handles_labels()
            legend_handles = handles1 + handles2
            legend_labels = labels1 + labels2
    figure.suptitle(
        "Baseline accessibility does not behave as a simple permissive gate",
        y=0.97,
        fontsize=13,
    )
    figure.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.88),
        ncol=3,
        frameon=False,
    )
    figure.subplots_adjust(left=0.07, right=0.94, bottom=0.15, top=0.68, wspace=0.65)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("reports/chromatin_mechanism/candidate_peak_E2_state_robustness.csv.gz"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/chromatin_mechanism/chromatin_trajectory_diagnostic.png"),
    )
    args = parser.parse_args()
    plot(args.input, args.output)


if __name__ == "__main__":
    main()
