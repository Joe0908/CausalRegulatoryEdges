"""Add sequence-level TF motif support to provisional chromatin evidence."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import time
import urllib.request

import numpy as np
import pandas as pd
import yaml


# JASPAR 2024 CORE vertebrate position-frequency matrices.  Values are kept in
# source form so the score calculation and database version remain auditable.
JASPAR_PFMS = {
    "GATA1": {
        "matrix_id": "MA0035.5",
        "A": [3953, 1314, 49692, 67550, 2206, 2567, 7397],
        "C": [62419, 595, 710, 1292, 1238, 65937, 3025],
        "G": [2712, 652, 856, 988, 618, 1471, 765],
        "T": [2744, 69267, 20570, 1998, 67766, 1853, 60641],
    },
    "NFE2": {
        "matrix_id": "MA0841.2",
        "A": [18831, 2, 0, 18831, 37, 0, 1, 18831, 0, 487],
        "C": [9, 8, 0, 11, 14360, 20, 18831, 13, 2689, 9860],
        "G": [1306, 0, 18831, 0, 4472, 3, 0, 8, 10, 5067],
        "T": [5, 18831, 2, 4, 17, 18831, 3, 5, 18831, 3417],
    },
}

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def pwm_from_pfm(pfm: dict[str, list[int]], pseudocount: float = 0.8) -> np.ndarray:
    matrix = np.array([pfm[base] for base in BASES], dtype=float)
    probability = (matrix + pseudocount) / (
        matrix.sum(axis=0, keepdims=True) + 4 * pseudocount
    )
    return np.log2(probability / 0.25)


def scan_sequence(sequence: str, pwm: np.ndarray, threshold: float) -> dict:
    """Scan both strands and return the best relative PWM score and hit count."""
    sequence = sequence.upper()
    width = pwm.shape[1]
    minimum = float(pwm.min(axis=0).sum())
    maximum = float(pwm.max(axis=0).sum())
    denominator = maximum - minimum
    best = -np.inf
    best_position = -1
    best_strand = "."
    hits = 0
    for strand, scanned in (("+", sequence), ("-", sequence.translate(COMPLEMENT)[::-1])):
        for position in range(len(scanned) - width + 1):
            word = scanned[position : position + width]
            if any(base not in BASES for base in word):
                continue
            score = sum(pwm[BASES.index(base), j] for j, base in enumerate(word))
            relative = (score - minimum) / denominator
            if relative >= threshold:
                hits += 1
            if relative > best:
                best = relative
                best_position = position
                best_strand = strand
    return {
        "motif_best_relative_score": float(best),
        "motif_best_position": int(best_position),
        "motif_best_strand": best_strand,
        "motif_hit_count": int(hits),
        "motif_support": bool(best >= threshold),
    }


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(pieces)
                name = line[1:]
                pieces = []
            else:
                pieces.append(line)
    if name is not None:
        sequences[name] = "".join(pieces)
    return sequences


def fetch_peak_sequences(peaks: pd.DataFrame, fasta_path: Path) -> dict[str, str]:
    """Fetch GRCh38 sequences from Ensembl, converting BED to 1-based coordinates."""
    if fasta_path.exists():
        return read_fasta(fasta_path)
    records = []
    for row in peaks.itertuples(index=False):
        chromosome = str(row.chromosome).removeprefix("chr")
        records.append(
            (row.peak_id, f"{chromosome}:{int(row.start) + 1}..{int(row.end)}:1")
        )
    sequences: dict[str, str] = {}
    url = "https://rest.ensembl.org/sequence/region/homo_sapiens"
    for offset in range(0, len(records), 20):
        batch = records[offset : offset + 20]
        body = json.dumps({"regions": [region for _, region in batch]}).encode()
        request = urllib.request.Request(
            url,
            data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            method="POST",
        )
        last_error = None
        for attempt in range(3):
            try:
                response = json.loads(
                    urllib.request.urlopen(request, timeout=45).read().decode()
                )
                if len(response) != len(batch):
                    raise ValueError("Ensembl returned an unexpected sequence count")
                for (peak_id, _), item in zip(batch, response, strict=True):
                    sequences[peak_id] = str(item["seq"]).upper()
                last_error = None
                break
            except Exception as error:  # network retry is intentionally narrow
                last_error = error
                time.sleep(2 * (attempt + 1))
        if last_error is not None:
            raise RuntimeError(f"Ensembl sequence request failed: {last_error}")
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(fasta_path, "wt", encoding="utf-8") as handle:
        for peak_id, sequence in sequences.items():
            handle.write(f">{peak_id}\n{sequence}\n")
    return sequences


def run(config_path: Path, mechanism_dir: Path, sequence_path: Path) -> dict:
    config = load_config(config_path)
    settings = config["chromatin_mechanism"]
    threshold = float(settings["motif_relative_score_threshold"])
    evidence = pd.read_csv(
        mechanism_dir / "candidate_peak_mechanism_evidence.csv.gz"
    )
    sequences = fetch_peak_sequences(evidence, sequence_path)
    records = []
    for row in evidence.itertuples(index=False):
        model = JASPAR_PFMS[row.TF]
        configured_id = settings["motif_models"][row.TF]
        if model["matrix_id"] != configured_id:
            raise ValueError(f"Configured and embedded motif differ for {row.TF}")
        pwm = pwm_from_pfm(model)
        scan = scan_sequence(sequences[row.peak_id], pwm, threshold)
        scan.update(
            {
                "peak_index": int(row.peak_index),
                "motif_matrix_id": model["matrix_id"],
                "peak_sequence_length": len(sequences[row.peak_id]),
            }
        )
        records.append(scan)
    motif = pd.DataFrame(records)
    final = evidence.merge(motif, on="peak_index", validate="one_to_one")
    final["E2_peak"] = final.provisional_E2_peak & final.motif_support
    final["chromatin_gating_E2_peak"] = final.E2_peak & final.state_dependent
    final.to_csv(
        mechanism_dir / "candidate_peak_final_E2_evidence.csv.gz",
        index=False,
        compression="gzip",
    )

    counts = final.groupby(["TF", "target", "candidate_role"]).agg(
        local_peaks=("peak_index", "size"),
        linked_peaks=("link_pass", "sum"),
        atac_sensitive_peaks=("perturbation_sensitive", "sum"),
        provisional_E2_peaks=("provisional_E2_peak", "sum"),
        motif_supported_provisional_peaks=(
            "E2_peak",
            "sum",
        ),
        chromatin_gating_E2_peaks=("chromatin_gating_E2_peak", "sum"),
    ).reset_index()
    best = final.sort_values(
        ["target", "E2_peak", "provisional_E2_peak", "motif_best_relative_score"],
        ascending=[True, False, False, False],
    ).groupby(["TF", "target"], as_index=False).first()
    edge_summary = counts.merge(
        best[
            [
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
                "atac_rna_effect_pattern_correlation",
                "motif_matrix_id",
                "motif_best_relative_score",
                "motif_hit_count",
                "E2_peak",
                "chromatin_gating_E2_peak",
            ]
        ],
        on=["TF", "target"],
    )
    edge_summary.to_csv(mechanism_dir / "candidate_edge_final_E2_summary.csv", index=False)
    summary = {
        "peaks_scanned": int(len(final)),
        "motif_threshold_relative_score": threshold,
        "motif_supported_peaks": int(final.motif_support.sum()),
        "provisional_E2_peaks": int(final.provisional_E2_peak.sum()),
        "final_E2_peaks": int(final.E2_peak.sum()),
        "edges_with_final_E2": int(
            edge_summary.motif_supported_provisional_peaks.gt(0).sum()
        ),
        "chromatin_gating_E2_peaks": int(final.chromatin_gating_E2_peak.sum()),
        "assembly": "GRCh38",
        "coordinate_conversion": "10x BED start + 1 through BED end",
        "motif_source": "JASPAR 2024 CORE vertebrates",
    }
    with (mechanism_dir / "motif_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/mvp.yaml"))
    parser.add_argument(
        "--mechanism-dir", type=Path, default=Path("reports/chromatin_mechanism")
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        default=Path("data/processed/targeted_multiome/candidate_peak_sequences.fa.gz"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.mechanism_dir, args.sequences), indent=2))


if __name__ == "__main__":
    main()
