# Reproducibility guide

## Analysis contract

- Observational discovery uses `NT_1`, `NT_2`, and `NT_3` cells only.
- Perturbed cells enter only after rankings, thresholds, and candidate roles
  are frozen.
- Collection time is the primary context; author cell labels are secondary
  post-treatment annotations.
- Libraries are batch strata, not independent donors.
- The global random seed is `20260820`.
- Primary parameters are versioned in `config/mvp.yaml`; temporal and
  discovery-stage sensitivity analyses have separate configuration files.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test,external]'
pytest -q
```

The test suite uses synthetic fixtures and does not require public source
matrices.

## Input placement

Download the accessions and reference resources listed in
[`DATA_SOURCES.md`](DATA_SOURCES.md). The primary paths are encoded in
`config/mvp.yaml`, and `data/metadata/data_manifest.csv` records the byte sizes
and SHA-256 checksums used for the reference run.

## Execution order

Run commands from the repository root.

```bash
# Input audit and controls-only discovery
edge-audit --config config/mvp.yaml
edge-audit --config config/mvp.yaml --all-timepoints --output reports/all_timepoint_audit
edge-build-controls --config config/mvp.yaml --output data/processed/controls_day14
edge-build-controls --config config/mvp.yaml --output data/processed/day14_all_guides --all-guides
edge-score-guides --config config/mvp.yaml --input data/processed/day14_all_guides
edge-residual-grn --config config/mvp.yaml
edge-grnboost2 --config config/mvp.yaml --workers 8
edge-validate --config config/mvp.yaml
edge-evaluate-author --config config/mvp.yaml --output reports/author_truth

# Time, state, and differentiation-axis audits
edge-pseudobulk --config config/mvp.yaml --output data/processed/pseudobulk --all-timepoints-total
edge-state-dependence --config config/mvp.yaml --output reports/state_dependence
edge-response-persistence --config config/mvp.yaml --output reports/state_dependence
edge-time-resolved-support --mode strict --config config/mvp.yaml
edge-trajectory-shift --config config/mvp.yaml --output reports/state_dependence
edge-state-pseudobulk --config config/mvp.yaml --output data/processed/state_pseudobulk
edge-within-state --config config/mvp.yaml --output reports/state_dependence

# Targeted paired RNA/ATAC analysis
edge-build-targeted-multiome --config config/mvp.yaml
edge-chromatin-mechanism --config config/mvp.yaml
edge-motif-support --config config/mvp.yaml
edge-within-state-chromatin --config config/mvp.yaml
edge-e2-library-robustness --config config/mvp.yaml

# Guide-signature and threshold sensitivity
edge-guide-signature-sensitivity --config config/mvp.yaml
edge-robustness-audit --grid config/robustness.yaml --output reports/robustness

# External temporal validation
edge-external-validation --config config/mvp.yaml
edge-terminal-validation --config config/mvp.yaml

# Frozen cross-time and discovery-stage audits
edge-build-controls --config config/mvp.yaml --output data/processed/controls_all_timepoints --all-timepoints
edge-temporal-transportability --config config/temporal_transportability.yaml --output reports/temporal_transportability
edge-discovery-selection-bias --config config/discovery_selection_bias.yaml --mvp-config config/mvp.yaml --output reports/discovery_selection_bias

```

## Output verification

Reference outputs are organized by analysis stage under `reports/`; see
[`reports/README.md`](../reports/README.md). Large raw HDF5 files, processed
cell-level matrices, and external source files are ignored by Git.

When reproducing results, verify:

1. source checksums and genome builds;
2. the installed package version and Python version;
3. unchanged YAML configuration files;
4. test-suite completion;
5. expected row counts and decision fields in the compact JSON/CSV summaries.

Non-significant interaction tests are not equivalence tests, post-treatment
state adjustment is not a substitute for the collection-time estimand, and
motif/peak-gene evidence does not by itself establish direct physical binding.
