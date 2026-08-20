# Reproducibility guide

## Design safeguards

- Observational discovery uses NT_1/2/3 cells only.
- Perturbed cells enter only after edge rankings and candidate roles are frozen.
- Collection time is the primary context because it is assigned before
  perturbation; author cell labels are post-treatment and are used only for
  persistence sensitivity analyses.
- The primary TF panel, thresholds, candidate roles, genome mapping and external
  temporal rule are recorded in `config/mvp.yaml`.
- AAVS1 guides are intervention controls, not observational discovery controls.
- Libraries are treated as batch strata, not as independent donors.
- Random seed: `20260820`.

## Installation and tests

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test,external]'
pytest -q
```

The tests use compact synthetic inputs and do not require the raw public
matrices.

## Execution order

After placing public data at the paths in `DATA_SOURCES.md`:

```bash
# Audit and controls-only discovery
edge-audit --config config/mvp.yaml
edge-audit --config config/mvp.yaml --all-timepoints --output reports/all_timepoint_audit
edge-build-controls --config config/mvp.yaml --output data/processed/controls_day14
edge-score-guides --config config/mvp.yaml
edge-residual-grn --config config/mvp.yaml
edge-grnboost2 --config config/mvp.yaml --workers 8
edge-validate --config config/mvp.yaml
edge-evaluate-author --config config/mvp.yaml --output reports/author_truth

# Collection time, state and differentiation-axis diagnostics
edge-pseudobulk --config config/mvp.yaml --output data/processed/pseudobulk --all-timepoints-total
edge-state-dependence --config config/mvp.yaml --output reports/state_dependence
edge-trajectory-shift --config config/mvp.yaml --output reports/state_dependence
edge-state-pseudobulk --config config/mvp.yaml --output data/processed/state_pseudobulk
edge-within-state --config config/mvp.yaml --output reports/state_dependence

# Targeted paired RNA/ATAC mechanism
edge-build-targeted-multiome --config config/mvp.yaml
edge-chromatin-mechanism --config config/mvp.yaml
edge-motif-support --config config/mvp.yaml
edge-within-state-chromatin --config config/mvp.yaml
edge-plot-chromatin

# Independent temporal validation
edge-external-validation --config config/mvp.yaml
edge-terminal-validation --config config/mvp.yaml
edge-plot-external-validation

# Manuscript summary figures
edge-plot-manuscript-summary
```

## Interpretation boundaries

The edge labels correspond to different estimands. E1 is a state-matched RNA
response. A collection-time interaction tests context dependence of the total
response. Differentiation-axis movement tests whether a TF perturbation changes
developmental progression. State-adjusted E2 asks whether a chromatin response
persists after conditioning on a post-treatment label and can remove mediated
effects or introduce selection bias. These outputs should not be treated as an
automatic ladder from “false” to “true.”

## Compact checkpoints

The following files summarize every chapter without requiring raw matrices:

- `reports/MVP_FIRST_RESULTS.md`
- `reports/CHAPTER2_FIRST_RESULTS.md`
- `reports/CHAPTER3_FIRST_RESULTS.md`
- `reports/CHAPTER4_DATASET_AUDIT.md`
- `reports/CHAPTER4_FIRST_RESULTS.md`

Large raw HDF5, processed cell-level matrices and external source files remain
ignored by Git.
