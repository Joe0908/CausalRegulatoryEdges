# Interventional validation of single-cell regulatory networks

[![Tests](https://github.com/Joe0908/CausalRegulatoryEdges/actions/workflows/tests.yml/badge.svg)](https://github.com/Joe0908/CausalRegulatoryEdges/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

This repository contains a leakage-controlled workflow for testing whether
TF-target relationships inferred from unperturbed single-cell data remain
supported after genetic perturbation in human erythropoiesis.

## Question

Which observational regulatory edges retain intervention support, and how does
that support depend on developmental context and the evidence layer being
measured?

## Why it matters

Observational gene-regulatory networks can reflect shared cell state or
differentiation rather than a response to perturbing the nominated regulator.
The workflow therefore keeps observational discovery, RNA response, paired
chromatin evidence, and external temporal validation as distinct analyses.

## Approach

1. Infer candidate TF-target edges using non-targeting cells only.
2. Freeze edge rankings and decision rules before examining perturbed cells.
3. Test guide-consistent target-RNA responses with library-aware pseudobulk
   models and leave-one-guide-out checks.
4. Audit collection-time, cell-state, guide-signature, threshold, and
   inferential-unit sensitivity.
5. Evaluate targeted paired RNA/ATAC evidence and stage-resolved external
   datasets without treating association, perturbation response, and chromatin
   support as interchangeable.

Operational evidence labels used by the code are:

| Label | Operational definition |
|---|---|
| E0 | Stable controls-only TF-target association |
| E1 | E0 edge with a guide-consistent target-RNA response |
| E2 | Targeted local chromatin evidence linked to a compatible RNA response |

These labels are specific to this analysis and are not universal causal
categories.

## Data

The primary dataset is the public human haematopoietic Perturb-multiome atlas
GSE274113. External checks use E-MTAB-9067/E-MTAB-9068 and
GSE115672/GSE115678. Raw matrices are not redistributed.

Accessions, genome builds, source links, expected paths, and file checksums are
documented in [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). Frozen parameters
are stored in `config/`.

## Repository structure

| Path | Contents |
|---|---|
| `src/edge_causality/` | Executable discovery, intervention, state, chromatin, validation, and audit modules |
| `config/` | Frozen panels, thresholds, seeds, and analysis contracts |
| `data/metadata/` | Versioned manifests and public metadata needed to locate inputs |
| `reports/` | Machine-readable reference outputs and compact audit summaries |
| `tests/` | Synthetic unit tests that do not require the raw atlases |
| `docs/DATA_SOURCES.md` | Data provenance and acquisition details |
| `docs/REPRODUCIBILITY.md` | Installation and execution order |

## Reproduction

Python 3.10 or newer is required.

```bash
git clone https://github.com/Joe0908/CausalRegulatoryEdges.git
cd CausalRegulatoryEdges
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test,external]'
pytest -q
```

After downloading the public inputs, follow the ordered commands in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). The tests use compact
synthetic fixtures and can be run without the source matrices.

## Headline findings

- The controls-only screen evaluated 43,740 candidate edges from six erythroid
  TFs and retained 619 stable E0 associations.
- The top 1% of the residualized ranking was enriched 2.44-fold for published
  intervention-responsive targets under the prespecified matched-null audit.
- Strict state-matched day-14 testing retained 13 E1 edges; detailed estimates
  and sensitivity results are kept as machine-readable outputs rather than
  narrated in this README.

## Scope and status

This is a completed computational re-analysis of public data. It does not claim
that every supported edge is direct, that motif occurrence proves binding, or
that library batches are independent biological donors. Post-treatment state
adjustment changes the estimand, and the external datasets provide temporal
support rather than an independent causal intervention.

The public repository intentionally contains methods, provenance,
reproducibility assets, and limited reference results. Manuscript drafts,
extended biological interpretation, unpublished hypotheses, and future project
plans are not maintained here.

## Citation and license

Until an archival DOI is available, use [`CITATION.cff`](CITATION.cff). Code is
released under the [MIT License](LICENSE); source data remain subject to their
original terms.
