# Reference outputs

This directory contains machine-readable outputs from the frozen workflow.
They are retained as numerical checkpoints and audit artifacts, not as a prose
results or discussion section.

| Directory | Output class |
|---|---|
| `metadata_audit/` | Input coverage and representation checks |
| `residualized_grn/`, `grnboost2/` | Controls-only edge rankings |
| `validation/`, `perturbation_score/` | Guide scoring and RNA intervention support |
| `state_dependence/`, `time_resolved_support/` | Context and response-pattern audits |
| `chromatin_mechanism/` | Targeted paired RNA/ATAC estimates |
| `external_validation/` | Stage-resolved external checks |
| `robustness/` | Threshold and inferential-unit sensitivity |
| `temporal_transportability/`, `discovery_selection_bias/` | Frozen post hoc temporal audits |

Files ending in `summary.json` or compact `.csv` tables are the preferred
checkpoints. Larger tables preserve row-level auditability and can be regenerated
from the public inputs using `docs/REPRODUCIBILITY.md`.
