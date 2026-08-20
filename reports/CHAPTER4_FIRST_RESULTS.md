# Chapter 4 first results: external stem-cell differentiation validation

## Main result

The Chapter-3 chromatin result transports to an independent adult human
erythroid differentiation system: **8 of 11 frozen E2 peaks show reproducible
ATAC activation before target-RNA activation, compared with 0 of 10 linked
non-E2 peaks.**

At the edge level, 3 of 4 E2 edges pass versus 0 of 6 linked non-E2 edge
comparisons (two-sided Fisher exact p = 0.0333).  The peak-level Fisher result is
p = 0.00103, but the edge-level result is the primary enrichment statistic because
peaks assigned to the same gene are not independent.

## The fetal single-cell atlas was informative but temporally incomplete

In the Ranzoni fetal hematopoiesis atlas:

- 4/11 E2 peaks and 5/10 linked non-E2 peaks had a direct overlap in the atlas
  peak set.
- 0/4 mapped E2 peaks met the frozen MPP-to-Comm-Prog accessibility threshold.
- Nevertheless, mean GATA1 chromVAR activity rose from -0.300 in MPP to 1.418 in
  Comm-Prog (standardized effect 1.005; batch-stratified bootstrap 95% CI
  0.794-1.219).
- ALAS2, SLC25A37, and SPTA1 showed strong delayed RNA induction in annotated
  erythroid cells; the atlas lacks a matching mature erythroid ATAC state.

This pattern is internally coherent: the data detect activation of the GATA
regulatory program, but their ATAC sampling stops before the terminal interval in
which most frozen loci establish accessibility.  We therefore treat the 0/4 result
as a coverage-limited non-replication, not as positive evidence for the mechanism.

## Terminal erythroid validation

All 21 frozen peaks lifted from hg38 to hg19 and overlapped a Ludwig ATAC peak.

| Chapter-3 class | Peaks passing | Edges passing | Interpretation |
|---|---:|---:|---|
| E2 | 8/11 | 3/4 | ATAC-before-RNA is common |
| Linked non-E2 | 0/10 | 0/6 | Association alone does not predict the temporal pattern |

Passing E2 edges:

- **GATA1 -> ALAS2:** the E2 peak activates at P2 (CFU-E), whereas ALAS2 RNA
  reaches its activation midpoint at P4 (ProE2); bootstrap lead support = 0.998.
- **GATA1 -> SLC25A37:** four of five E2 peaks pass.  Their ATAC midpoints occur
  at P2-P3 and the RNA midpoint at P4; support = 0.998-1.000.  The fifth peak is a
  genuine exception, opening at P5, after the RNA midpoint.
- **GATA1 -> CPEB4:** all three E2 peaks activate at P2 and CPEB4 RNA at P5
  (BasoE); support = 1.000 for all three.

Near miss:

- **GATA1 -> SPTA1:** both E2 peaks have the expected point estimate (P2 ATAC,
  P3 RNA), but bootstrap support is 0.717 and 0.742, below the pre-registered 0.80
  threshold.  The edge is therefore not counted as externally validated.

The linked non-E2 controls fail for mechanistically distinct reasons: OSBP2 and a
CPEB4 comparison peak do not have sufficient accessibility dynamic range; the
ALAS2 comparison peak opens at the same stage as RNA; and SPTA1/SLC25A37
comparison peaks open at or after their RNA midpoint.

## Biological interpretation

The combined result supports an **establishment-versus-maintenance model**:

1. GATA1-sensitive regulatory elements establish accessibility early in committed
   erythropoiesis.
2. Their target genes often cross their main expression threshold one to three
   differentiation stages later.
3. Once a locus is broadly accessible, later perturbation can have a smaller
   chromatin effect even while the regulatory program is maintained.

This explains the Chapter-3 inverse relationship between baseline accessibility
and perturbation-effect magnitude more precisely than a simple "open chromatin is
required" gate.  The causal edge is most detectable during establishment, not
necessarily after the mature regulatory state has already formed.

## What this chapter does and does not establish

It establishes temporal transport of a perturbation-defined mechanism across:

- fetal tissue and adult ex vivo differentiation,
- single-cell and bulk stage-resolved assays,
- hg38 and independently called hg19 peak sets,
- and positive E2 versus linked non-E2 frozen controls.

It does not turn either external observational dataset into a second causal
intervention.  The causal component remains the Perturb-seq/multiome result; the
external datasets test whether its temporal chromatin signature recurs in an
independent biological context.

## Reproducible outputs

- `reports/external_validation/external_validation_summary.json`
- `reports/external_validation/external_peak_evidence.csv`
- `reports/external_validation/ludwig_validation_summary.json`
- `reports/external_validation/ludwig_peak_temporal_evidence.csv`
- `reports/external_validation/ludwig_edge_summary.csv`
- `reports/external_validation/terminal_temporal_validation.png`

Commands:

```bash
PYTHONPATH=src python -m edge_causality.external_validation
PYTHONPATH=src python -m edge_causality.terminal_validation
PYTHONPATH=src python -m edge_causality.plot_external_validation
python -m pytest -q
```

