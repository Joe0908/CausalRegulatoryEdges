# MVP first results: when does a regulatory edge survive intervention?

## Executive result

The first analysis supports a selective, top-of-ranking causal signal rather
than broad causal validity of observational networks.

Residualized signed TF-target associations contain intervention-supported edges
at their extreme top ranks, but the initial GRNBoost2 ranking does not. Across
43,740 possible edges from six erythroid TFs to 7,291 expressed genes, the top 1%
of residualized associations has a 12.0% rate of published intervention support
with |log2FC| >= 0.25, compared with 4.9% in expression- and detection-matched
null edges: 2.44-fold enrichment, empirical p = 0.001. GRNBoost2 top 1% shows
1.02-fold enrichment (p = 0.490).

This distinction is local rather than global. The residualized ranking has only
modest genome-wide discrimination (AUROC 0.526; AUPRC 0.048 for effect-size-
qualified intervention support), yet its strongest edges are reproducibly
enriched. That is already a useful answer to the project question: observational
edge confidence is concentrated in a narrow upper tail and depends strongly on
how observational structure is estimated.

## Data integrity

- GSE274113 only; GSE274110 is not mixed into the single-cell MVP.
- 14/14 multiome libraries downloaded and opened successfully.
- 137,604/137,604 author-QC cells recovered from their corresponding H5 files.
- All libraries have an identical ordered set of 36,601 RNA features.
- Day-14 late erythroid analysis: 22,437 cells.
- Strict observational discovery set: 1,304 `NT_1`, `NT_2`, and `NT_3` cells.
- AAVS1 cutting controls are reserved for intervention comparisons.

## Cross-guide perturbation scoring

PCA-matched AAVS1_1 neighbors were used to subtract local control profiles. For
each held-out guide, its 100-gene perturbation vector was learned only from the
other two guides targeting the same TF, projected in the Mixscale style, and
standardized against control-guide cells.

| TF | Effective guides | Passes >=2-guide rule |
|---|---:|:---:|
| NFE2 | 3/3 | Yes |
| KLF1 | 3/3 | Yes |
| BCL11A | 2/3 | Yes |
| GATA1 | 1/3 | No |
| GFI1B | 1/3 | No |
| TAL1 | 0/3 | No |

GATA1, GFI1B, and TAL1 remain informative heterogeneity cases; they are not
silently removed from descriptive comparisons.

## E0 observational network

The controls-only residualized network adjusts for log library depth, library
batch, and late erythroid cell state. One hundred stratified bootstraps were run.
An edge is stable when it falls in the top 5% for its TF in at least 70 of 100
bootstraps.

- Candidate edges tested: 43,740.
- Stable residualized edges: 619.
- Stable by TF: NFE2 171, BCL11A 135, KLF1 83, GATA1 82, GFI1B 81, TAL1 67.
- Stable residualized edges also in GRNBoost2 top 1%: 41.
- Stable residualized edges also in GRNBoost2 top 5%: 192.
- Stable residualized edges also in GRNBoost2 top 10%: 293.

The initial GRNBoost2 network is complete, but its 100-bootstrap stability run is
not yet complete. This is recorded as a remaining analysis rather than implied to
have been done.

## Strict within-day-14 E1 result

Using guide x library x late-state pseudobulks and the frozen criteria (FDR <
0.05, |log2FC| >= 0.25, at least two effective guides, guide-direction agreement,
and leave-one-guide-out direction agreement):

- E1-supported edges: 13.
- Directionally concordant E1 edges: 8.
- All 13 currently arise from NFE2.

This does not meet the pre-specified target of >=200 E1 edges. The threshold has
not been changed. The low count may reflect genuine causal sparsity, the restriction
to terminal erythroid states, and the absence of independent donor replicates in
the public day-14 pseudobulk structure. It must not be interpreted as a definitive
failure until the full timepoint-weighted Mixscale regression is reproduced.

## Published intervention truth evaluation

The authors' Table S3 contains genes significant after Bonferroni correction in
their perturbation-score-weighted, leave-one-feature-out regression across the
atlas. Within the 43,740-edge observational universe:

- 11,468 edges link to a published TF-sensitive gene.
- 1,833 also have |published log2FC| >= 0.25.
- 1,041 of these are directionally concordant with the observational edge sign.

### Expression/detection-matched top-k enrichment

| Method | Top | Effect-qualified enrichment | Empirical p | Direction-concordant enrichment | Empirical p |
|---|---:|---:|---:|---:|---:|
| Residualized association | 1% | 2.44x | 0.001 | 2.75x | 0.001 |
| Residualized association | 5% | 1.82x | 0.001 | 2.02x | 0.001 |
| Residualized association | 10% | 1.51x | 0.001 | 1.66x | 0.001 |
| GRNBoost2 | 1% | 1.02x | 0.490 | 0.94x | 0.618 |
| GRNBoost2 | 5% | 0.98x | 0.623 | 0.96x | 0.621 |
| GRNBoost2 | 10% | 0.88x | 0.963 | 0.76x | 0.998 |
| Consensus rank | 1% | 1.47x | 0.017 | 1.20x | 0.265 |
| Consensus rank | 5% | 1.27x | 0.008 | 1.12x | 0.182 |
| Consensus rank | 10% | 1.23x | 0.001 | 1.22x | 0.015 |

The simple consensus dilutes the strongest residualized signal because GRNBoost2
is uninformative in this first context. We therefore retain methods separately
rather than choosing a combined score after observing the result.

## Interpretation boundary

The published truth is atlas-level and the strict pseudobulk result is day-14
specific. Their difference is not a nuisance: it is the entry point to chapter 2,
which will ask whether edges supported across the atlas disappear specifically in
terminal erythroid states. Chromatin accessibility will then be tested as an
explanation, not merely added as another predictor.

No donor-level replication is claimed from `rep13`, `rep14`, and `rep16`.
