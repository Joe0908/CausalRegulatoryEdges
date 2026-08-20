# Week 1 data audit and MVP lock

## Decision

Proceed with GSE274113 as the primary Perturb-multiome atlas. Do not mix
GSE274110 into the single-cell MVP: it is a supporting PRO-seq comparison of
GATA1- and AAVS1-edited cells, not the main Perturb-multiome experiment.

The day-14 primary TF panel is:

1. GATA1
2. NFE2
3. GFI1B
4. BCL11A
5. TAL1
6. KLF1

This is a purposefully chosen mechanistic panel, not the six targets with the
largest cell counts. All six are central to late erythroid regulation and have at
least two guides with adequate day-14 representation. GATA2 and SPI1 are held out
as pre-specified positive cases for cell-state dependence: the source paper
reports strong knockdown but low perturbation scores in erythroid lineages because
their strongest effects occur in non-erythroid lineages.

## Verified facts

- The author metadata contains 137,604 QC-passing cells. The CSV also contains
  137,604 fully blank interleaved rows; the loader removes these explicitly.
- Day 14 contains 22,974 QC-passing cells: rep13 = 6,611, rep14 = 7,743,
  rep16 = 8,620.
- The locked late-erythroid subset contains 22,437 cells:
  13,317 orthochromatic and 9,120 polychromatic erythroblasts.
- Each H5 contains 36,601 gene features plus 93,574–100,809 peak features.
- Every author-QC barcode is present in the corresponding H5. The H5 files also
  contain 4,067–4,976 cells per library that the author metadata excludes, so the
  metadata—not the H5 alone—is the analysis whitelist.
- The main atlas uses 19 targeted TFs, three sgRNAs per TF, and six controls
  (three NT and three AAVS1 guides).

## Critical interpretation constraints

1. `rep13`, `rep14`, and `rep16` are handled as library/batch strata. Public GEO
   metadata does not identify them as independent donors, so they must not be used
   to claim donor-level replication.
2. Target-gene mRNA is not a valid editing-efficiency proxy. The source paper
   explicitly used single-cell genotyping because coding disruption can leave
   transcript abundance unchanged.
3. Day 14 is a strong test bed for within-late-erythroid causal edge validation,
   but it cannot by itself establish dependence across the full differentiation
   trajectory. Days 7, 9, and 11 enter the second phase.
4. Candidate edge discovery remains strictly control-only. Perturbed cells are
   sequestered until effect estimation and mechanism testing.
5. Conditioning on post-perturbation cell state can introduce post-treatment
   selection bias. The final analysis must report both total effects and
   within-state effects, and interpret their difference rather than substituting
   one for the other.

## Primary-panel coverage at day 14

| TF | Cells | Guides | Smallest guide | Minimum guide count in a library |
|---|---:|---:|---:|---:|
| NFE2 | 1,880 | 3 | 458 | 123 |
| KLF1 | 1,391 | 3 | 134 | 29 |
| GFI1B | 980 | 3 | 223 | 67 |
| TAL1 | 814 | 3 | 216 | 59 |
| BCL11A | 755 | 3 | 220 | 64 |
| GATA1 | 602 | 3 | 61 | 19 |

GATA1's weakest guide is sparse but still exceeds the pre-specified minimum of 15
cells per library. The two stronger GATA1 guides provide the required redundancy.

## Locked analysis boundary

- Discovery: NT/AAVS1 controls only; infer and bootstrap candidate TF-target edges.
- Validation: targeted cells only after edge ranking is frozen.
- E1 call: FDR < 0.05, |log2FC| >= 0.25, at least two effective guides with the
  same direction, and unchanged direction in leave-one-guide-out analyses.
- State dependence: TF × state interaction plus effect heterogeneity across the
  control-defined trajectory.
- Chromatin mechanism: compare supported and unsupported edges using baseline
  accessibility, TF motif evidence, peak-to-gene linkage, and perturbation-induced
  accessibility change.

## Immediate next executable milestone

Download the RNA/multiome matrices for days 7, 9, and 11; reproduce the authors'
nearest-neighbor perturbation scoring logic; freeze effective guides; then infer
the control-only candidate GRN. The pseudobulk builder already produces 378
replicate × guide × late-erythroid-state profiles across 36,601 genes and serves
as a fast QA layer, not as a replacement for single-cell inference.
