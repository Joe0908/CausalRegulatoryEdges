# Chapter 2 first results: cell-state dependence

## Result in one sentence

Collection-time interactions are detectable, but much of the apparent rescue of
atlas-supported/day-14-unsupported edges is associated with perturbation-induced
movement along the differentiation axis; after that diagnostic, only two E0
edges show both timepoint dependence and within-erythroid-state heterogeneity.

## Frozen analysis hierarchy

1. **Primary:** perturbation-by-collection-time interaction across day 7, 9, 11,
   and 14. Collection time is exogenous to the perturbation.
2. **Fate-shift diagnostic:** train a temporal RNA axis using `NT_1/2/3` profiles
   only, then project all guide pseudobulks and compare targeting guides with
   `AAVS1_1/2/3` within library.
3. **Secondary:** perturbation-by-author-erythroid-state interaction. These labels
   are observed after treatment, so this is a persistence diagnostic rather than
   a primary causal test.

No E1 threshold was changed. `rep1`–`rep16` remain library/batch strata, not
independent donors.

## Data layer

- 14 libraries and all 137,604 author-QC cells were recovered.
- Primary pseudobulk: 882 guide-by-library profiles, 36,601 RNA features.
- Secondary erythroid-state pseudobulk: 4,016 profiles, 115,517 cells, and the
  419 unique genes represented by the 619 E0 edges.
- BFU-E was retained in the raw state aggregation but excluded from the model
  because no AAVS1 BFU-E pseudobulk met the frozen 15-cell minimum.

## Primary perturbation-by-timepoint test

Among 619 E0 edges, 53 (8.6%) had a joint interaction FDR below 0.05.

| TF | E0 edges | Time-dependent | Fraction |
|---|---:|---:|---:|
| GATA1 | 82 | 20 | 24.4% |
| NFE2 | 171 | 16 | 9.4% |
| GFI1B | 81 | 8 | 9.9% |
| TAL1 | 67 | 4 | 6.0% |
| BCL11A | 135 | 3 | 2.2% |
| KLF1 | 83 | 2 | 2.4% |

The prespecified shape classifier assigned 29 edges as gated, one as reversed,
23 as other state-dependent patterns, 15 as guide-unstable, one as constitutive,
and 550 as having no detected interaction. “No detected interaction” is not proof
of invariance.

The 62 author-atlas-effect edges that failed strict day-14 E1 included 15
time-dependent edges (24.2%), versus 38/557 (6.8%) among other edges: odds ratio
4.36, Fisher p = 6.32e-5.

## The important negative control: differentiation-axis movement

The controls-only RNA temporal axis explained 96.4% of variance in ten PCs and
recovered collection time under leave-one-library-out validation with R2 = 0.987
and mean absolute error = 0.23 days. The 95th-percentile null shift was 0.547 days.

Eight TF-by-timepoint shifts passed FDR, null-effect-size, and guide-consistency
requirements: GATA1 at days 7/9/11; GFI1B at days 7/9/11; NFE2 at day 11; and
TAL1 at day 11.

Of the 53 time-dependent E0 edges, 25 were associated with a supported TF-level
trajectory shift at the edge's strongest timepoint. All 15 of the
atlas-supported/day-14-unsupported time-dependent edges were in this group.
After excluding fate-shift-associated edges, the original enrichment disappeared
(odds ratio 0, Fisher p = 0.100).

Therefore the 4.36-fold result cannot currently be presented as evidence of a
cell-intrinsic chromatin gate. It is evidence that observational edges can fail a
late intervention because perturbations change developmental progression and the
edge's total effect is concentrated in an earlier context.

## Secondary within-erythroid-state test

The model included CFU-E, pro-erythroblast, basophilic, polychromatic, and
orthochromatic erythroblast profiles, with library and state main effects. Nine
of 619 edges had a perturbation-by-state interaction FDR below 0.05. Only two also
belonged to the 53 primary time-dependent edges:

| Edge | Time FDR | State FDR | Total-effect pattern | Within-state pattern | Key caveat |
|---|---:|---:|---|---|---|
| GATA1→OSBP2 | 0.00380 | 0.0203 | -0.867 at day 7 to +0.037 at day 14 | -0.685 in CFU-E to +0.088 in orthochromatic | Also fate-shift-associated; day-14 E1 absent |
| NFE2→DCAF11 | 0.0300 | 0.0340 | +0.028 at day 7 to +0.398 at day 14 | -0.058 in CFU-E to +0.454 in orthochromatic | Perturbation sign conflicts with the observational edge |

GATA1→OSBP2 is the cleaner state-dependent causal candidate because the early
knockout effect opposes the positive observational association and agrees with
the author's direction. It still cannot be called a pure cell-intrinsic edge:
both trajectory movement and within-state heterogeneity are present.

NFE2→DCAF11 is a reproducible state-dependent perturbation response, but not a
surviving observational edge under the current sign criterion.

## Positive controls

The held-out GATA2 and SPI1 perturbations produced 126 and 59 transcriptome-wide
time-interaction genes, respectively. Neither TF passed the stronger temporal-axis
shift criterion, so their context dependence is not explained by a large global
differentiation delay in this diagnostic.

## Decision for the project

The original project should continue; no method-packaging pivot is needed. The
second chapter now has a sharper result than “chromatin accessibility matters”:

> An observational edge can disappear under intervention for at least two
> distinct reasons: the perturbation changes which developmental states are
> occupied, or the molecular response changes within a matched erythroid state.

Chapter 3 should test chromatin gating on the small, explicitly separated set of
within-state candidates, led by GATA1→OSBP2, while using fate-shift-associated
erythroid genes such as SPTA1, ALAS2, SLC25A37, and EPOR as a distinct comparison
class rather than mixing both mechanisms.
