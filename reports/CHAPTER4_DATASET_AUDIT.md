# Chapter 4 dataset audit: independent erythroid differentiation

## Validation target

Chapter 4 does **not** ask whether an observational atlas independently proves a
causal edge.  It asks whether the temporal mechanism learned from the perturbation
atlas transports to independent human differentiation systems:

> A Chapter-3 perturbation-sensitive GATA1 peak should become accessible before
> the corresponding target gene undergoes its principal erythroid induction.

The candidate sets and thresholds in `config/mvp.yaml` were frozen before target
profiles were inspected.

## Dataset 1: Ranzoni et al. fetal hematopoiesis

- Primary paper: Ranzoni et al., *Cell Stem Cell* (2021),
  https://doi.org/10.1016/j.stem.2020.11.015
- Official analysis/data repository:
  https://gitlab.com/cvejic-group/integrative-scrna-scatac-human-foetal
- RNA: 4,504 annotated cells in the repository PAGA object; 4,463 remain after
  removal of endothelial cells in the authors' trajectory object.
- ATAC: 3,611 cells across three batches; 2,264 predicted MPP and 133 predicted
  Comm-Prog cells form the non-cycling primary comparison.
- The modalities are unpaired.  The dataset can test trajectory-level temporal
  transport, not cell-level peak-to-gene coupling.
- Limitation discovered by the frozen analysis: the ATAC collection is strongly
  HSPC-enriched and does not contain an annotated mature erythroid ATAC state.

Primary state comparison:

- RNA: HSC-MPPs -> MEMPs -> Erythroid cells.
- ATAC: MPP -> Comm-Prog.
- Cycling states are excluded from the primary test.
- Only direct genomic interval overlap counts as primary peak mapping; a nearest
  peak within 2 kb is retained only as sensitivity information.

## Dataset 2: Ludwig et al. adult terminal erythropoiesis

- Primary paper: Ludwig et al., *Cell Reports* (2019),
  https://doi.org/10.1016/j.celrep.2019.05.046
- ATAC: GEO GSE115672; RNA: GEO GSE115678.
- 28 paired RNA/ATAC libraries spanning eight sorted populations, with three or
  four replicates per population from two or three healthy adult donors.
- Stage order: P1 MyP, P2 CFU-E, P3 ProE1, P4 ProE2, P5 BasoE, P6 PolyE,
  P7 OrthoE, P8 Orth/Ret.
- The authors aligned to hg19.  Frozen Chapter-3 hg38 intervals are lifted with
  the UCSC hg38-to-hg19 chain, after which direct interval overlap is required.

Frozen temporal rule:

1. For each feature, define activation as the first post-P1 population reaching
   50% of its increase from P1 to its maximum.
2. Require at least 1.0 log2-CPM dynamic range for both ATAC and RNA.
3. Require ATAC activation at least one population before RNA activation.
4. Require this lead in at least 80% of 1,000 within-population bootstrap draws.
5. An edge passes if at least one of its frozen Chapter-3 E2 peaks passes.

This second dataset is bulk and stage-sorted rather than single-cell.  Its value is
orthogonal replication across donors and full terminal-stage coverage.

