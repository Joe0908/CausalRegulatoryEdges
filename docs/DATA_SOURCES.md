# Public data sources

Raw source matrices are intentionally excluded from version control. The code
expects the following public data and preserves the original coordinate systems
until an explicitly documented mapping step.

## Primary perturbation atlas

| Field | Value |
|---|---|
| Study | Martin-Rufino et al., *Science* (2025) |
| Assay | Human haematopoietic Perturb-multiome |
| GEO | [GSE274113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274113) |
| Paper | [doi:10.1126/science.ads7951](https://doi.org/10.1126/science.ads7951) |
| Source repository | [sankaranlab/perturb_multiome](https://github.com/sankaranlab/perturb_multiome) |
| Genome build | GRCh38 |

Expected local paths:

```text
data/metadata/GSE274113_annotated_metadata.csv.gz
data/raw/all_timepoints/GSE274113_rep1_filtered_feature_bc_matrix.h5
...
data/raw/all_timepoints/GSE274113_rep16_filtered_feature_bc_matrix.h5
```

The exact 14-library manifest is encoded in `config/mvp.yaml`. GSE274110 is a
supporting PRO-seq experiment and is not mixed into this single-cell analysis.
Exact byte sizes and SHA-256 checksums for the metadata and source HDF5 files
used here are recorded in `data/metadata/data_manifest.csv`.

## Fetal haematopoiesis validation

| Field | Value |
|---|---|
| Study | Ranzoni et al., *Cell Stem Cell* (2021) |
| Assay | Unpaired single-cell RNA-seq and ATAC-seq |
| ArrayExpress | E-MTAB-9067 and E-MTAB-9068 |
| Paper | [doi:10.1016/j.stem.2020.11.015](https://doi.org/10.1016/j.stem.2020.11.015) |
| Source repository | [cvejic-group/integrative-scrna-scatac-human-foetal](https://gitlab.com/cvejic-group/integrative-scrna-scatac-human-foetal) |
| Genome build | GRCh38 |

Expected local root: `external/ranzoni_fetal/`. The exact files are listed under
`external_validation.fetal` in `config/mvp.yaml`.

The fetal ATAC atlas is HSPC-enriched and lacks an annotated mature erythroid
ATAC state. It is therefore used as a coverage-limited trajectory test, not as a
complete test of terminal chromatin ordering.

## Adult terminal erythropoiesis validation

| Field | Value |
|---|---|
| Study | Ludwig et al., *Cell Reports* (2019) |
| Assay | Stage-sorted bulk RNA-seq and ATAC-seq |
| ATAC | [GSE115672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115672) |
| RNA | [GSE115678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115678) |
| Paper | [doi:10.1016/j.celrep.2019.05.046](https://doi.org/10.1016/j.celrep.2019.05.046) |
| Genome build | hg19 |

Expected local root: `external/ludwig_erythropoiesis/`. Frozen GRCh38 candidate
intervals are lifted to hg19 with the UCSC chain specified in `config/mvp.yaml`;
both endpoints must map compatibly and the mapped interval must overlap an
independently called Ludwig ATAC peak.

## Reference resources

- Motifs: JASPAR 2024 CORE vertebrates; GATA1 MA0035.5 and NFE2 MA0841.2.
- Peak sequences and coordinates: Ensembl REST resources for GRCh38.
- Coordinate conversion: UCSC `hg38ToHg19.over.chain.gz`.

Users are responsible for complying with the terms attached to each source
dataset and reference resource.
