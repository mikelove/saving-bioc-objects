# Project

This is a Bioconductor workflow package whose sole output is a vignette:
`vignettes/savingBiocObjects.qmd`. The package contains no functions —
it exists to document and explain strategies for saving and sharing
Bioconductor objects, covering:

- R serialization with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and cross-release
  stability via `updateObject()`
- HDF5-backed storage for large assay data (`HDF5Array`)
- Converting `SingleCellExperiment` to/from Python AnnData (`.h5ad`) via
  `anndataR` or `zellkonverter`
- Reading RDS files in Python via `rds2py`
- The `alabaster` / ArtifactDB ecosystem for cross-language, long-term
  serialization
- BED format for genomic ranges, including sidecar files for metadata
  and Seqinfo
- Preserving object-level metadata as JSON sidecars

The target audience is Bioconductor users who need to share objects with
collaborators, deposit data in repositories, or ensure objects remain
readable across Bioconductor releases.

The vignette is rendered into a GitHub Pages site via pkgdown (see
`.github/workflows/pkgdown.yaml`). To render locally, run:

``` r

pkgdown::build_site()
```

## DESCRIPTION Suggests

`Suggests` in DESCRIPTION must include every package used in an
**evaluated** code chunk in the vignette — these are the packages that
get installed in CI and when building the pkgdown site. Packages
referenced only in `eval: false` chunks do not need to be listed.
