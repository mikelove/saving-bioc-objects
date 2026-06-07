# Saving Bioconductor Objects for Sharing

## Introduction

A common question in Bioinformatics workflows is: how should I save a
Bioconductor object so that a collaborator can load and use it, or so
that it persists reliably across time? The answer depends on several
factors:

- Longevity: will the file need to be readable in 5 or 10 years?
- Language interoperability: does the recipient use R, Python, or both?
- Object size: is the object small enough to serialize entirely, or does
  it contain large on-disk arrays?
- Reproducibility: should the saved form capture provenance and metadata
  alongside data?

This vignette walks through the main options, their trade-offs, and our
recommendations for common scenarios.

## R serialization

The simplest approach is R’s built-in binary serialization.
[`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) saves a
single object to an `.rds` file;
[`save()`](https://rdrr.io/r/base/save.html) bundles one or more named
objects into an `.RData` (or `.rda`) file.

``` r

library(SummarizedExperiment)

se <- SummarizedExperiment(
  assays = list(counts = matrix(1:12, nrow = 3)),
  colData = DataFrame(condition = c("A", "A", "B", "B"))
)

# Single object
tmp_rds <- tempfile(fileext = ".rds")
saveRDS(se, file = tmp_rds)
se2 <- readRDS(tmp_rds)

# Multiple objects in one file
tmp_rda <- tempfile(fileext = ".RData")
save(se, file = tmp_rda)
load(tmp_rda)
```

Prefer [`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html)
for most use cases. [`save()`](https://rdrr.io/r/base/save.html) /
[`load()`](https://rdrr.io/r/base/load.html) silently overwrites any
object in the calling environment that shares a name, which is a common
source of confusion. The main remaining use case for
[`save()`](https://rdrr.io/r/base/save.html) is `.rda` data files
shipped inside R packages (under `data/`).

Advantages:

- Zero setup — works with any R object.
- Compression applied automatically (gzip by default).
- Round-trips perfectly: the loaded object is identical to the saved
  one.

Disadvantages:

- R-only: Python or other languages cannot read these files without a
  bridge library.
- Large objects are loaded entirely into RAM; there is no lazy / on-disk
  access.

When to use it: quick sharing between R users on the same project,
saving intermediate objects in a pipeline, anything under ~1 GB.

If your workflow downloads `.rds` files from a remote URL, consider
using [`BiocFileCache`](https://bioconductor.org/packages/BiocFileCache)
to cache them locally so they are only fetched once.

### Cross-release stability

There is no guarantee that an S4 object serialized today will work with
older versions of Bioconductor. The most common reasons are that the S4
class is not defined at all in the earlier version, or that it exists
but its definition has changed — new slots added, slots renamed or
removed, or infrastructure moved between packages.

A concrete example: in a recent Bioconductor release, `Seqinfo` was
moved out of `GenomicRanges` into its own package. A `GRanges`
serialized on a machine with the newer setup could not be loaded on a
machine with the older `GenomicRanges` because the class definition for
the embedded `Seqinfo` slot was not found. The reverse was also true:
older `GRanges` objects loaded into a newer session sometimes required
[`updateObject()`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html)
to migrate the internal representation:

``` r

# not evaluated — requires an object saved under an older Bioconductor release
gr <- readRDS("old_granges.rds")
gr <- updateObject(gr, verbose = TRUE)
```

The general lesson is: if you can save your data in a format that is not
tied to a particular version of Bioconductor — or better, not tied to R
at all — you should. For a `GRanges`, for instance, a simple TSV is
often enough:

``` r

library(GenomicRanges)

gr <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = c(100, 200, 300), width = 50)
)

tmp_tsv <- tempfile(fileext = ".tsv")
write.table(as.data.frame(gr), tmp_tsv, sep = "\t", quote = FALSE)

gr2 <- makeGRangesFromDataFrame(
  read.table(tmp_tsv, header = TRUE, sep = "\t"),
  keep.extra.columns = TRUE
)
```

This round-trip survives any Bioconductor version and is readable from
Python or the command line. The `alabaster` ecosystem (described below)
applies the same principle more systematically and with better support
for complex objects, using HDF5 and JSON as the underlying storage
formats.

## HDF5-backed storage

For large assay matrices (e.g., single-cell count matrices with millions
of cells), it is impractical to hold the entire object in RAM. The
[`HDF5Array`](https://bioconductor.org/packages/HDF5Array) package
provides array classes backed by HDF5 files, enabling lazy loading and
out-of-memory computation.

``` r

library(HDF5Array)

tmp_hdf5 <- tempfile()
saveHDF5SummarizedExperiment(se, dir = tmp_hdf5, replace = TRUE)

# Assay data remains on disk until accessed
se_h5 <- loadHDF5SummarizedExperiment(tmp_hdf5)
```

The saved directory contains an HDF5 file with the assay data and an
`.rds` file for the non-assay metadata.

Advantages:

- Assay data is stored on disk; only the chunks you access are read into
  RAM.
- HDF5 is a widely used binary format with readers in Python (`h5py`,
  `anndata`), Julia, C/C++, and more.
- Good for objects with tens of gigabytes of assay data.

Disadvantages:

- The output is a directory, not a single file, which complicates
  transfer (use `tar` or `zip` before sharing).
- The `.rds` envelope for metadata is still R-specific.
- Write performance can be slower than
  [`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) for
  small objects.

When to use it: large single-cell or spatial datasets where you want
on-disk access; workflows shared between R users who need memory
efficiency.

### SummarizedExperiment and AnnData

For single-cell workflows that move between R and Python, the
[zellkonverter](https://bioconductor.org/packages/zellkonverter) package
converts directly between `SingleCellExperiment` and the Python
`AnnData` format (`.h5ad` files used by scanpy and related tools). Like
`HDF5Array`, `.h5ad` is an HDF5-based format.

``` r

library(zellkonverter)
library(SingleCellExperiment)

sce <- as(se, "SingleCellExperiment")

tmp_h5ad <- tempfile(fileext = ".h5ad")
writeH5AD(sce, file = tmp_h5ad)
```

    Installing pyenv ...
    Done! pyenv has been installed to '/home/runner/.local/share/r-reticulate/pyenv/bin/pyenv'.
    Using Python: /home/runner/.pyenv/versions/3.14.0/bin/python3.14
    Creating virtual environment '/home/runner/.cache/R/basilisk/1.24.0/zellkonverter/1.22.0/zellkonverterAnnDataEnv-0.12.3' ... 

    Done!
    Installing packages: pip, wheel, setuptools

    Installing packages: 'anndata==0.12.3', 'h5py==3.15.1', 'natsort==8.4.0', 'numpy==2.3.4', 'pandas==2.3.3', 'scipy==1.16.2'

    Virtual environment '/home/runner/.cache/R/basilisk/1.24.0/zellkonverter/1.22.0/zellkonverterAnnDataEnv-0.12.3' successfully created.

``` r

sce2 <- readH5AD(tmp_h5ad)
```

This is often the most convenient path when collaborators are working in
scanpy, as `.h5ad` is the native format on that side. The trade-off
relative to alabaster is that `.h5ad` is AnnData-specific rather than a
general Bioconductor serialization format.

## The alabaster ecosystem

The [`alabaster`](https://bioconductor.org/packages/alabaster.base)
family of packages provides a language-agnostic, file-based
serialization format. Objects are saved as a directory of standard files
(HDF5, JSON, CSV) with a schema that can be read from R, Python, or any
other language.

``` r

library(alabaster.base)
library(alabaster.se)

tmp_alabaster <- tempfile()
saveObject(se, path = tmp_alabaster)

se3 <- readObject(tmp_alabaster)
```

The `alabaster` umbrella package pulls in support for the most common
Bioconductor classes. Individual sub-packages cover specific classes:
`alabaster.se` for `SummarizedExperiment`, `alabaster.sce` for
`SingleCellExperiment`, and so on.

Advantages:

- Truly language-agnostic: Python readers exist via the `dolomite`
  family of packages, enabling seamless R ↔︎ Python interoperability.
- Built on open standards (HDF5, JSON); inspectable without R.
- Forward-designed for long-term reproducibility and schema versioning.
- A good choice for archives or data portals.

Disadvantages:

- Newer ecosystem; not all Bioconductor classes have `alabaster` support
  yet.
- Requires installing the relevant `alabaster.*` sub-package for each
  class.
- Like HDF5Array, the output is a directory.

When to use it: archival storage, data portal submissions,
cross-language workflows, or any situation where you want the saved
format to be readable without R.

## Exporting genomic ranges to BED

When the object is a `GRanges` or similar ranges object and the goal is
interoperability with other tools (genome browsers, Python, command-line
utilities), exporting to BED format is often more useful than R-specific
serialization. Two options:

``` r

library(rtracklayer)

tmp_bed <- tempfile(fileext = ".bed")
export(gr, tmp_bed)
```

``` r

library(plyranges)

tmp_bed2 <- tempfile(fileext = ".bed")
write_bed(gr, tmp_bed2)
```

### Saving Seqinfo separately

BED files do not store chromosome lengths or genome build information,
so `Seqinfo` is silently dropped on export. To preserve it, write it out
alongside the BED file and restore it on load:

``` r

tmp_seqinfo <- tempfile(fileext = ".csv")
write.csv(as.data.frame(seqinfo(gr)), tmp_seqinfo)

df <- read.csv(tmp_seqinfo, row.names = 1)
si <- Seqinfo(
  seqnames   = rownames(df),
  seqlengths = as.integer(df$seqlengths),
  isCircular = as.logical(df$isCircular),
  genome     = as.character(df$genome)
)
si
```

    Seqinfo object with 1 sequence from an unspecified genome; no seqlengths:
      seqnames seqlengths isCircular genome
      chr1             NA         NA   <NA>

## Saving metadata to JSON

Object-level metadata stored in `metadata(object)` — things like
processing parameters, provenance notes, or experiment descriptors — is
lost in any format that only encodes the ranges or assay data. This
applies whether you are writing a BED file, an HDF5 matrix, or any other
non-R format. Write it to a JSON sidecar file so it travels with the
data:

``` r

library(jsonlite)

metadata(se) <- list(genome = "hg38", pipeline = "v2.1", n_samples = 4L)

tmp_json <- tempfile(fileext = ".json")
writeLines(toJSON(metadata(se), pretty = TRUE, auto_unbox = TRUE), tmp_json)

metadata(se2) <- fromJSON(tmp_json)
```

`toJSON` handles simple R types (lists, vectors, data frames) well, but
complex objects (S4 instances, environments) need to be simplified or
omitted before serializing.

## Summary and recommendations

| Scenario | Recommended approach |
|----|----|
| Quick sharing between R users, same Bioc release | [`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) |
| Shipping small data in a package | [`save()`](https://rdrr.io/r/base/save.html) → `data/*.rda` |
| Large assay matrices, R-only | [`saveHDF5SummarizedExperiment()`](https://rdrr.io/pkg/HDF5Array/man/saveHDF5SummarizedExperiment.html) |
| Cross-language (R + Python) | `alabaster::saveObject()` |
| Long-term archive / data portal | `alabaster::saveObject()` |
| Sharing genomic ranges with non-R tools | `rtracklayer::export()` or [`plyranges::write_bed()`](https://tidyomics.github.io/plyranges/reference/io-bed-write.html) |
| Loading an object from an older Bioc release | [`updateObject()`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html) after [`readRDS()`](https://rdrr.io/r/base/readRDS.html) |

In most new projects we recommend defaulting to
[`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) for
convenience and upgrading to `alabaster` when cross-language access or
archival stability becomes a priority. Be aware that any R-serialized
Bioconductor object may require
[`updateObject()`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html)
when loaded under a different Bioconductor release.

## Session info

``` r

sessionInfo()
```

    R version 4.6.0 (2026-04-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.4 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods
    [8] base

    other attached packages:
     [1] jsonlite_2.0.0              plyranges_1.32.0
     [3] dplyr_1.2.1                 rtracklayer_1.72.0
     [5] alabaster.se_1.12.0         alabaster.base_1.12.0
     [7] SingleCellExperiment_1.34.0 zellkonverter_1.22.0
     [9] HDF5Array_1.40.0            h5mread_1.4.0
    [11] rhdf5_2.56.0                DelayedArray_0.38.2
    [13] SparseArray_1.12.2          S4Arrays_1.12.0
    [15] abind_1.4-8                 Matrix_1.7-5
    [17] SummarizedExperiment_1.42.0 Biobase_2.72.0
    [19] GenomicRanges_1.64.0        Seqinfo_1.2.0
    [21] IRanges_2.46.0              S4Vectors_0.50.1
    [23] BiocGenerics_0.58.1         generics_0.1.4
    [25] MatrixGenerics_1.24.0       matrixStats_1.5.0

    loaded via a namespace (and not attached):
     [1] dir.expiry_1.20.0        rjson_0.2.23             xfun_0.58
     [4] lattice_0.22-9           vctrs_0.7.3              rhdf5filters_1.24.0
     [7] tools_4.6.0              bitops_1.0-9             curl_7.1.0
    [10] parallel_4.6.0           tibble_3.3.1             pkgconfig_2.0.3
    [13] cigarillo_1.2.0          lifecycle_1.0.5          compiler_4.6.0
    [16] Rsamtools_2.28.0         Biostrings_2.80.1        codetools_0.2-20
    [19] htmltools_0.5.9          RCurl_1.98-1.19          alabaster.matrix_1.12.0
    [22] yaml_2.3.12              pillar_1.11.1            crayon_1.5.3
    [25] BiocParallel_1.46.0      basilisk_1.24.0          tidyselect_1.2.1
    [28] digest_0.6.39            restfulr_0.0.16          fastmap_1.2.0
    [31] grid_4.6.0               cli_3.6.6                magrittr_2.0.5
    [34] XML_3.99-0.23            withr_3.0.2              filelock_1.0.3
    [37] rappdirs_0.3.4           rmarkdown_2.31           XVector_0.52.0
    [40] httr_1.4.8               otel_0.2.0               reticulate_1.46.0
    [43] png_0.1-9                evaluate_1.0.5           knitr_1.51
    [46] BiocIO_1.22.0            rlang_1.2.0              Rcpp_1.1.1-1.1
    [49] glue_1.8.1               alabaster.ranges_1.12.0  alabaster.schemas_1.12.0
    [52] R6_2.6.1                 Rhdf5lib_2.0.0           GenomicAlignments_1.48.0
