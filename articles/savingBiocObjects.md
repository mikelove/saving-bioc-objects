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

### saveRDS and readRDS

The simplest approach is R’s built-in binary serialization.
[`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) saves a
single object to an `.rds` file;
[`save()`](https://rdrr.io/r/base/save.html) bundles one or more named
objects into an `.RData` (or `.rda`) file.

``` r

library(SummarizedExperiment)

se <- SummarizedExperiment(
  assays = list(counts = matrix(1:12, nrow = 3)),
  colData = DataFrame(
    condition = c("A", "A", "B", "B"), 
    row.names=1:4
    ),
  rowData = DataFrame(
    gene = c("gene1","gene2","gene3"), 
    row.names=1:3
    )
)

# Single object
tmp_rds <- tempfile(fileext = ".rds")
saveRDS(se, file = tmp_rds)
se_from_rds <- readRDS(tmp_rds)
se_from_rds
```

    class: SummarizedExperiment
    dim: 3 4
    metadata(0):
    assays(1): counts
    rownames(3): 1 2 3
    rowData names(1): gene
    colnames(4): 1 2 3 4
    colData names(1): condition

``` r

# Multiple objects in one file
tmp_rda <- tempfile(fileext = ".RData")
save(se, file = tmp_rda)
load(tmp_rda)  # restores 'se' by name into the current environment
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
names(gr) <- c("peak1", "peak2", "peak3")
gr$score  <- c(500, 800, 300)   # standard BED score column
gr$log2fc <- c(1.2, -0.5, 2.1) # extra metadata column

tmp_tsv <- tempfile(fileext = ".tsv")
write.table(as.data.frame(gr), tmp_tsv, sep = "\t", quote = FALSE)

gr_from_tsv <- makeGRangesFromDataFrame(
  read.table(tmp_tsv, header = TRUE, sep = "\t"),
  keep.extra.columns = TRUE
)
gr_from_tsv
```

    GRanges object with 3 ranges and 2 metadata columns:
            seqnames    ranges strand |     score    log2fc
               <Rle> <IRanges>  <Rle> | <integer> <numeric>
      peak1     chr1   100-149      * |       500       1.2
      peak2     chr1   200-249      * |       800      -0.5
      peak3     chr1   300-349      * |       300       2.1
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

This round-trip survives any Bioconductor version and is readable from
Python or the command line. The `alabaster` ecosystem (described below)
applies the same principle more systematically and with better support
for complex objects, using HDF5 and JSON as the underlying storage
formats.

## HDF5-backed storage

### Saving with HDF5Array

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
se_from_hdf5 <- loadHDF5SummarizedExperiment(tmp_hdf5)
se_from_hdf5
```

    class: SummarizedExperiment
    dim: 3 4
    metadata(0):
    assays(1): counts
    rownames(3): 1 2 3
    rowData names(1): gene
    colnames(4): 1 2 3 4
    colData names(1): condition

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

# not evaluated — zellkonverter installs a full Python environment via basilisk
# on first use, which takes too long in CI
library(zellkonverter)
library(SingleCellExperiment)

sce <- as(se, "SingleCellExperiment")

tmp_h5ad <- tempfile(fileext = ".h5ad")
writeH5AD(sce, file = tmp_h5ad)

sce_from_h5ad <- readH5AD(tmp_h5ad)
sce_from_h5ad
```

This is often the most convenient path when collaborators are working in
scanpy, as `.h5ad` is the native format on that side. The trade-off
relative to alabaster is that `.h5ad` is AnnData-specific rather than a
general Bioconductor serialization format.

## The alabaster ecosystem

### saveObject and readObject

The [`alabaster`](https://bioconductor.org/packages/alabaster.base)
family of packages is part of the broader
[ArtifactDB](https://github.com/ArtifactDB) project, which provides a
multi-language system for storing and retrieving analysis-ready
Bioconductor objects. The core idea is to save objects as directories of
standard files (HDF5, JSON, CSV) whose format is defined by explicit,
versioned specifications — meaning the saved form is readable without R,
and can evolve over time without breaking previously saved objects.

``` r

library(alabaster.base)
library(alabaster.se)

tmp_alabaster <- tempfile()
saveObject(se, path = tmp_alabaster)

se_from_alabaster <- readObject(tmp_alabaster)
se_from_alabaster
```

    class: SummarizedExperiment
    dim: 3 4
    metadata(0):
    assays(1): counts
    rownames(3): 1 2 3
    rowData names(1): gene
    colnames(4): 1 2 3 4
    colData names(1): condition

The `alabaster` umbrella package pulls in support for the most common
Bioconductor classes. Individual sub-packages cover specific classes:
`alabaster.se` for `SummarizedExperiment`, `alabaster.sce` for
`SingleCellExperiment`, and so on.

### Validation with takane

A key part of the ArtifactDB design is that saved directories can be
independently validated against the format specification. This is
handled by [takane](https://github.com/ArtifactDB/takane), a C++ library
that maintains separate, versioned specifications for 30+ Bioconductor
object types. Calling `takane::validate()` on a saved directory checks
that all files conform to the expected layout and types, which means a
collaborator or downstream tool can verify the integrity of a saved
object without needing to load it into R. This makes alabaster
directories suitable for deposition in data repositories where format
conformance needs to be auditable.

The Python counterpart to `alabaster` is the
[dolomite](https://github.com/ArtifactDB/dolomite-base) family of
packages, which reads and writes the same on-disk format. An object
saved with `alabaster` in R can be read with `dolomite` in Python, and
vice versa, with no conversion step.

Advantages:

- Truly language-agnostic: Python readers exist via the `dolomite`
  family of packages, enabling seamless R ↔︎ Python interoperability.
- Built on open standards (HDF5, JSON); inspectable without R.
- Versioned format specifications with independent validation via
  takane.
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

## BED format for ranges

### Writing BED files

When the object is a `GRanges` or similar ranges object and the goal is
interoperability with other tools (genome browsers, Python, command-line
utilities), exporting to BED format is often more useful than R-specific
serialization. Our `gr` has range names, a standard BED score column,
and an extra metadata column `log2fc`:

``` r

gr
```

    GRanges object with 3 ranges and 2 metadata columns:
            seqnames    ranges strand |     score    log2fc
               <Rle> <IRanges>  <Rle> | <numeric> <numeric>
      peak1     chr1   100-149      * |       500       1.2
      peak2     chr1   200-249      * |       800      -0.5
      peak3     chr1   300-349      * |       300       2.1
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

Both `rtracklayer` and `plyranges` can write BED files:

``` r

library(rtracklayer)
library(plyranges)

tmp_bed <- tempfile(fileext = ".bed")
export(gr, tmp_bed)

tmp_bed2 <- tempfile(fileext = ".bed")
write_bed(gr, tmp_bed2)
```

When reading back, names and the standard `score` column are preserved,
but `log2fc` is silently dropped — BED has no mechanism to carry
arbitrary metadata columns:

``` r

gr_rtracklayer <- import(tmp_bed)
gr_rtracklayer
```

    GRanges object with 3 ranges and 2 metadata columns:
          seqnames    ranges strand |        name     score
             <Rle> <IRanges>  <Rle> | <character> <numeric>
      [1]     chr1   100-149      * |       peak1       500
      [2]     chr1   200-249      * |       peak2       800
      [3]     chr1   300-349      * |       peak3       300
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r

gr_plyranges <- read_bed(tmp_bed2)
gr_plyranges
```

    GRanges object with 3 ranges and 2 metadata columns:
          seqnames    ranges strand |        name     score
             <Rle> <IRanges>  <Rle> | <character> <numeric>
      [1]     chr1   100-149      * |       peak1       500
      [2]     chr1   200-249      * |       peak2       800
      [3]     chr1   300-349      * |       peak3       300
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

### Saving metadata columns

To round-trip extra mcols, write them to a sidecar file. Here using
plyranges to read the BED back, then reattach from the sidecar:

``` r

tmp_meta <- tempfile(fileext = ".tsv")
write.table(
  data.frame(name = names(gr), log2fc = gr$log2fc),
  tmp_meta, sep = "\t", quote = FALSE, row.names = FALSE
)

gr_restored <- read_bed(tmp_bed2)
meta <- read.table(tmp_meta, header = TRUE, sep = "\t")
gr_restored$log2fc <- meta$log2fc
gr_restored
```

    GRanges object with 3 ranges and 3 metadata columns:
          seqnames    ranges strand |        name     score    log2fc
             <Rle> <IRanges>  <Rle> | <character> <numeric> <numeric>
      [1]     chr1   100-149      * |       peak1       500       1.2
      [2]     chr1   200-249      * |       peak2       800      -0.5
      [3]     chr1   300-349      * |       peak3       300       2.1
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

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

## Propagating object metadata

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

metadata(se_from_rds) <- fromJSON(tmp_json)
metadata(se_from_rds)
```

    $genome
    [1] "hg38"

    $pipeline
    [1] "v2.1"

    $n_samples
    [1] 4

`toJSON` handles simple R types (lists, vectors, data frames) well, but
complex objects (S4 instances, environments) need to be simplified or
omitted before serializing.

## Summary and recommendations

| Scenario | Recommended approach |
|----|----|
| Quick sharing between R users, same Bioc release | [`saveRDS()`](https://rdrr.io/pkg/BiocGenerics/man/saveRDS.html) |
| Loading an object from an older Bioc release | [`updateObject()`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html) after [`readRDS()`](https://rdrr.io/r/base/readRDS.html) |
| Large assay matrices, R-only | [`saveHDF5SummarizedExperiment()`](https://rdrr.io/pkg/HDF5Array/man/saveHDF5SummarizedExperiment.html) |
| SingleCellExperiment ↔︎ Python AnnData | [`zellkonverter::writeH5AD()`](https://rdrr.io/pkg/zellkonverter/man/writeH5AD.html) / [`readH5AD()`](https://rdrr.io/pkg/zellkonverter/man/readH5AD.html) |
| Cross-language (R + Python), general | `alabaster::saveObject()` |
| Long-term archive / data portal | `alabaster::saveObject()` |

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
     [7] HDF5Array_1.40.0            h5mread_1.4.0
     [9] rhdf5_2.56.0                DelayedArray_0.38.2
    [11] SparseArray_1.12.2          S4Arrays_1.12.0
    [13] abind_1.4-8                 Matrix_1.7-5
    [15] SummarizedExperiment_1.42.0 Biobase_2.72.0
    [17] GenomicRanges_1.64.0        Seqinfo_1.2.0
    [19] IRanges_2.46.0              S4Vectors_0.50.1
    [21] BiocGenerics_0.58.1         generics_0.1.4
    [23] MatrixGenerics_1.24.0       matrixStats_1.5.0

    loaded via a namespace (and not attached):
     [1] rjson_0.2.23             xfun_0.58                lattice_0.22-9
     [4] rhdf5filters_1.24.0      vctrs_0.7.3              tools_4.6.0
     [7] bitops_1.0-9             curl_7.1.0               parallel_4.6.0
    [10] tibble_3.3.1             pkgconfig_2.0.3          cigarillo_1.2.0
    [13] lifecycle_1.0.5          compiler_4.6.0           Rsamtools_2.28.0
    [16] Biostrings_2.80.1        codetools_0.2-20         htmltools_0.5.9
    [19] RCurl_1.98-1.19          alabaster.matrix_1.12.0  yaml_2.3.12
    [22] pillar_1.11.1            crayon_1.5.3             BiocParallel_1.46.0
    [25] tidyselect_1.2.1         digest_0.6.39            restfulr_0.0.16
    [28] fastmap_1.2.0            grid_4.6.0               cli_3.6.6
    [31] magrittr_2.0.5           XML_3.99-0.23            rmarkdown_2.31
    [34] XVector_0.52.0           httr_1.4.8               otel_0.2.0
    [37] evaluate_1.0.5           knitr_1.51               BiocIO_1.22.0
    [40] rlang_1.2.0              Rcpp_1.1.1-1.1           glue_1.8.1
    [43] alabaster.ranges_1.12.0  alabaster.schemas_1.12.0 R6_2.6.1
    [46] Rhdf5lib_2.0.0           GenomicAlignments_1.48.0
