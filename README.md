# savingBiocObjects

Documentation and guidance on strategies for saving and sharing Bioconductor
objects.

Topics covered:

- `saveRDS()` / `readRDS()` for simple R-to-R sharing, and cross-release stability
- HDF5-backed storage with `HDF5Array` for large assay data
- The `alabaster` ecosystem for language-agnostic serialization
- Converting between `SingleCellExperiment` and `AnnData` with `zellkonverter`
- Exporting genomic ranges to BED with `rtracklayer` or `plyranges`
- Preserving `metadata(object)` as a JSON sidecar file
- A summary table of recommendations by use case

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) for
guidelines on adding new content, fixing errors, or proposing new sections.

## License

MIT
