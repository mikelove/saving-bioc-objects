# Contributing to savingBiocObjects

Thank you for your interest in improving this documentation. Contributions of
all kinds are welcome: fixing typos, improving explanations, adding new
sections, or updating code examples as the Bioconductor ecosystem evolves.

## How to contribute

1. **Open an issue first** for anything beyond a minor fix — describe what you
   think is missing or wrong, and we can agree on the approach before you
   write.

2. **Fork the repository** and create a branch from `main`.

3. **Edit the vignette** at `vignettes/savingBiocObjects.qmd`. The file is
   standard Quarto markdown; any code chunks should be reproducible.

4. **Preview locally** with pkgdown:
   ```r
   pkgdown::build_site()
   # or, for just the vignette:
   pkgdown::build_article("savingBiocObjects")
   ```

5. **Open a pull request** against `main`. The GitHub Actions workflow will
   render the site automatically; the deployed preview is only published on
   merge to `main`.

## Content guidelines

- **Be concrete.** Prefer runnable code examples over prose-only descriptions.
- **State trade-offs.** Each approach should include honest disadvantages, not
  just advantages.
- **Keep recommendations up to date.** If a new Bioconductor package supersedes
  an older approach, note it and update the summary table.
- **Cite packages properly.** Use `citation("pkgname")` and link to the
  Bioconductor landing page.

## Code of conduct

This project follows the
[Bioconductor Code of Conduct](https://bioconductor.org/about/code-of-conduct/).
Please be respectful and constructive in all interactions.

## Questions

Open a [GitHub issue](https://github.com/mikelove/saving-bioc-objects/issues)
or reach out via the [Bioconductor support site](https://support.bioconductor.org).
