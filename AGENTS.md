# Repository Guidelines

## Project Structure & Module Organization
- `R/` holds exported workflow modules (`SCP-analysis.R`, `SCP-cellqc.R`) plus Seurat V5 helpers under `R/fixes/`.
- `src/` contains the Rcpp sources; artifacts rebuild automatically during `devtools::install()`.
- `inst/python/` ships optional reticulate hooks (`SCP_analysis.py`) for advanced features.
- `data/` and `man/` store example objects and generated documentation; avoid manual edits.
- `tests/testthat/` provides edition-3 specs and snapshots under `_snaps/`, while `docs/` is the pkgdown output.

## Build, Test, and Development Commands
- `R -q -e "devtools::load_all()"` loads package code for interactive iteration.
- `R -q -e "devtools::document(); devtools::install()"` refreshes roxygen docs and reinstalls locally.
- `R -q -e "devtools::check()"` runs linting, tests, and metadata checks; execute before opening a PR.
- `R -q -e "testthat::test_file('tests/testthat/test-workflow-comprehensive.R')"` targets the workflow regression suite.
- `R -q -e "pkgdown::build_site()"` regenerates the reference site inside `docs/`.

## Coding Style & Naming Conventions
- Follow tidyverse-aligned R style: 2-space indents, spaces around `=` in named args, and `snake_case` helpers like `get_seurat_data`.
- Exported user APIs stay UpperCamelCase (`RunKNNPredict`, `AnnotateFeatures`); keep filenames aligned with exported objects.
- Maintain concise roxygen2 blocks that call out Seurat layer usage and compatibility notes.
- Run `styler::style_pkg()` before committing R changes; never hand-edit generated `.Rd` files.

## Testing Guidelines
- Tests rely on `testthat` (edition 3) with heavy snapshot coverage in `_snaps/`; review updates when `devtools::test()` rewrites snapshots.
- Name new specs `test-<feature>.R` and scope helper functions within the test file.
- For Seurat I/O changes, add paired V4/V5 fixtures or guarded skips to preserve backwards support.
- Extend `test-workflow-comprehensive.R` when altering end-to-end behaviour or integration branches.

## Commit & Pull Request Guidelines
- Write imperative commit subjects under 72 characters; prefixes like `feat:` and `fix:` mirror recent history.
- Reference issues with `#<id>` and describe Seurat-impacting behaviour or API changes in the body.
- PRs should list validation commands run (`devtools::check()`, targeted tests) and add screenshots for plotting or pkgdown updates.
- Confirm optional dependencies (Python, Harmony, etc.) used or skipped so reviewers can reproduce your setup.
