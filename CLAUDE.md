# SCP-SeuratV5 Development Guide

## Build/Test/Lint Commands
- Build: `devtools::build()` or `R CMD build .`
- Documentation: `devtools::document()` 
- Check: `devtools::check()` or `R CMD check --no-manual .`
- Test: `devtools::test()` or `testthat::test_package("SCP")`
- Test single file: `testthat::test_file("tests/testthat/test-file.R")`
- Lint: `styler::style_pkg()`
- Build docs site: `pkgdown::build_site()`

## Code Style Guidelines
- Function names: CamelCase with initial capital (e.g., `RunPAGA`)
- Helper functions: snake_case (e.g., `try_get`)
- Variables: snake_case
- Parameters: period-separated (e.g., `group.by`, `reduction`)
- Indentation: 2 spaces
- Documentation: roxygen2 with comprehensive @param, @return, @examples
- Comments: `#` with space following
- Imports: Use `@importFrom` tags or `package::function` notation
- Error handling: Use `tryCatch()` for robust error handling
- Line length: Keep under 120 characters
- Default parameters: Provide sensible defaults when possible