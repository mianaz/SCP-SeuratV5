# Internal helper functions for macOS compatibility

#' Apply basic macOS fixes  
#'
#' Simple function to apply essential macOS compatibility fixes
#' 
#' @return NULL invisibly
#' @keywords internal
apply_macos_fixes <- function() {
  if (Sys.info()["sysname"] == "Darwin") {
    # Essential OpenMP fix for macOS
    Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  }
  invisible(NULL)
}

#' macOS Detection Helper
#'
#' Helper function to detect if running on macOS
#'
#' @return Boolean indicating if running on macOS
#' @keywords internal
is_macos <- function() {
  return(Sys.info()["sysname"] == "Darwin")
}

#' Apple Silicon Detection Helper
#'
#' Helper function to detect if running on Apple Silicon (M1/M2/M3)
#'
#' @return Boolean indicating if running on Apple Silicon
#' @keywords internal
is_apple_silicon <- function() {
  if (!is_macos()) {
    return(FALSE)
  }
  
  # Check via machine architecture
  arch <- Sys.info()["machine"]
  return(grepl("arm64", arch, ignore.case = TRUE))
}