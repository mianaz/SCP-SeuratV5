#' Cache Management for External Data Files
#'
#' Internal functions for managing cached external data files to reduce
#' network dependencies during package checks and testing.
#'
#' @name SCP-cache
#' @keywords internal

#' Get cache directory path
#' @keywords internal
.get_cache_dir <- function() {
  # Use R.cache if available, otherwise use package-specific cache
  cache_dir <- tryCatch({
    R.cache::getCacheRootPath()
  }, error = function(e) {
    file.path(tempdir(), "SCP_cache")
  })

  # Create subdirectory for SCP
  scp_cache <- file.path(cache_dir, "SCP", "external_data")
  if (!dir.exists(scp_cache)) {
    dir.create(scp_cache, recursive = TRUE, showWarnings = FALSE)
  }
  return(scp_cache)
}

#' Check if cached file exists and is recent
#' @param cache_file Path to cached file
#' @param max_age Maximum age in days (default: 30)
#' @keywords internal
.is_cache_valid <- function(cache_file, max_age = 30) {
  if (!file.exists(cache_file)) return(FALSE)

  file_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
  return(as.numeric(file_age) < max_age && file.size(cache_file) > 0)
}

#' Download file with caching and SSL fix
#'
#' Downloads a file with automatic caching, SSL certificate workaround,
#' and fallback mechanisms.
#'
#' @param url URL to download from
#' @param destfile Destination file path
#' @param mode Download mode ("w" or "wb")
#' @param use_cache Whether to use caching (default: TRUE)
#' @param max_cache_age Maximum cache age in days (default: 30)
#' @param verbose Print messages (default: FALSE)
#' @return Logical indicating success
#' @keywords internal
.download_with_cache <- function(url,
                                destfile,
                                mode = ifelse(.Platform$OS.type == "windows", "wb", "w"),
                                use_cache = TRUE,
                                max_cache_age = 30,
                                verbose = FALSE) {

  # Generate cache filename based on URL
  cache_name <- paste0(
    digest::digest(url, algo = "md5"),
    "_",
    basename(url)
  )
  cache_file <- file.path(.get_cache_dir(), cache_name)

  # Check cache first
  if (use_cache && .is_cache_valid(cache_file, max_cache_age)) {
    if (verbose) message("Using cached file: ", basename(cache_file))
    file.copy(cache_file, destfile, overwrite = TRUE)
    return(TRUE)
  }

  # Download with SSL workaround
  success <- FALSE

  # Method 1: curl with SSL bypass
  if (!success) {
    success <- tryCatch({
      if (verbose) message("Downloading with curl (SSL bypass)...")
      download.file(url = url,
                   destfile = destfile,
                   mode = mode,
                   method = "curl",
                   extra = "-k",
                   quiet = !verbose)
      file.exists(destfile) && file.size(destfile) > 0
    }, error = function(e) {
      if (verbose) message("Curl method failed: ", e$message)
      FALSE
    })
  }

  # Method 2: Default method
  if (!success) {
    success <- tryCatch({
      if (verbose) message("Trying default download method...")
      download.file(url = url,
                   destfile = destfile,
                   mode = mode,
                   quiet = !verbose)
      file.exists(destfile) && file.size(destfile) > 0
    }, error = function(e) {
      if (verbose) message("Default method failed: ", e$message)
      FALSE
    })
  }

  # Method 3: httr with SSL bypass (if available)
  if (!success && requireNamespace("httr", quietly = TRUE)) {
    success <- tryCatch({
      if (verbose) message("Trying httr with SSL bypass...")
      response <- httr::GET(url, httr::config(ssl_verifypeer = 0))
      if (httr::status_code(response) == 200) {
        writeBin(httr::content(response, "raw"), destfile)
        TRUE
      } else {
        FALSE
      }
    }, error = function(e) {
      if (verbose) message("httr method failed: ", e$message)
      FALSE
    })
  }

  # If successful, update cache
  if (success && use_cache) {
    if (verbose) message("Updating cache...")
    file.copy(destfile, cache_file, overwrite = TRUE)
  }

  return(success)
}

#' Clear SCP cache
#'
#' Remove cached external data files
#'
#' @param older_than Remove only files older than this many days (NULL = all files)
#' @export
clear_SCP_cache <- function(older_than = NULL) {
  cache_dir <- .get_cache_dir()

  if (!dir.exists(cache_dir)) {
    message("No cache directory found.")
    return(invisible(NULL))
  }

  cache_files <- list.files(cache_dir, full.names = TRUE)

  if (length(cache_files) == 0) {
    message("Cache is already empty.")
    return(invisible(NULL))
  }

  if (!is.null(older_than)) {
    # Filter by age
    file_ages <- sapply(cache_files, function(f) {
      as.numeric(difftime(Sys.time(), file.mtime(f), units = "days"))
    })
    cache_files <- cache_files[file_ages > older_than]
  }

  if (length(cache_files) > 0) {
    unlink(cache_files)
    message("Removed ", length(cache_files), " cached file(s).")
  } else {
    message("No files matched criteria.")
  }

  invisible(NULL)
}

#' Show SCP cache information
#'
#' Display information about cached files
#'
#' @export
show_SCP_cache <- function() {
  cache_dir <- .get_cache_dir()

  if (!dir.exists(cache_dir)) {
    message("No cache directory found.")
    return(invisible(NULL))
  }

  cache_files <- list.files(cache_dir, full.names = TRUE)

  if (length(cache_files) == 0) {
    message("Cache is empty.")
    return(invisible(NULL))
  }

  cache_info <- data.frame(
    file = basename(cache_files),
    size_mb = round(file.size(cache_files) / 1024^2, 2),
    age_days = round(as.numeric(sapply(cache_files, function(f) {
      difftime(Sys.time(), file.mtime(f), units = "days")
    })), 1),
    stringsAsFactors = FALSE
  )

  message("SCP Cache Directory: ", cache_dir)
  message("Total size: ", round(sum(cache_info$size_mb), 2), " MB")
  message("\nCached files:")
  print(cache_info, row.names = FALSE)

  invisible(cache_info)
}