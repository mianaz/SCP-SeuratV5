.onLoad <- function(libname, pkgname) {
  # Set environment variables to suppress known warnings
  # These need to be set before Python is initialized
  Sys.setenv(PYTHONWARNINGS = "ignore::DeprecationWarning:pkg_resources")

  # Note: "Registered S3 method overwritten" for FoldChange.default is intentional.
  # SCPNext provides an enhanced version with NA handling (na.rm = TRUE).
  # This message cannot be suppressed as it occurs during namespace loading.
}

.onAttach <- function(libname, pkgname) {
  # Get auto_activate option (default FALSE)
  auto_activate <- getOption("SCPNext.auto_activate", FALSE)

  # Build and show startup message
  msg <- .build_startup_message(auto_activate)
  packageStartupMessage(msg)
}

.build_startup_message <- function(auto_activate = FALSE) {
  # Color codes
  cyan <- "\033[36m"
  green <- "\033[32m"
  yellow <- "\033[33m"
  red <- "\033[31m"
  bold <- "\033[1m"
  reset <- "\033[0m"

  # ASCII Logo
  logo <- paste0(
    cyan, bold,
    "\n   _____ _____ _____  _   _           _   ",
    "\n  / ____/ ____|  __ \\| \\ | |         | |  ",
    "\n | (___| |    | |__) |  \\| | _____  _| |_ ",
    "\n  \\___ \\ |    |  ___/| . ` |/ _ \\ \\/ / __|",
    "\n  ____) | |____| |    | |\\  |  __/>  <| |_ ",
    "\n |_____/\\_____|_|    |_| \\_|\\___/_/\\_\\\\__|",
    reset, "\n"
  )

  # Get version information
  scpnext_ver <- as.character(packageVersion("SCPNext"))

  versions <- character(0)
  if (requireNamespace("Seurat", quietly = TRUE)) {
    seurat_ver <- as.character(packageVersion("Seurat"))
    versions <- c(versions, paste0("Seurat ", seurat_ver))
  }
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    so_ver <- as.character(packageVersion("SeuratObject"))
    versions <- c(versions, paste0("SeuratObject ", so_ver))
  }

  version_line <- if (length(versions) > 0) {
    paste0("  ", bold, "Version: ", reset, scpnext_ver,
           " (", paste(versions, collapse = " | "), ")")
  } else {
    paste0("  ", bold, "Version: ", reset, scpnext_ver)
  }

  # Check Seurat version and add performance note
  seurat_msg <- ""
  if (requireNamespace("Seurat", quietly = TRUE)) {
    seurat_version <- packageVersion("Seurat")
    is_v5 <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1]) >= 5

    if (!is_v5) {
      seurat_msg <- paste0(
        "\n  ", yellow, "\u26A0 Seurat v4 detected - consider upgrading to v5 for better performance", reset
      )
    }
  }

  # Check Python environment
  py_status <- .get_python_status(auto_activate)

  # Additional info
  info_msg <- paste0(
    "\n  ", bold, "Tip: ", reset,
    "Use suppressPackageStartupMessages(library(SCPNext)) for quiet loading"
  )

  # Combine all parts
  paste0(
    logo,
    version_line,
    seurat_msg,
    py_status$message,
    info_msg,
    "\n"
  )
}

.get_python_status <- function(auto_activate = FALSE) {
  # Color codes
  green <- "\033[32m"
  yellow <- "\033[33m"
  reset <- "\033[0m"
  bold <- "\033[1m"

  if (!uv_env_exists()) {
    return(list(
      status = "not_found",
      message = paste0(
        "\n  ", yellow, "\u26A0 Python environment not found", reset,
        " - Run ", bold, "PrepareEnv()", reset, " to set up"
      )
    ))
  }

  # Environment exists
  env_activated <- FALSE

  # Check if already activated
  if (reticulate::py_available()) {
    tryCatch({
      reticulate::py_config()
      env_activated <- TRUE
    }, error = function(e) {
      env_activated <<- FALSE
    })
  }

  if (env_activated) {
    return(list(
      status = "activated",
      message = paste0(
        "\n  ", green, "\u2713 Python environment active", reset
      )
    ))
  }

  # Try auto-activation if enabled
  if (auto_activate) {
    tryCatch({
      use_uv_env()
      return(list(
        status = "auto_activated",
        message = paste0(
          "\n  ", green, "\u2713 Python environment auto-activated", reset
        )
      ))
    }, error = function(e) {
      return(list(
        status = "auto_activate_failed",
        message = paste0(
          "\n  ", yellow, "\u26A0 Python environment found but activation failed", reset,
          " - Run ", bold, "use_uv_env()", reset, " manually"
        )
      ))
    })
  } else {
    return(list(
      status = "needs_activation",
      message = paste0(
        "\n  ", green, "\u2713 Python environment ready", reset,
        " - Run ", bold, "use_uv_env()", reset, " to activate"
      )
    ))
  }
}

