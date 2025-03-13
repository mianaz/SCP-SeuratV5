#' This function prepares the SCP Python environment by installing the required dependencies and setting up the environment.
#'
#' @param miniconda_repo Repositories for miniconda. Default is \code{https://repo.anaconda.com/miniconda}
#' @param force Whether to force a new environment to be created. If \code{TRUE}, the existing environment will be recreated. Default is \code{FALSE}
#' @param version A character vector specifying the version of the environment (default is "3.10-1").
#' @param use_conda_install Whether to use conda instead of pip to install Python packages. This can be faster but may have less package compatibility. Default is \code{FALSE}
#' @param prompt Whether to prompt for confirmation before creating or recreating environments. Default is \code{TRUE}
#' @inheritParams check_Python
#' @details This function prepares the SCP Python environment by checking if conda is installed, creating a new conda environment if needed, installing the required packages, and setting up the Python environment for use with SCP.
#'
#' Prerequisites:
#' \itemize{
#'   \item Internet connection to download conda packages
#'   \item Sufficient disk space (approximately 3-4GB)
#'   \item Windows users: Windows 7 or newer
#'   \item Mac users: macOS 10.13 or newer
#'   \item Linux users: glibc 2.17 or newer
#' }
#'
#' Process:
#' \enumerate{
#'   \item The function first checks if conda is installed, installing miniconda if necessary
#'   \item It then checks for an existing SCP environment (default name: "SCP_env")
#'   \item If the environment exists and force=FALSE, it will use the existing environment
#'   \item If force=TRUE, it will recreate the environment
#'   \item Required Python packages are installed as specified by the version parameter
#' }
#'
#' In order to create the environment, this function requires the path to the conda binary. If \code{conda} is set to \code{"auto"}, it will attempt to automatically find the conda binary.
#' If a conda environment with the specified name already exists and \code{force} is set to \code{FALSE}, the function will use the existing environment. If \code{force} set to \code{TRUE}, the existing environment will be recreated. Note that recreating the environment will remove any existing data in the environment.
#' The function also checks if the package versions in the environment meet the requirements specified by the \code{version} parameter. The default is \code{3.10-1}.
#'
#' @return A list containing information about the created/used environment, including:
#' \itemize{
#'   \item conda_path: Path to the conda executable
#'   \item env_path: Path to the conda environment
#'   \item python_path: Path to the Python executable in the environment
#'   \item python_version: The Python version used
#'   \item packages: The packages installed
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic usage: create or use existing SCP environment
#' PrepareEnv()
#' 
#' # Force recreation of the environment
#' PrepareEnv(force = TRUE)
#' 
#' # Use a specific Python version
#' PrepareEnv(version = "3.11-1")
#' 
#' # Use conda instead of pip for faster package installation
#' PrepareEnv(use_conda_install = TRUE)
#' }
#' 
#' @export
PrepareEnv <- function(conda = "auto", miniconda_repo = "https://repo.anaconda.com/miniconda",
                       envname = NULL, version = "3.10-1", force = FALSE, 
                       use_conda_install = FALSE, prompt = TRUE, ...) {
  # Check for adequate disk space (rough estimate)
  if (!is.null(getwd())) {
    free_space <- tryCatch({
      if (.Platform$OS.type == "windows") {
        wmic_cmd <- system2("wmic", args = c("logicaldisk", "where", paste0("DeviceID='", substr(getwd(), 1, 2), "'"), "get", "freespace"), stdout = TRUE)
        as.numeric(wmic_cmd[2]) / (1024^3)  # Convert to GB
      } else if (.Platform$OS.type == "unix") {
        df_cmd <- system2("df", args = c("-k", getwd()), stdout = TRUE)
        df_parts <- strsplit(df_cmd[2], "\\s+")[[1]]
        as.numeric(df_parts[4]) / (1024^2)  # Convert to GB
      } else {
        NA
      }
    }, error = function(e) NA)
    
    if (!is.na(free_space) && free_space < 5) {
      warning("Low disk space detected (", round(free_space, 1), " GB free). SCP environment requires ~3-4GB of space.", 
              " Installation may fail if disk space is insufficient.")
    }
  }
  
  envname <- get_envname(envname)
  
  message("Preparing SCP Python environment...")
  message("- Environment name: ", envname)
  
  requirements <- Env_requirements(version = version)
  python_version <- requirements[["python"]]
  packages <- requirements[["packages"]]
  
  message("- Python version: ", python_version)
  
  # Find conda binary
  if (!is.null(conda)) {
    if (identical(conda, "auto")) {
      conda <- find_conda()
      if (!is.null(conda)) {
        message("- Found conda at: ", conda)
      }
    } else {
      options(reticulate.conda_binary = conda)
      conda <- find_conda()
      message("- Using specified conda at: ", conda)
    }
  }
  
  # Check if conda is available, otherwise env is FALSE
  if (is.null(conda)) {
    message("- Conda not found")
    env <- FALSE
  } else {
    # Check if environment exists
    envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    env <- env_exist(conda = conda, envname = envname, envs_dir = envs_dir)
    
    if (isTRUE(env)) {
      env_path <- paste0(envs_dir, "/", envname)
      message("- Found existing SCP environment at: ", env_path)
      
      # Check Python version in existing environment
      python_path <- conda_python(conda = conda, envname = envname)
      installed_python_version <- tryCatch({
        reticulate:::python_version(python_path)
      }, error = function(e) {
        numeric_version("0.0.0")
      })
      
      python_version_numeric <- numeric_version(python_version)
      
      # Version compatibility check
      compatible <- installed_python_version >= numeric_version("3.9.0") && 
                   installed_python_version < numeric_version("3.13.0")
      
      if (!compatible) {
        message("⚠️ WARNING: Existing environment has Python ", as.character(installed_python_version), 
                " but SCP requires Python 3.9-3.12")
        message("⚠️ Consider recreating the environment with force=TRUE")
        
        if (prompt && interactive()) {
          recreate <- utils::menu(c("Yes", "No"), 
                             title = paste0("Would you like to recreate the environment with Python ", python_version, "?"))
          if (recreate == 1) {
            force <- TRUE
            message("- Will recreate environment with Python ", python_version)
          } else {
            message("- Using existing environment with incompatible Python version")
          }
        }
      } else if (as.character(installed_python_version) != python_version) {
        message("- Existing environment has Python ", as.character(installed_python_version), 
                ", requested version is ", python_version)
        if (!force && prompt && interactive()) {
          recreate <- utils::menu(c("Yes", "No"), 
                            title = paste0("Would you like to recreate the environment with Python ", python_version, "?"))
          if (recreate == 1) {
            force <- TRUE
            message("- Will recreate environment with Python ", python_version)
          } else {
            message("- Using existing environment with Python ", as.character(installed_python_version))
          }
        }
      }
      
      # If force is TRUE, delete the existing environment
      if (isTRUE(force) && isTRUE(env)) {
        message("- Removing existing environment...")
        unlink(env_path, recursive = TRUE)
        env <- FALSE
      }
    } else {
      message("- No existing SCP environment found")
    }
  }
  
  # If environment doesn't exist or was deleted, create it
  if (isFALSE(env)) {
    force <- TRUE
    
    # If conda is not found, install miniconda
    if (is.null(conda)) {
      message("- Installing miniconda...")
      
      # Increase timeout for large downloads
      options(timeout = 360)
      version <- "3"
      info <- as.list(Sys.info())
      
      # Special handling for Apple Silicon
      if (info$sysname == "Darwin" && info$machine == "arm64") {
        base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
        name <- "Miniforge3-MacOSX-arm64.sh"
        url <- file.path(base, name)
        message("- Detected Apple Silicon, using Miniforge3")
      } else {
        # Standard miniconda setup for other platforms
        base <- miniconda_repo
        info <- as.list(Sys.info())
        arch <- reticulate:::miniconda_installer_arch(info)
        version <- as.character(version)
        
        name <- if (reticulate:::is_windows()) {
          sprintf("Miniconda%s-latest-Windows-%s.exe", version, arch)
        } else if (reticulate:::is_osx()) {
          sprintf("Miniconda%s-latest-MacOSX-%s.sh", version, arch)
        } else if (reticulate:::is_linux()) {
          sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
        } else {
          stop("Unsupported platform: ", shQuote(Sys.info()[["sysname"]]))
        }
        url <- file.path(base, name)
      }
      
      message("- Downloading from: ", url)
      options(reticulate.miniconda.url = url)
      
      # Handle USER variable in path
      if (!is.na(Sys.getenv("USER", unset = NA))) {
        miniconda_path <- gsub(pattern = "\\$USER", replacement = Sys.getenv("USER"), reticulate::miniconda_path())
      } else {
        miniconda_path <- reticulate::miniconda_path()
      }
      
      # Remove existing miniconda if present
      unlink(miniconda_path, recursive = TRUE)
      
      # Install miniconda
      message("- Installing miniconda to: ", miniconda_path)
      reticulate::install_miniconda(path = miniconda_path, force = TRUE, update = FALSE)
      conda <- reticulate:::miniconda_conda(miniconda_path)
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
      message("- Miniconda installation complete")
    }
    
    # Validate Python version compatibility
    if (python_version < numeric_version("3.9.0") || python_version >= numeric_version("3.13.0")) {
      stop("SCP currently only supports Python version 3.9-3.12!")
    }
    
    message("- Creating new conda environment '", envname, "' with Python ", python_version, "...")
    
    # Create conda environment with Python
    python_path <- reticulate::conda_create(
      conda = conda, 
      envname = envname, 
      python_version = python_version, 
      packages = "pytables"
    )
    
    env_path <- paste0(envs_dir, "/", envname)
    env <- file.exists(env_path)
    
    if (isFALSE(env)) {
      message("- Failed to create conda environment. Debug information:")
      print(reticulate:::conda_info(conda = conda))
      print(reticulate::conda_list(conda = conda))
      stop(
        "Unable to find SCP environment under the expected path: ", env_path, "\n",
        "conda: ", conda, "\n",
        "SCP python: ", python_path, "\n"
      )
    }
    
    message("- Environment created successfully at: ", env_path)
  }
  
  # Install packages
  message("- Installing required Python packages...")
  if (use_conda_install) {
    message("- Using conda for package installation (faster but may have compatibility issues)")
  } else {
    message("- Using pip for package installation (default, better compatibility)")
  }
  
  check_result <- check_Python(
    packages = packages, 
    envname = envname, 
    conda = conda, 
    force = force, 
    pip = !use_conda_install,
    ...
  )
  
  # Configure the Python environment
  message("- Configuring Python environment...")
  Sys.unsetenv("RETICULATE_PYTHON")
  python_path <- conda_python(conda = conda, envname = envname)
  
  # Check if reticulate is already loaded, if not, use conda env before loading
  if (!isNamespaceLoaded("reticulate")) {
    # Set the Python path before loading reticulate
    Sys.setenv(RETICULATE_PYTHON = python_path)
    Sys.setenv(RETICULATE_PYTHON_ENV = envname)
    # Now load reticulate with the correct Python path
    requireNamespace("reticulate", quietly = TRUE)
  }
  
  # Use the specified Python regardless
  reticulate::use_python(python_path, required = TRUE)
  
  # Display environment information
  pyinfo <- utils::capture.output(reticulate::py_config())
  pyinfo_mesg <- c(
    "====================== SCP conda environment ======================",
    paste0("conda: ", conda),
    paste0("environment: ", paste0(envs_dir, "/", get_envname())),
    "======================== SCP python config ========================",
    pyinfo,
    "==================================================================="
  )
  invisible(lapply(pyinfo_mesg, packageStartupMessage))
  
  # Test importing key packages
  message("- Testing Python module imports...")
  import_results <- list()
  for (module in c("matplotlib", "scanpy", "numpy", "pandas")) {
    import_results[[module]] <- tryCatch({
      if (module == "matplotlib") {
        invisible(run_Python(command = "import matplotlib", envir = .GlobalEnv))
        if (!interactive()) {
          invisible(run_Python(command = "matplotlib.use('pdf')", envir = .GlobalEnv))
        }
        invisible(run_Python(command = "import matplotlib.pyplot as plt", envir = .GlobalEnv))
      } else {
        invisible(run_Python(command = paste0("import ", module), envir = .GlobalEnv))
      }
      TRUE
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })
  }
  
  # Check if imports were successful
  all_imports_ok <- all(sapply(import_results, function(x) isTRUE(x)))
  
  if (all_imports_ok) {
    message("✓ Successfully loaded all Python dependencies.")
  } else {
    message("⚠️ Warning: Some Python modules could not be imported:")
    for (module in names(import_results)) {
      if (!isTRUE(import_results[[module]])) {
        message("  - ", module, ": ", import_results[[module]]$error)
      }
    }
    message("⚠️ Some Python-dependent functionality may not work correctly.")
  }
  
  # Return information about the environment
  message("\n✓ SCP Python environment setup complete!")
  message("- Environment location: ", env_path)
  message("- Python path: ", python_path)
  message("- Python version: ", reticulate:::python_version(python_path))
  message("\nNext steps:")
  message("1. Test functionality with: run_Python(\"import scanpy; scanpy.__version__\")")
  message("2. Run SCP functions that depend on Python, such as RunPAGA() or RunSCVELO()")
  message("3. If you encounter issues, try reinstalling with force=TRUE")
  message("4. To clean up this environment later, use CleanupEnv()")
  
  invisible(list(
    conda_path = conda,
    env_path = env_path,
    python_path = python_path,
    python_version = reticulate:::python_version(python_path),
    packages = packages
  ))
}

#' Clean up the SCP Python environment
#'
#' This function removes the conda environment created by PrepareEnv
#' to free up disk space when it's no longer needed.
#'
#' @param envname Name of the conda environment to remove. Default is NULL, which uses the 
#'               environment name set in options("SCP_env_name") or "SCP_env".
#' @param conda Path to the conda executable. Default is "auto" which will attempt to 
#'             automatically find conda.
#' @param prompt Whether to prompt for confirmation before removing the environment. Default is TRUE.
#'
#' @return Invisible logical indicating whether the environment was successfully removed.
#'
#' @examples
#' \dontrun{
#' # Remove the default SCP environment
#' CleanupEnv()
#' 
#' # Remove a specific environment without prompting
#' CleanupEnv(envname = "SCP_test_env", prompt = FALSE)
#' }
#'
#' @export
CleanupEnv <- function(envname = NULL, conda = "auto", prompt = TRUE) {
  envname <- get_envname(envname)
  
  # Find conda binary
  if (!is.null(conda)) {
    if (identical(conda, "auto")) {
      conda <- find_conda()
    } else {
      options(reticulate.conda_binary = conda)
      conda <- find_conda()
    }
  }
  
  # Check if conda is available
  if (is.null(conda)) {
    message("Conda not found. Cannot remove environment.")
    return(invisible(FALSE))
  }
  
  # Check if environment exists
  envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
  env <- env_exist(conda = conda, envname = envname, envs_dir = envs_dir)
  
  if (!isTRUE(env)) {
    message("Environment '", envname, "' not found. Nothing to remove.")
    return(invisible(FALSE))
  }
  
  env_path <- paste0(envs_dir, "/", envname)
  
  # Get environment information
  env_info <- tryCatch({
    python_path <- conda_python(conda = conda, envname = envname)
    python_version <- reticulate:::python_version(python_path)
    list(path = env_path, python = python_path, version = python_version)
  }, error = function(e) {
    list(path = env_path, python = NA, version = NA)
  })
  
  # Display environment information and prompt for confirmation
  message("About to remove SCP Python environment:")
  message("- Environment name: ", envname)
  message("- Environment path: ", env_info$path)
  if (!is.na(env_info$version)) {
    message("- Python version: ", env_info$version)
  }
  message("- Disk space to be freed: ", format_size(dir_size(env_info$path)))
  
  proceed <- TRUE
  if (prompt && interactive()) {
    proceed <- utils::menu(c("Yes", "No"), 
                       title = "Are you sure you want to remove this environment?") == 1
  }
  
  if (proceed) {
    message("Removing environment '", envname, "'...")
    
    # Try using conda remove first (cleaner)
    success <- tryCatch({
      result <- system2(conda, c("env", "remove", "-n", envname, "--yes"), stdout = TRUE, stderr = TRUE)
      # Check if the command was successful
      !any(grepl("error", result, ignore.case = TRUE))
    }, error = function(e) FALSE)
    
    # If conda remove failed, try directly removing the directory
    if (!success) {
      message("Conda removal failed, attempting direct directory removal...")
      success <- tryCatch({
        unlink(env_path, recursive = TRUE)
        !file.exists(env_path)
      }, error = function(e) FALSE)
    }
    
    if (success) {
      message("✓ Successfully removed environment '", envname, "'")
      # Unset Python environment variables to avoid pointing to deleted environment
      Sys.unsetenv("RETICULATE_PYTHON")
      return(invisible(TRUE))
    } else {
      message("⚠️ Failed to remove environment '", envname, "'")
      message("You may need to manually remove the directory: ", env_path)
      return(invisible(FALSE))
    }
  } else {
    message("Environment removal cancelled.")
    return(invisible(FALSE))
  }
}

#' Get the size of a directory in bytes
#' 
#' @param path Directory path
#' @return Directory size in bytes, or 0 if directory doesn't exist
#' @keywords internal
dir_size <- function(path) {
  if (!dir.exists(path)) return(0)
  
  # For Windows, use PowerShell
  if (.Platform$OS.type == "windows") {
    cmd <- paste0('(Get-ChildItem -Path "', path, '" -Recurse | Measure-Object -Property Length -Sum).Sum')
    size <- as.numeric(system2("powershell", c("-Command", cmd), stdout = TRUE, stderr = FALSE))
    return(if (is.na(size)) 0 else size)
  } 
  
  # For Unix-like systems, use du
  if (.Platform$OS.type == "unix") {
    result <- system2("du", c("-sb", path), stdout = TRUE, stderr = FALSE)
    # Parse the result (format: "<size> <path>")
    size <- as.numeric(strsplit(result, "\\s+")[[1]][1])
    return(if (is.na(size)) 0 else size)
  }
  
  # Fallback method if system commands don't work
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  sum(file.size(files), na.rm = TRUE)
}

#' Format a size in bytes to a human-readable format
#' 
#' @param bytes Size in bytes
#' @return Human-readable string (e.g., "1.2 GB")
#' @keywords internal
format_size <- function(bytes) {
  if (bytes < 1024) return(paste0(round(bytes, 1), " B"))
  if (bytes < 1024^2) return(paste0(round(bytes/1024, 1), " KB"))
  if (bytes < 1024^3) return(paste0(round(bytes/1024^2, 1), " MB"))
  return(paste0(round(bytes/1024^3, 1), " GB"))
}

#' Env_requirements function
#'
#' This function provides the SCP Python environment requirements for a specific version.
#'
#' @param version A character vector specifying the version of the environment (default is "3.10-1").
#'                Valid options are "3.9-1", "3.10-1", "3.11-1", and "3.12-1", representing
#'                different Python versions.
#' @return A list containing:
#'         \item{python}{The Python version to use (e.g., "3.10")}
#'         \item{packages}{A character vector of required package names}
#' @details All supported Python versions use the same set of required packages.
#'          The function simply maps the version parameter to the appropriate Python version
#'          and returns a standardized list of package requirements.
#' @examples
#' # Get requirements for version "3.10-1"
#' Env_requirements("3.10-1")
#'
#' @export
Env_requirements <- function(version = "3.10-1") {
  version <- match.arg(version, choices = c("3.9-1", "3.10-1", "3.11-1", "3.12-1"))
  
  # Extract Python version from the parameter
  python_version <- switch(version,
    "3.9-1" = "3.9",
    "3.10-1" = "3.10",
    "3.11-1" = "3.11",
    "3.12-1" = "3.12"
  )
  
  # Define the common packages needed for all Python versions
  common_packages <- c(
    "leidenalg" = "leidenalg",
    "matplotlib" = "matplotlib",
    "numba" = "numba",
    "numpy" = "numpy",
    "palantir" = "palantir",
    "pandas" = "pandas",
    "python-igraph" = "python-igraph",
    "scanpy" = "scanpy",
    "scikit-learn" = "scikit-learn",
    "scipy" = "scipy",
    "scvelo" = "scvelo",
    "wot" = "wot",
    "trimap" = "trimap",
    "pacmap" = "pacmap",
    "phate" = "phate",
    "bbknn" = "bbknn",
    "scanorama" = "scanorama",
    "scvi-tools" = "scvi-tools"
  )
  
  # Create and return the requirements list
  requirements <- list(
    python = python_version,
    packages = common_packages
  )
  
  return(requirements)
}

#' Show all the python packages in the environment
#'
#' @inheritParams check_Python
#' @export
installed_Python_pkgs <- function(envname = NULL, conda = "auto") {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- reticulate:::conda_list_packages(conda = conda, envname = envname, no_pip = FALSE)
  return(all_installed)
}

#' Check if the python package exists in the environment
#'
#' @inheritParams check_Python
#' @export
exist_Python_pkgs <- function(packages, envname = NULL, conda = "auto") {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- installed_Python_pkgs(envname = envname, conda = conda)
  packages_installed <- NULL
  for (i in seq_along(packages)) {
    pkg <- packages[i]
    if (grepl("==", pkg)) {
      pkg_info <- strsplit(pkg, split = "==")[[1]]
      pkg_name <- names(pkg) %||% pkg_info[1]
      pkg_version <- pkg_info[2]
    } else if (grepl("git+", pkg)) {
      pkg_info <- strsplit(pkg, "/")[[1]]
      pkg_name <- names(pkg) %||% pkg_info[length(pkg_info)]
      pkg_version <- NA
    } else {
      pkg_name <- names(pkg) %||% pkg
      pkg_version <- NA
    }
    if (pkg_name %in% all_installed$package) {
      if (!is.na(pkg_version)) {
        packages_installed[pkg] <- all_installed$version[all_installed$package == pkg_name] == pkg_version
      } else {
        packages_installed[pkg] <- TRUE
      }
    } else {
      packages_installed[pkg] <- FALSE
    }
  }
  return(packages_installed)
}

#' Check if a conda environment exists
#'
#' @param envs_dir Directories in which conda environments are located.
#' @inheritParams check_Python
env_exist <- function(conda = "auto", envname = NULL, envs_dir = NULL) {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  if (!is.null(conda)) {
    if (is.null(envs_dir)) {
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    }
    exist <- file.exists(paste0(envs_dir, "/", envname))
  } else {
    exist <- FALSE
  }
  return(exist)
}

get_envname <- function(envname = NULL) {
  if (is.character(envname)) {
    envname <- envname
  } else {
    envname <- getOption("SCP_env_name", default = "SCP_env")
  }
  return(envname)
}

#' Find an appropriate conda binary
#'
#' @export
find_conda <- function() {
  conda <- tryCatch(reticulate::conda_binary(conda = "auto"), error = identity)
  conda_exist <- !inherits(conda, "error")
  if (isFALSE(conda_exist)) {
    if (!is.na(Sys.getenv("USER", unset = NA))) {
      miniconda_path <- gsub(pattern = "\\$USER", replacement = Sys.getenv("USER"), reticulate::miniconda_path())
    } else {
      miniconda_path <- reticulate::miniconda_path()
    }
    conda_exist <- reticulate:::miniconda_exists(miniconda_path) && reticulate:::miniconda_test(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- reticulate:::miniconda_conda(miniconda_path)
    } else {
      conda <- NULL
    }
  }
  return(conda)
}

#' Installs a list of packages into a specified conda environment
#'
#' @inheritParams reticulate::conda_install
conda_install <- function(envname = NULL, packages, forge = TRUE, channel = character(),
                          pip = FALSE, pip_options = character(), pip_ignore_installed = FALSE,
                          conda = "auto", python_version = NULL, ...) {
  envname <- get_envname(envname)
  reticulate:::check_forbidden_install("Python packages")
  if (missing(packages)) {
    if (!is.null(envname)) {
      fmt <- paste("argument \"packages\" is missing, with no default",
        "- did you mean 'conda_install(<envname>, %1$s)'?",
        "- use 'py_install(%1$s)' to install into the active Python environment",
        sep = "\n"
      )
      stopf(fmt, deparse1(substitute(envname)), call. = FALSE)
    } else {
      packages
    }
  }
  conda <- reticulate::conda_binary(conda)
  envname <- reticulate:::condaenv_resolve(envname)
  python_package <- if (is.null(python_version)) {
    NULL
  } else if (grepl("[><=]", python_version)) {
    paste0("python", python_version)
  } else {
    sprintf("python=%s", python_version)
  }
  python <- tryCatch(conda_python(envname = envname, conda = conda), error = identity)
  if (inherits(python, "error") || !file.exists(python)) {
    reticulate::conda_create(envname = envname, packages = python_package %||% "python", forge = forge, channel = channel, conda = conda)
    python <- conda_python(envname = envname, conda = conda)
  }
  if (!is.null(python_version)) {
    args <- reticulate:::conda_args("install", envname, python_package)
    status <- reticulate:::system2t(conda, shQuote(args))
    if (status != 0L) {
      fmt <- "installation of '%s' into environment '%s' failed [error code %i]"
      msg <- sprintf(fmt, python_package, envname, status)
      stop(msg, call. = FALSE)
    }
  }
  if (pip) {
    # target_dir <- system2(command = python, args = c("-c \"import site; print(site.getsitepackages()[0])\""), stdout = TRUE)
    # pip_options <- c(pip_options, paste("--target", target_dir))
    result <- reticulate:::pip_install(
      python = python, packages = packages,
      pip_options = pip_options, ignore_installed = pip_ignore_installed,
      conda = conda, envname = envname
    )
    return(result)
  }
  args <- reticulate:::conda_args("install", envname)
  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }
  for (ch in channels) args <- c(args, "-c", ch)
  args <- c(args, python_package, packages)
  result <- reticulate:::system2t(conda, shQuote(args))
  if (result != 0L) {
    fmt <- "one or more Python packages failed to install [error code %i]"
    stopf(fmt, result)
  }
  invisible(packages)
}

#' Find the path to Python associated with a conda environment
#'
#' @inheritParams reticulate::conda_python
conda_python <- function(envname = NULL, conda = "auto", all = FALSE) {
  envname <- get_envname(envname)
  envname <- reticulate:::condaenv_resolve(envname)
  if (grepl("[/\\\\]", envname)) {
    suffix <- if (reticulate:::is_windows()) "python.exe" else "bin/python"
    path <- file.path(envname, suffix)
    if (file.exists(path)) {
      return(path)
    }
    fmt <- "no conda environment exists at path '%s'"
    stop(sprintf(fmt, envname))
  }
  conda_envs <- reticulate::conda_list(conda = conda)
  conda_envs <- conda_envs[grep(normalizePath(reticulate:::conda_info(conda = conda)$envs_dirs[1], mustWork = FALSE), x = normalizePath(conda_envs$python, mustWork = FALSE), fixed = TRUE), , drop = FALSE]
  env <- conda_envs[conda_envs$name == envname, , drop = FALSE]
  if (nrow(env) == 0) {
    stop("conda environment \"", envname, "\" not found")
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python), mustWork = FALSE))
}

#' Run Python code safely
#'
#' Execute Python code through reticulate with proper error handling
#'
#' @param command A string containing Python code to execute
#' @param envir The R environment where Python objects should be made available
#' @param stop_on_error Whether to stop execution when an error occurs. If FALSE, will return the error object instead.
#' @param use_scp_env Whether to ensure the SCP environment is used before running the command. If TRUE (default),
#'   this will attempt to configure reticulate to use the SCP conda environment before executing the command.
#'   When FALSE, it will use whatever Python environment reticulate is currently configured to use.
#'
#' @return Invisibly returns NULL on success, or an error object if stop_on_error is FALSE
#' @examples
#' \dontrun{
#' # Basic usage
#' run_Python("import numpy as np; print(np.array([1,2,3]))")
#' 
#' # Return error instead of stopping
#' error <- run_Python("import non_existent_module", stop_on_error = FALSE)
#' 
#' # Use a custom environment
#' my_env <- new.env()
#' run_Python("x = [1, 2, 3]", envir = my_env)
#' my_env$x  # Access the Python object from R
#' 
#' # Use a specific Python environment (not SCP)
#' run_Python("import sys; print(sys.executable)", use_scp_env = FALSE)
#' }
#' @export
run_Python <- function(command, envir = .GlobalEnv, stop_on_error = TRUE, use_scp_env = TRUE) {
  # Ensure reticulate is using the SCP environment if requested
  if (use_scp_env) {
    # Get the SCP environment name
    envname <- get_envname()
    
    # Check if reticulate is loaded
    if (!isNamespaceLoaded("reticulate")) {
      # Find conda and the Python path
      conda <- find_conda()
      if (!is.null(conda)) {
        envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
        env <- env_exist(conda = conda, envname = envname, envs_dir = envs_dir)
        
        if (isTRUE(env)) {
          python_path <- conda_python(conda = conda, envname = envname)
          
          # Set the Python path before loading reticulate
          Sys.setenv(RETICULATE_PYTHON = python_path)
          Sys.setenv(RETICULATE_PYTHON_ENV = envname)
        }
      }
    } else {
      # If reticulate is already loaded, try to use the SCP environment
      conda <- find_conda()
      if (!is.null(conda)) {
        python_path <- conda_python(conda = conda, envname = envname)
        reticulate::use_python(python_path, required = FALSE)
      }
    }
  }
  
  # Now execute the Python command
  result <- tryCatch(expr = {
    # Make sure reticulate is loaded
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("The reticulate package is required but could not be loaded")
    }
    
    eval(
      {
        reticulate::py_run_string(command)
      },
      envir = envir
    )
    NULL
  }, error = function(error) {
    if (stop_on_error) {
      message(error)
      stop("Failed to run \"", command, "\". Please check manually.")
    } else {
      error
    }
  })
  invisible(result)
}

#' Test Python version compatibility
#'
#' This function tests compatibility of different Python versions with SCP
#' by attempting to create environments for each supported version.
#'
#' @param versions Character vector of environment versions to test (e.g., "3.9-1", "3.10-1", etc.)
#' @param cleanup Logical, whether to remove the test environments after testing
#' @param test_prefix Prefix to add to test environment names (default: "test_")
#' @param verbose Logical, whether to print detailed output during tests
#'
#' @return A data frame with compatibility test results for each Python version
#' @export
#'
#' @examples
#' \dontrun{
#' # Test all supported Python versions
#' test_results <- TestPythonCompatibility()
#' print(test_results)
#' 
#' # Test only Python 3.10 and 3.11
#' test_results <- TestPythonCompatibility(versions = c("3.10-1", "3.11-1"))
#' }
TestPythonCompatibility <- function(versions = c("3.9-1", "3.10-1", "3.11-1", "3.12-1"), 
                                   cleanup = TRUE, 
                                   test_prefix = "test_", 
                                   verbose = TRUE) {
  results <- data.frame(
    version = character(),
    python_version = character(),
    success = logical(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  
  for (version in versions) {
    if (verbose) message("Testing Python compatibility for version: ", version)
    
    # Get Python version from requirements
    requirements <- try(Env_requirements(version = version), silent = TRUE)
    if (inherits(requirements, "try-error")) {
      results <- rbind(results, data.frame(
        version = version,
        python_version = NA,
        success = FALSE,
        error_message = "Invalid version specification",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    python_version <- requirements[["python"]]
    test_env_name <- paste0(test_prefix, "py", gsub("\\.", "", python_version))
    
    if (verbose) message("Creating test environment: ", test_env_name, " with Python ", python_version)
    
    # Try to create the environment
    env_result <- tryCatch({
      PrepareEnv(
        envname = test_env_name,
        version = version,
        force = TRUE
      )
      TRUE
    }, error = function(e) {
      list(success = FALSE, message = as.character(e))
    })
    
    if (is.logical(env_result) && env_result) {
      # Test importing key packages
      test_imports <- tryCatch({
        run_Python("import scanpy", stop_on_error = FALSE)
        run_Python("import matplotlib", stop_on_error = FALSE)
        run_Python("import numpy", stop_on_error = FALSE)
        NULL
      }, error = function(e) {
        as.character(e)
      })
      
      import_success <- is.null(test_imports)
      
      results <- rbind(results, data.frame(
        version = version,
        python_version = python_version,
        success = import_success,
        error_message = if (import_success) "All packages imported successfully" else test_imports,
        stringsAsFactors = FALSE
      ))
    } else {
      results <- rbind(results, data.frame(
        version = version,
        python_version = python_version,
        success = FALSE,
        error_message = if (is.list(env_result)) env_result$message else "Unknown error",
        stringsAsFactors = FALSE
      ))
    }
    
    # Clean up test environment if requested
    if (cleanup) {
      if (verbose) message("Cleaning up test environment: ", test_env_name)
      conda <- find_conda()
      if (!is.null(conda)) {
        tryCatch({
          envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
          env_path <- paste0(envs_dir, "/", test_env_name)
          if (file.exists(env_path)) {
            unlink(env_path, recursive = TRUE)
          }
        }, error = function(e) {
          if (verbose) message("Warning: Could not clean up environment: ", as.character(e))
        })
      }
    }
  }
  
  return(results)
}

#' Check and install python packages
#'
#' @param packages A character vector, indicating package names which should be installed or removed. Use \code{⁠<package>==<version>}⁠ to request the installation of a specific version of a package.
#' @param envname The name of a conda environment.
#' @param conda The path to a conda executable. Use \code{"auto"} to allow SCP to automatically find an appropriate conda binary.
#' @param force Whether to force package installation. Default is \code{FALSE}.
#' @param pip Whether to use pip for package installation. If TRUE (default), packages are installed with pip.
#'          If FALSE, packages are installed with conda, which may be faster but have less package compatibility.
#' @param pip_options An optional character vector of additional command line arguments to be passed to \code{pip}. Only relevant when \code{pip = TRUE}.
#' @param use_parallel Whether to use parallel installation for packages. May speed up installation but can cause issues on some systems. Default is FALSE.
#' @param verbose Whether to display detailed messages during installation. Default is TRUE.
#' @param ... Other arguments passed to \code{\link[reticulate]{conda_install}}
#'
#' @examples
#' check_Python(packages = c("bbknn", "scanorama"))
#' \dontrun{
#' # Install a specific version using pip
#' check_Python(packages = "scvi-tools==0.20.0", pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple")
#' 
#' # Use conda for faster installation (but possible compatibility issues)
#' check_Python(packages = c("scanpy", "leidenalg"), pip = FALSE)
#' }
#' @export
check_Python <- function(packages, envname = NULL, conda = "auto", force = FALSE, 
                        pip = TRUE, pip_options = character(), use_parallel = FALSE,
                        verbose = TRUE, ...) {
  start_time <- Sys.time()
  
  # Resolve environment name
  envname <- get_envname(envname)
  
  # Find conda
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  
  # Check if environment exists
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    warning(envname, " python environment does not exist. Create it with the PrepareEnv function...", immediate. = TRUE)
    PrepareEnv()
  }
  
  # Determine which packages need installation
  if (isTRUE(force)) {
    # Force reinstallation of all packages
    pkg_installed <- setNames(rep(FALSE, length(packages)), packages)
    if (pip) {
      pip_options <- c(pip_options, "--force-reinstall")
    }
  } else {
    # Check which packages are already installed
    if (verbose) message("Checking installed packages...")
    pkg_installed <- exist_Python_pkgs(packages = packages, envname = envname, conda = conda)
  }
  
  # Install packages if needed
  missing_pkgs <- sum(!pkg_installed)
  if (missing_pkgs > 0) {
    pkgs_to_install <- names(pkg_installed)[!pkg_installed]
    
    if (verbose) {
      if (missing_pkgs == 1) {
        message("Installing 1 package: ", pkgs_to_install)
      } else {
        message("Installing ", missing_pkgs, " packages: ", paste0(pkgs_to_install, collapse = ", "))
      }
      message("Installation method: ", if (pip) "pip" else "conda")
    }
    
    # Add pip itself if using pip
    if (isTRUE(pip)) {
      pkgs_to_install <- c("pip", pkgs_to_install)
    }
    
    # Use parallel installation if requested
    if (use_parallel && missing_pkgs > 1) {
      # Split packages into batches for parallel installation
      if (verbose) message("Using parallel installation...")
      
      # For conda installation, we don't split as conda handles dependencies better in one batch
      if (!pip) {
        installation_result <- tryCatch(expr = {
          conda_install(
            conda = conda, 
            packages = pkgs_to_install, 
            envname = envname, 
            pip = pip, 
            pip_options = pip_options,
            ...
          )
          TRUE
        }, error = function(e) {
          if (verbose) message("Installation error: ", e$message)
          FALSE
        })
      } else {
        # For pip, we can install in parallel batches
        batch_size <- min(5, ceiling(length(pkgs_to_install) / 2))
        pkg_batches <- split(pkgs_to_install, ceiling(seq_along(pkgs_to_install) / batch_size))
        
        installation_results <- sapply(pkg_batches, function(batch) {
          if (verbose) message("Installing batch: ", paste(batch, collapse = ", "))
          tryCatch(expr = {
            conda_install(
              conda = conda, 
              packages = batch, 
              envname = envname, 
              pip = pip, 
              pip_options = pip_options,
              ...
            )
            TRUE
          }, error = function(e) {
            if (verbose) message("Batch installation error: ", e$message)
            FALSE
          })
        })
        
        installation_result <- all(installation_results)
      }
    } else {
      # Sequential installation
      installation_result <- tryCatch(expr = {
        conda_install(
          conda = conda, 
          packages = pkgs_to_install, 
          envname = envname, 
          pip = pip, 
          pip_options = pip_options,
          ...
        )
        TRUE
      }, error = function(e) {
        if (verbose) message("Installation error: ", e$message)
        FALSE
      })
    }
    
    # If conda/bulk installation failed, try installing packages one by one
    if (!installation_result && verbose) {
      message("Bulk installation failed, trying individual package installation...")
      for (pkg in unique(pkgs_to_install)) {
        if (pkg == "pip") next  # Skip pip itself in individual installation
        
        if (verbose) message("Installing ", pkg, "...")
        tryCatch(expr = {
          conda_install(
            conda = conda, 
            packages = pkg, 
            envname = envname, 
            pip = pip, 
            pip_options = pip_options,
            ...
          )
        }, error = function(e) {
          if (verbose) message("Failed to install ", pkg, ": ", e$message)
        })
      }
    }
  } else if (verbose) {
    message("All required packages are already installed.")
  }
  
  # Final check that packages are installed
  pkg_installed <- exist_Python_pkgs(packages = packages, envname = envname, conda = conda)
  missing_pkgs <- sum(!pkg_installed)
  
  if (missing_pkgs > 0) {
    failed_pkgs <- names(pkg_installed)[!pkg_installed]
    stop("Failed to install ", missing_pkgs, " package(s): ", paste0(failed_pkgs, collapse = ", "), 
         " into the environment \"", envname, "\".\n",
         "Try installing manually or consider using ", if(pip) "conda (pip=FALSE)" else "pip (pip=TRUE)", 
         " instead.")
  } else {
    end_time <- Sys.time()
    if (verbose) {
      message("✓ All packages installed successfully in ", round(difftime(end_time, start_time, units = "secs"), 1), " seconds")
    }
    return(invisible(TRUE))
  }
}

#' Check and install R packages
#'
#' @param packages Package to be installed. Package source can be CRAN, Bioconductor or Github, e.g. scmap, quadbiolab/simspec.
#' By default, the package name is extracted according to the \code{packages} parameter.
#' @param install_methods Functions used to install R packages.
#' @param lib  The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages. Default is \code{FALSE}.
#'
#' @importFrom utils packageVersion
#' @export
check_R <- function(packages, install_methods = c("BiocManager::install", "install.packages", "devtools::install_github"), lib = .libPaths()[1], force = FALSE) {
  status_list <- list()
  for (pkg in packages) {
    version <- NULL
    if (grepl("/", pkg)) {
      # github package
      pkg_name <- strsplit(pkg, split = "/|@|==", perl = TRUE)[[1]][[2]]
    } else {
      pkg_name <- strsplit(pkg, split = "@|==", perl = TRUE)[[1]][[1]]
      if (length(strsplit(pkg, split = "@|==", perl = TRUE)[[1]]) > 1) {
        version <- strsplit(pkg, split = "@|==", perl = TRUE)[[1]][[2]]
      }
    }
    dest <- gsub("@.*|==.*|>=.*", "", pkg)
    if (is.null(version)) {
      force_update <- isTRUE(force)
    } else {
      force_update <- isTRUE(packageVersion(pkg_name) < package_version(version)) || isTRUE(force)
    }
    if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE)) || isTRUE(force_update)) {
      message("Install package: \"", pkg_name, "\" ...")
      status_list[[pkg]] <- FALSE
      i <- 1
      while (isFALSE(status_list[[pkg]])) {
        tryCatch(expr = {
          if (grepl("BiocManager", install_methods[i])) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", lib = lib)
            }
            eval(str2lang(paste0(install_methods[i], "(\"", dest, "\", lib=\"", lib, "\", update = FALSE, upgrade = \"never\", ask = FALSE, force = TRUE)")))
          } else if (grepl("devtools", install_methods[i])) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
              install.packages("devtools", lib = lib)
            }
            if (!requireNamespace("withr", quietly = TRUE)) {
              install.packages("withr", lib = lib)
            }
            eval(str2lang(paste0("withr::with_libpaths(new = \"", lib, "\", ", install_methods[i], "(\"", dest, "\", upgrade = \"never\", force = TRUE))")))
          } else {
            eval(str2lang(paste0(install_methods[i], "(\"", dest, "\", lib=\"", lib, "\", force = TRUE)")))
          }
        }, error = function(e) {
          status_list[[pkg]] <- FALSE
        })
        status_list[[pkg]] <- requireNamespace(pkg_name, quietly = TRUE)
        i <- i + 1
        if (i > length(install_methods)) {
          break
        }
      }
    } else {
      status_list[[pkg]] <- TRUE
    }
  }
  out <- sapply(status_list, isTRUE)
  out <- out[!out]
  if (length(out) > 0) {
    stop("Failed to install the package(s): ", paste0(names(out), collapse = ","), ". Please install manually.")
  }
}

#' Try to evaluate an expression a set number of times before failing
#'
#' The function is used as a fail-safe if your R code sometimes works and sometimes
#' doesn't, usually because it depends on a resource that may be temporarily
#' unavailable. It tries to evaluate the expression `max_tries` times. If all the
#' attempts fail, it throws an error; if not, the evaluated expression is returned.
#'
#' @param expr The expression to be evaluated.
#' @param max_tries The maximum number of attempts to evaluate the expression before giving up. Default is set to 5.
#' @param error_message a string, additional custom error message you would like to be displayed when an error occurs.
#' @param retry_message a string, a message displayed when a new try to evaluate the expression would be attempted.
#'
#' @return This function returns the evaluated expression if successful, otherwise it throws an error if all attempts are unsuccessful.
#' @export
#'
#' @examples
#' f <- function() {
#'   value <- runif(1, min = 0, max = 1)
#'   if (value > 0.5) {
#'     message("value is larger than 0.5")
#'     return(value)
#'   } else {
#'     stop("value is smaller than 0.5")
#'   }
#' }
#' f_evaluated <- try_get(expr = f())
#' print(f_evaluated)
#'
try_get <- function(expr, max_tries = 5, error_message = "", retry_message = "Retrying...") {
  out <- simpleError("start")
  ntry <- 0
  while (inherits(out, "error")) {
    ntry <- ntry + 1
    # print(paste0("ntry: ", ntry, collapse = ""))
    out <- tryCatch(
      expr = eval.parent(substitute(expr)),
      error = function(error) {
        message(error)
        message("")
        message(error_message)
        Sys.sleep(1)
        return(error)
      }
    )
    if (inherits(out, "error") && ntry >= max_tries) {
      stop(out)
    } else {
      if (!inherits(out, "error")) {
        break
      } else {
        message(retry_message)
      }
    }
  }
  return(out)
}

#' Install all dependencies required for SCP functionality
#'
#' This function installs all required and recommended R packages needed for SCP to function properly.
#' It is particularly useful when setting up a fresh installation on a new machine.
#'
#' @param bioc_deps Whether to install Bioconductor dependencies. Default is TRUE.
#' @param optional_deps Whether to install optional dependencies that enable additional functionality. Default is TRUE.
#' @param update Whether to update packages that are already installed. Default is FALSE.
#' @param lib The library location to install packages into. Defaults to .libPaths()[1].
#' 
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Install all dependencies
#' install_all_dependencies()
#' 
#' # Install only required dependencies
#' install_all_dependencies(bioc_deps = FALSE, optional_deps = FALSE)
#' }
install_all_dependencies <- function(bioc_deps = TRUE, optional_deps = TRUE, update = FALSE, lib = .libPaths()[1]) {
  # Core dependencies - essential for basic functionality
  core_deps <- c(
    "Matrix", 
    "Seurat", 
    "SeuratObject", 
    "reticulate", 
    "rlang",
    "dplyr",
    "ggplot2", 
    "future",
    "methods",
    "FNN",
    "irlba",
    "igraph"
  )
  
  # Bioconductor dependencies - required packages from Bioconductor
  bioc_required <- c(
    "AnnotationDbi",
    "biomaRt",
    "BiocParallel",
    "ComplexHeatmap",
    "GO.db",
    "GOSemSim",
    "HDF5Array",
    "rhdf5",
    "clusterProfiler",
    "slingshot"
  )
  
  # Bioconductor dependencies - suggested packages from Bioconductor
  bioc_suggested <- c(
    "AUCell",
    "batchelor",
    "Biobase",
    "BiocGenerics",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "PFAM.db",
    "reactome.db",
    "scds",
    "scDblFinder",
    "scmap",
    "SingleR",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors"
  )
  
  # Development tools
  dev_tools <- c(
    "devtools",
    "withr",
    "testthat"
  )
  
  # Optional packages for extended functionality
  optional_pkgs <- c(
    "patchwork",
    "plotly",
    "ggnewscale",
    "ggrepel",
    "ggforce",
    "pbapply",
    "data.table",
    "intrinsicDimension",
    "uwot",
    "simplifyEnrichment",
    "circlize",
    "Signac",
    "R.cache",
    "R.utils",
    "proxyC",
    "reshape2",
    "scales",
    "mgcv",
    "png",
    "future.apply"
  )
  
  install_status <- list(
    cran = list(),
    bioc = list()
  )
  
  message("Installing core dependencies...")
  for (pkg in core_deps) {
    if (!requireNamespace(pkg, quietly = TRUE) || update) {
      tryCatch({
        utils::install.packages(pkg, lib = lib)
        install_status$cran[[pkg]] <- TRUE
      }, error = function(e) {
        message("Failed to install ", pkg, ": ", e$message)
        install_status$cran[[pkg]] <- FALSE
      })
    } else {
      install_status$cran[[pkg]] <- TRUE
    }
  }
  
  # Ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE) || update) {
    tryCatch({
      utils::install.packages("BiocManager", lib = lib)
      install_status$cran[["BiocManager"]] <- TRUE
    }, error = function(e) {
      message("Failed to install BiocManager: ", e$message)
      install_status$cran[["BiocManager"]] <- FALSE
    })
  } else {
    install_status$cran[["BiocManager"]] <- TRUE
  }
  
  # Install Bioconductor packages if BiocManager is available
  if (bioc_deps && requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing required Bioconductor dependencies...")
    for (pkg in bioc_required) {
      if (!requireNamespace(pkg, quietly = TRUE) || update) {
        tryCatch({
          BiocManager::install(pkg, update = FALSE, ask = FALSE, lib = lib)
          install_status$bioc[[pkg]] <- TRUE
        }, error = function(e) {
          message("Failed to install Bioconductor package ", pkg, ": ", e$message)
          install_status$bioc[[pkg]] <- FALSE
        })
      } else {
        install_status$bioc[[pkg]] <- TRUE
      }
    }
    
    if (optional_deps) {
      message("Installing optional Bioconductor dependencies...")
      for (pkg in bioc_suggested) {
        if (!requireNamespace(pkg, quietly = TRUE) || update) {
          tryCatch({
            BiocManager::install(pkg, update = FALSE, ask = FALSE, lib = lib)
            install_status$bioc[[pkg]] <- TRUE
          }, error = function(e) {
            message("Failed to install optional Bioconductor package ", pkg, ": ", e$message)
            install_status$bioc[[pkg]] <- FALSE
          })
        } else {
          install_status$bioc[[pkg]] <- TRUE
        }
      }
    }
  } else if (bioc_deps) {
    message("IMPORTANT: BiocManager could not be installed. Bioconductor packages cannot be installed.")
  }
  
  message("Installing development tools...")
  for (pkg in dev_tools) {
    if (!requireNamespace(pkg, quietly = TRUE) || update) {
      tryCatch({
        utils::install.packages(pkg, lib = lib)
        install_status$cran[[pkg]] <- TRUE
      }, error = function(e) {
        message("Failed to install ", pkg, ": ", e$message)
        install_status$cran[[pkg]] <- FALSE
      })
    } else {
      install_status$cran[[pkg]] <- TRUE
    }
  }
  
  if (optional_deps) {
    message("Installing optional CRAN packages...")
    for (pkg in optional_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE) || update) {
        tryCatch({
          utils::install.packages(pkg, lib = lib)
          install_status$cran[[pkg]] <- TRUE
        }, error = function(e) {
          message("Failed to install ", pkg, ": ", e$message)
          install_status$cran[[pkg]] <- FALSE
        })
      } else {
        install_status$cran[[pkg]] <- TRUE
      }
    }
  }
  
  # Summarize installation results
  failed_cran <- names(install_status$cran)[!unlist(install_status$cran)]
  failed_bioc <- names(install_status$bioc)[!unlist(install_status$bioc)]
  
  if (length(failed_cran) > 0 || length(failed_bioc) > 0) {
    message("\n==== INSTALLATION SUMMARY ====")
    message("Some packages could not be installed automatically.")
    
    if (length(failed_cran) > 0) {
      message("\nFailed CRAN packages:")
      message(paste(" -", failed_cran, collapse = "\n"))
    }
    
    if (length(failed_bioc) > 0) {
      message("\nFailed Bioconductor packages:")
      message(paste(" -", failed_bioc, collapse = "\n"))
    }
    
    message("\nTo install failed packages manually:")
    
    if (length(failed_cran) > 0) {
      cran_cmd <- sprintf('install.packages(c("%s"))', paste(failed_cran, collapse = '", "'))
      message("\n# Install CRAN packages")
      message(cran_cmd)
    }
    
    if (length(failed_bioc) > 0) {
      bioc_cmd <- sprintf('BiocManager::install(c("%s"))', paste(failed_bioc, collapse = '", "'))
      message("\n# Install Bioconductor packages")
      message(bioc_cmd)
    }
  } else {
    message("All dependencies were installed successfully!")
  }
  
  invisible(NULL)
}

#' Download File from the Internet
#'
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files. The default is to try different download methods in turn until the download is successfully completed.
#' @param max_tries Number of tries for each download method.
#' @param ... Other arguments passed to \code{\link[utils]{download.file}}
#'
#' @importFrom utils download.file
#' @export
download <- function(url, destfile, methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"), quiet = FALSE, ..., max_tries = 2) {
  if (missing(url) || missing(destfile)) {
    stop("'url' and 'destfile' must be both provided.")
  }
  ntry <- 0
  status <- NULL
  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(expr = {
        suppressWarnings(download.file(url, destfile = destfile, method = method, quiet = quiet, ...))
        status <- 1
      }, error = function(error) {
        message(error)
        message("Cannot download from the url: ", url)
        message("Failed to download using \"", method, "\". Retry...\n")
        Sys.sleep(1)
        return(NULL)
      })
      if (!is.null(status)) {
        break
      }
    }
    ntry <- ntry + 1
    if (is.null(status) && ntry >= max_tries) {
      stop("Download failed.")
    }
  }
  return(invisible(NULL))
}

kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(url = url, destfile = temp)
  content <- as.data.frame(do.call(rbind, strsplit(readLines(temp), split = "\t")))
  return(content)
}

rescale <- function(x, from = range(x, na.rm = TRUE, finite = TRUE), to = c(0, 1)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  } else {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  }
}

zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) {
    stop("x must be length 1 or 2")
  }
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2]) / m) < tol
}

#' @importFrom grDevices col2rgb rgb
col2hex <- function(cname) {
  colMat <- col2rgb(cname)
  rgb(red = colMat[1, ] / 255, green = colMat[2, ] / 255, blue = colMat[3, ] / 255)
}

#' Invoke a function with a list of arguments
#' @param .fn A function, or function name as a string.
#' @param .args A list of arguments.
#' @param ... Other arguments passed to the function.
#' @param .env Environment in which to evaluate the call. This will be most useful if .fn is a string, or the function has side-effects.
#' @importFrom rlang caller_env is_null is_scalar_character is_character is_function set_names env env_get env_bind syms call2
#' @export
invoke <- function(.fn, .args = list(), ..., .env = caller_env()) {
  args <- c(.args, list(...))
  .bury <- c(".fn", "")
  if (is_null(.bury) || !length(args)) {
    if (is_scalar_character(.fn)) {
      .fn <- env_get(.env, .fn, inherit = TRUE)
    }
    call <- call2(.fn, !!!args)
    return(.External2(rlang:::ffi_eval, call, .env))
  }
  if (!is_character(.bury, 2L)) {
    abort("`.bury` must be a character vector of length 2")
  }
  arg_prefix <- .bury[[2]]
  fn_nm <- .bury[[1]]
  buried_nms <- paste0(arg_prefix, seq_along(args))
  buried_args <- set_names(args, buried_nms)
  .env <- env(.env, !!!buried_args)
  args <- set_names(buried_nms, names(args))
  args <- syms(args)
  if (is_function(.fn)) {
    env_bind(.env, `:=`(!!fn_nm, .fn))
    .fn <- fn_nm
  }
  call <- call2(.fn, !!!args)
  .External2(rlang:::ffi_eval, call, .env)
}

#' Implement similar functions to the \code{unnest} function in the tidyr package
#' @param data A data frame.
#' @param cols Columns to unnest.
#' @param keep_empty By default, you get one row of output for each element of the list your unchopping/unnesting. This means that if there's a size-0 element (like \code{NULL} or an empty data frame), that entire row will be dropped from the output. If you want to preserve all rows, use \code{keep_empty = TRUE} to replace size-0 elements with a single row of missing values.
#' @export
#'
#' @importFrom Seurat GetAssayData LayerData DefaultAssay
#' 
#' #' Get data from Seurat object in a version-agnostic way
#' #'
#' #' This function retrieves data from a Seurat object, handling both V4 and V5 versions
#' #' appropriately by using either GetAssayData or LayerData.
#' #'
#' #' @param srt A Seurat object
#' #' @param slot The data slot/layer to access ("counts", "data", or "scale.data")
#' #' @param assay Name of assay to use (defaults to DefaultAssay)
#' #'
#' #' @return A matrix or sparse matrix of the requested data
#' #' @export
get_seurat_data <- function(srt, slot = "data", assay = NULL) {
  is_v5 <- is_seurat_v5(srt)
  assay <- assay %||% DefaultAssay(srt)
  
  if (is_v5) {
    return(LayerData(srt[[assay]], layer = slot))
  } else {
    return(GetAssayData(srt, slot = slot, assay = assay))
  }
}
unnest <- function(data, cols, keep_empty = FALSE) {
  if (nrow(data) == 0 || length(cols) == 0) {
    return(data)
  }
  for (col in cols) {
    col_expand <- unlist(data[[col]])
    expand_times <- sapply(data[[col]], length)
    if (isTRUE(keep_empty)) {
      data[[col]][expand_times == 0] <- NA
      col_expand <- unlist(data[[col]])
      expand_times[expand_times == 0] <- 1
    }
    data <- data[rep(seq_len(nrow(data)), times = expand_times), ]
    data[, col] <- col_expand
  }
  rownames(data) <- NULL
  return(data)
}

#' Attempts to turn a dgCMatrix into a dense matrix
#'
#' @examples
#' data("pancreas_sub")
#' system.time(mat1 <- as.matrix(slot(pancreas_sub[["RNA"]], "counts")))
#' system.time(mat2 <- as_matrix(slot(pancreas_sub[["RNA"]], "counts")))
#' identical(mat1, mat2)
#'
#' @param x A matrix.
#' @useDynLib SCP
#' @importFrom Matrix as.matrix
#' @export
as_matrix <- function(x) {
  if (!inherits(matrix, "dgCMatrix")) {
    return(as.matrix(x))
  } else {
    row_pos <- x@i
    col_pos <- findInterval(seq_along(x@x) - 1, x@p[-1])
    out <- asMatrix(rp = row_pos, cp = col_pos, z = x@x, nrows = x@Dim[1], ncols = x@Dim[2])
    attr(out, "dimnames") <- list(x@Dimnames[[1]], x@Dimnames[[2]])
    return(out)
  }
}

#' Capitalizes the characters
#' Making the first letter uppercase
#'
#' @examples
#' x <- c("dna methylation", "rRNA processing", "post-Transcriptional gene silencing")
#' capitalize(x)
#' @param x A vector of character strings to be capitalized.
#' @param force_tolower Whether to force the remaining letters to be lowercase.
#' @export
capitalize <- function(x, force_tolower = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  if (!inherits(x, "character")) {
    stop("x must be the type of character.")
  }
  if (isTRUE(force_tolower)) {
    x <- paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
  } else {
    first_word <- sapply(strsplit(x, "\\s|-"), function(s) s[1])
    index <- which(first_word == tolower(first_word))
    x[index] <- paste(toupper(substr(x[index], 1, 1)), substr(x[index], 2, nchar(x[index])), sep = "")
  }
  return(x)
}

str_wrap <- function(x, width = 80) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(lapply(x, function(i) paste0(strwrap(i, width = width), collapse = "\n")))
  return(x_wrap)
}

#' Split a vector into the chunks
#'
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' x <- 1:10
#' names(x) <- letters[1:10]
#' tochunks(x, nchunks = 3)
#' @export
tochunks <- function(x, nchunks) {
  split(x, cut(seq_along(x), nchunks, labels = FALSE))
}

#' Generate a iterator along chunks of a vector
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' \dontrun{
#' library(BiocParallel)
#' x <- 1:100
#' BPPARAM <- bpparam()
#' bpprogressbar(BPPARAM) <- TRUE
#' bpworkers(BPPARAM) <- 10
#' slow_fun <- function(x) {
#'   out <- NULL
#'   for (i in seq_along(x)) {
#'     Sys.sleep(0.5)
#'     out[[i]] <- x[[i]] + 3
#'   }
#'   return(out)
#' }
#' system.time({
#'   res0 <- lapply(x, FUN = slow_fun)
#' })
#' unlist(res0, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res1 <- bplapply(x, FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res1, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res2 <- bplapply(tochunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res2, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res3 <- bpiterate(ITER = iterchunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res3, recursive = FALSE, use.names = FALSE)[71:73]
#' }
#' @export
iterchunks <- function(x, nchunks) {
  chunks <- tochunks(x, nchunks)
  i <- 0L
  function() {
    if (i >= length(chunks)) {
      return(NULL)
    }
    i <<- i + 1L
    x[chunks[[i]]]
  }
}
