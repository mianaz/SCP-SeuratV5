#' Validation utility functions for SCP package
#'
#' These internal helper functions provide consistent validation across the package,
#' reducing code duplication and ensuring uniform error messages.
#'
#' @name validation-utils
#' @keywords internal
NULL

#' Validate and standardize group parameter
#'
#' This function handles the common pattern of validating group parameters
#' that can be either a vector of values or a column name in metadata.
#'
#' @param srt Seurat object
#' @param group Group parameter value (vector or column name)
#' @param param_name Name of the parameter for error messages
#' @param target_slot Target slot name to assign the group to (default: param_name)
#' @return Modified Seurat object with group information added to metadata
#' @keywords internal
validate_group_parameter <- function(srt, group, param_name = "group",
                                    target_slot = NULL) {

  if (is.null(target_slot)) {
    target_slot <- param_name
  }

  if (!is.null(group)) {
    n_cells <- ncol(srt)

    if (length(group) == n_cells) {
      # Group is a vector with one value per cell
      srt[[target_slot]] <- group
    } else if (length(group) == 1) {
      # Group is a column name in metadata
      if (!group %in% colnames(srt@meta.data)) {
        stop(sprintf("'%s' must be one of the column names in the meta.data: %s",
                    param_name, paste(colnames(srt@meta.data), collapse = ", ")),
             call. = FALSE)
      }
      srt[[target_slot]] <- srt[[group, drop = TRUE]]
    } else {
      stop(sprintf("Length of '%s' must be 1 (column name) or %d (number of cells), not %d",
                  param_name, n_cells, length(group)),
           call. = FALSE)
    }
  }

  return(srt)
}


#' Require one or more packages
#'
#' Checks that required packages are installed, provides helpful error message.
#' This consolidates the 94 instances of package checking across the codebase.
#'
#' @param packages Character vector of package names
#' @param context Optional context message explaining why packages are needed
#' @param quietly Whether to check quietly
#' @return Invisible TRUE if all packages available
#' @keywords internal
require_packages <- function(packages, context = NULL, quietly = TRUE) {
  missing_packages <- character(0)

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = quietly)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }

  if (length(missing_packages) > 0) {
    context_msg <- if (!is.null(context)) {
      sprintf(" %s requires", context)
    } else {
      ""
    }

    install_msg <- if (length(missing_packages) == 1) {
      sprintf("\n\nInstall with: install.packages('%s')", missing_packages[1])
    } else {
      sprintf("\n\nInstall with: install.packages(c(%s))",
             paste0("'", missing_packages, "'", collapse = ", "))
    }

    stop(sprintf("Missing required package%s:%s %s%s",
                if (length(missing_packages) > 1) "s" else "",
                context_msg,
                paste(missing_packages, collapse = ", "),
                install_msg),
         call. = FALSE)
  }

  invisible(TRUE)
}

