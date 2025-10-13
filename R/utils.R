#' @importFrom BiocParallel bpparam
#' @importFrom Matrix t
#' @importFrom reticulate py_module_available conda_list use_python import_from_path py_available
#' @importFrom stats as.dist as.formula median model.matrix prcomp sd var aggregate kmeans na.omit setNames quantile weighted.mean
#' @importFrom utils askYesNo head installed.packages menu modifyList packageVersion read.csv read.table setTxtProgressBar tail txtProgressBar download.file
#' @importFrom grDevices col2rgb colorRampPalette dev.cur dev.interactive dev.new dev.off devAskNewPage palette pdf recordPlot rgb
#' @importFrom SeuratObject DefaultAssay CreateAssayObject CreateDimReducObject
#' @importFrom SeuratObject LayerData
#' @importFrom Seurat as.CellDataSet
#' @importFrom tidyr unnest
#' @importFrom stringr str_wrap
#' @import ggplot2
#' @import Seurat
#' @import future
#' @import future.apply
#' @import methods

# Declare global variables used in NSE contexts (ggplot, dplyr, tidyr, etc.)
utils::globalVariables(c(
  # Base plotting and aesthetics
  "x", "y", "z", "value", "variable", "Var1", "Var2",
  "xend", "yend", "hjust", "vjust", "angle", "color", "fill",
  "size", "shape", "alpha", "linetype", "linewidth",
  # Grouping and labeling
  ".", "group", "label", "name", "count", "frequency", "percent",
  # Common single-cell analysis variables
  "gene", "error", "from_dim1", "to_dim1", "from_dim2", "to_dim2",
  "feature", "cell", "cluster", "celltype", "sample", "condition",
  # Data manipulation variables
  "row", "col", "id", "type", "category", "class", "level",
  # Statistics and metrics
  "mean", "median", "sd", "var", "min", "max", "sum",
  "pct", "ratio", "score", "weight", "rank",
  # Enrichment and pathway analysis
  "pathway", "term", "geneID", "description", "pvalue", "qvalue",
  "padj", "p.adjust", "GeneRatio", "BgRatio", "richFactor",
  # Dimensionality reduction
  "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2", "PC_1", "PC_2",
  "dim1", "dim2", "component", "embedding",
  # Trajectory and pseudotime
  "pseudotime", "trajectory", "branch", "state", "fate",
  # Gene expression
  "expression", "logFC", "avg_log2FC", "pct.1", "pct.2",
  "avg.exp", "avg.exp.scaled"
))

#' Check if Python packages are available
#'
#' This function checks if required Python packages are available in the current Python environment.
#'
#' @param packages Character vector of Python package names to check
#' @return Invisible NULL. Throws an error if packages are not available.
#' @export
check_Python <- function(packages) {
  if (!reticulate::py_available()) {
    stop("Python is not available. Please set up Python environment with PrepareEnv()")
  }

  for (pkg in packages) {
    if (!reticulate::py_module_available(pkg)) {
      stop("Python package '", pkg, "' is required but not installed.\n",
           "Please install it in your Python environment or run PrepareEnv()")
    }
  }

  return(invisible(NULL))
}

#' @export
palette_default <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#8DA0CB", "#FC8D62", "#66C2A5", "#E6F5C9", "#FFF2AE", "#F4CAE4", "#F1E2CC"
)

#' @export
palette_category <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999",
  "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#E5C494", "#E78AC3", "#A6D854", "#FFD92F", "#FC8D62", "#66C2A5", "#BC80BD", "#CCEBC5"
)

#' @export
palette_dimplot <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
  "#E78AC3", "#A6D854", "#FFD92F", "#FC8D62", "#66C2A5", "#8DA0CB", "#E5C494",
  "#B3B3B3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
  "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#F4CAE4", "#F1E2CC", "#CCCCCC"
)

# Color palette for visualization
# Credit to Erica for the color palette (used in VoxHunt)
zissou <- c(
  "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00", "#A50026", "#860308", "#631109", "#4E190B",
  "#85B8D2", "#97C3CF", "#A8CDCD", "#BAD8CA", "#CBE2C8", "#DDEDC5", "#EEF7C3", "#FFFFC0", "#F8EEB3",
  "#F1DDA6", "#EBCC99", "#E4BB8D", "#DDAA80", "#D79973", "#D08866", "#CA775A", "#9D0014", "#B20D11",
  "#C61A0F", "#DB280C", "#EF350A", "#FF5207", "#FF6E05", "#FF8A02", "#FFA600", "#E29C11", "#C59321",
  "#A88932", "#8B8042", "#6E7653", "#526D63", "#356374", "#185A84", "#215E88", "#2A628C", "#336690",
  "#3B6B93", "#446F97", "#4D739B", "#56779F", "#5E7CA3", "#6780A7", "#7084AB", "#7888AF", "#818DB2",
  "#8A91B6", "#9395BA", "#9B99BE", "#A49EC2", "#ADA2C6", "#B6A6CA", "#BEAACE", "#C7AFD1", "#D0B3D5",
  "#D9B7D9", "#E1BBDD", "#E7BFE0", "#EEC4E3", "#F4C8E7", "#FBCCEA", "#FFD0ED", "#FFD5EF", "#FFD9F0",
  "#FFDDF2", "#FFE1F3", "#FFE6F5", "#FFEAF6", "#FFEEF8", "#FFF2F9", "#FFF7FB", "#FFFBFC", "#FFFFFE"
)


scp_colorname <- function() {
  cnames <- c(
    "default", "default2", "category", "dimplot", "lace",
    "npr", "nejm", "jama", "lancet", "jco", "ucscgb", "d3", "igv", "uchicago", "simpsons", "rickandmorty",
    "alphabet", "alphabet2", "paired", "set1", "set2", "set3", "dark2", "accent", "pastel1", "pastel2",
    grDevices::colors()[!grepl(pattern = "\\d", x = grDevices::colors())],
    grDevices::hcl.pals()
  )
  return(cnames)
}

CheckMatrix <- function(matrix, cells = NULL, features = NULL, n = NULL) {
  if (!inherits(matrix, c("matrix", "Matrix", "data.frame"))) {
    stop("'matrix' is not a matrix or data.frame.")
  }
  if (!is.null(cells)) {
    if (!any(cells %in% colnames(matrix))) {
      stop("None of the cells provided is in the colnames of the matrix.")
    }
    matrix <- matrix[, cells[cells %in% colnames(matrix)], drop = FALSE]
  }
  if (!is.null(features)) {
    if (!any(features %in% rownames(matrix))) {
      stop("None of the features provided is in the rownames of the matrix.")
    }
    matrix <- matrix[features[features %in% rownames(matrix)], , drop = FALSE]
  }
  if (!is.null(n)) {
    if (n > ncol(matrix)) {
      stop("'n' is larger than the number of columns in 'matrix'.")
    }
    matrix <- matrix[, sample(colnames(matrix), n), drop = FALSE]
  }
  return(matrix)
}

set_scp_palette <- function(x, palette = "default", palcolor = NULL, matched = FALSE, reverse = FALSE) {
  if (length(unique(x)) > 100) {
    stop("x has too many unique values. Unable to set colors.")
  }
  ggsci_db <- utils::getFromNamespace("ggsci_db", "ggsci")
  if (length(palette) > 1) {
    pal <- palette
  } else {
    if (is.null(palette)) {
      pal <- pal.def <- palette_default
    } else if (palette == "default") {
      pal <- pal.def <- palette_default
    } else if (palette == "default2") {
      pal <- pal.def <- c(
        "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100",
        "#6699CC", "#888888", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#CC79A7", "#999999"
      )
    } else if (palette == "category") {
      pal <- pal.def <- palette_category
    } else if (palette == "dimplot") {
      pal <- pal.def <- palette_dimplot
    } else if (palette == "lace") {
      pal <- pal.def <- c(
        "#E59060", "#2F4F4F", "#E4B39C", "#83C5A6", "#548687", "#E4E1A7", "#95A870", "#D19FA1", "#75AADB", "#FFFFB3", "#FFFACD",
        "#C12E34", "#D16B6A", "#E0301E", "#EB8C87", "#DC5034", "#DC9B7F", "#FFE3AA", "#FED5AD", "#32CD32", "#67C873", "#4CB877",
        "#228B22", "#8FBC8F", "#FFC125", "#B8860B", "#BC763C", "#9932CC", "#8B008B", "#E6E6FA", "#C89EC4", "#800080", "#E0B0FF"
      )
    } else if (palette == "hue") {
      pal <- grDevices::hcl(h = seq(15, 375, length = length(unique(x)) + 1), l = 65, c = 100)[seq_len(length(unique(x)))]
      pal.def <- pal
    } else if (palette == "grey") {
      pal <- grDevices::grey.colors(n = length(unique(x)), start = 0.3, end = 0.9)
      pal.def <- pal
    } else if (palette == "heat") {
      pal.raw <- grDevices::heat.colors(100, rev = TRUE)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else if (palette == "terrain") {
      pal.raw <- grDevices::terrain.colors(100)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else if (palette == "topo") {
      pal.raw <- grDevices::topo.colors(100)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else if (palette == "cm") {
      pal.raw <- grDevices::cm.colors(100)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else if (palette == "rainbow") {
      pal.raw <- grDevices::rainbow(100)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else if (toupper(palette) %in% names(ggsci_db)) {
      pal <- pal.def <- as.character(ggsci_db[[toupper(palette)]][[1]])
    } else if (tolower(palette) %in% names(ggsci_db)) {
      pal <- pal.def <- as.character(ggsci_db[[tolower(palette)]][[1]])
    } else if (tolower(palette) %in% tolower(names(ggsci_db))) {
      index <- which(tolower(names(ggsci_db)) == tolower(palette))
      pal <- pal.def <- as.character(ggsci_db[[index]][[1]])
    } else if (palette %in% grDevices::colors()[!grepl(pattern = "\\d", x = grDevices::colors())]) {
      fill <- palette
      # fill <- rep_len(palette, length.out = length(x))
      # names(fill) <- x
      return(fill)
    } else if (palette %in% grDevices::hcl.pals()) {
      pal.raw <- hcl.colors(100, palette = palette)
      pal <- pal.raw[round(seq.int(1, 100, length.out = length(unique(x))))]
      pal.def <- pal
    } else {
      if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
        pal.def <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[palette, "maxcolors"], palette)
        pal <- colorRampPalette(pal.def)(length(unique(x)))
      } else {
        if (all(strsplit(palette, split = ",")[[1]] %in% c(grDevices::colors(), paste0("#", c(0:9, letters[1:6]))))) {
          pal <- pal.def <- strsplit(palette, split = ",")[[1]]
        } else {
          warning("Cannot find the palette: ", palette, immediate. = TRUE)
          pal <- pal.def <- palette_default
        }
      }
    }
  }
  if (length(unique(x)) > length(pal)) {
    pal <- colorRampPalette(pal.def)(length(unique(x)))
  }
  if (isTRUE(reverse)) {
    pal <- rev(pal)
  }
  if (isTRUE(matched)) {
    fill <- pal[as.numeric(as.factor(x))]
    names(fill) <- x
  } else {
    if (inherits(x, "factor") && length(levels(x)) > 0) {
      fill <- pal[seq_along(levels(x))]
      names(fill) <- levels(x)
    } else {
      fill <- pal[seq_along(unique(x))]
      names(fill) <- unique(x)[order(unique(x))]
    }
  }
  return(fill)
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  } else {
    return(x)
  }
}

# check_DataType function has been moved to SCP-workflow.R to avoid duplication
# and improve V5 compatibility


as_matrix <- function(x, ...) {
  UseMethod(generic = "as_matrix", object = x)
}

#' @export
as_matrix.default <- function(x) {
  return(as.matrix(x = x))
}

#' @export
as_matrix.Matrix <- function(x) {
  # Force conversion for all Matrix-derived classes including V5 sparse matrices
  tryCatch({
    return(as.matrix(x = x))
  }, error = function(e) {
    # Fallback: try to convert via as() function for S4 classes
    return(as(x, "matrix"))
  })
}

#' @export
as_matrix.matrix <- function(x) {
  return(x)
}

#' @export
as_matrix.data.frame <- function(x) {
  return(as.matrix(x = x))
}

capitalize <- function(x, force_tolower = FALSE) {
  if (force_tolower) {
    x <- tolower(x)
  }
  gsub("(?<!^)(?=[A-Z])", " ", gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", tolower(x), perl = TRUE), perl = TRUE)
}

layer <- function(object, name = NULL) {
  if (!is.null(name)) {
    return(slot(object, name))
  }
  return(SeuratObject::Layers(object))
}

`layer<-` <- function(object, name = NULL, value) {
  if (!is.null(name)) {
    slot(object, name) <- value
  }
  return(object)
}

show_colors <- function(pal, label = TRUE) {
  n <- length(pal)
  if (n > 50) {
    show_palettes(palettes = split(pal, cut(1:n, ceiling(n / 50))), label = label, palettes_label = FALSE)
  } else {
    old <- graphics::par(mar = c(0.1, 0.1, 0.1, 0.1), bg = "white")
    on.exit(graphics::par(old))

    graphics::plot(1, 1,
      xlim = c(0, 1), ylim = c(0, 1), asp = 1,
      bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = ""
    )
    show_step <- 1 / n
    show_rect <- show_step * 0.85
    for (i in seq_along(pal)) {
      graphics::rect(
        xleft = (i - 1) * show_step, ybottom = 0, xright = (i - 1) * show_step + show_rect, ytop = 1,
        col = pal[i], lwd = 0.2, border = "black"
      )
      if (isTRUE(label)) {
        if (!is.null(names(pal)) && all(names(pal) != "")) {
          graphics::text((i - 1) * show_step + show_rect / 2, 0.5, names(pal)[i], srt = 90, adj = 0.5, cex = 0.8)
        } else {
          graphics::text((i - 1) * show_step + show_rect / 2, 0.5, pal[i], srt = 90, adj = 0.5, cex = 0.8)
        }
      }
    }
  }
}

show_palettes <- function(palettes = NULL, label = TRUE, palettes_label = TRUE) {
  if (is.null(palettes)) {
    palettes <- list(
      "palette_default" = palette_default,
      "palette_category" = palette_category,
      "palette_dimplot" = palette_dimplot,
      "ggplot_default" = rev(scales::hue_pal()(20)),
      "grey" = colorRampPalette(c("grey80", "grey20"))(20),
      "hue" = colorRampPalette(grDevices::hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[seq_len(3)])(20)
    )
    ggsci_db <- utils::getFromNamespace("ggsci_db", "ggsci")
    for (i in seq_along(ggsci_db)) {
      nm <- names(ggsci_db)[i]
      nm <- capitalize(tolower(nm))
      if (length(ggsci_db[[i]][[1]]) > 1) {
        palettes[[nm]] <- ggsci_db[[i]][[1]]
      }
    }
    for (i in unique(rownames(RColorBrewer::brewer.pal.info)[RColorBrewer::brewer.pal.info[, "colorblind"] == TRUE])) {
      palettes[[i]] <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[i, "maxcolors"], i)
    }
    for (i in c("cividis", "inferno", "magma", "plasma", "viridis", "mako", "rocket", "turbo")) {
      if (i %in% grDevices::hcl.pals()) {
        palettes[[capitalize(i)]] <- hcl.colors(20, palette = i)
      }
    }
  }
  n <- length(palettes)
  if (isTRUE(palettes_label)) {
    cat("All available palettes:", names(palettes), sep = "\n")
  }
  old <- graphics::par(mfrow = c(min(n, 30), 1), mar = c(0.2, 4, 0.2, 0.2), bg = "white")
  on.exit(graphics::par(old))
  if (n > 30) {
    devAskNewPage(ask = TRUE)
  }
  for (i in seq_along(palettes)) {
    if (length(palettes[[i]]) > 100) {
      pal <- palettes[[i]][round(seq.int(1, length(palettes[[i]]), length.out = 100))]
    } else {
      pal <- palettes[[i]]
    }
    m <- length(pal)
    graphics::plot(1, 1,
      xlim = c(0, 1), ylim = c(0, 1), asp = 1,
      bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = ""
    )
    graphics::mtext(names(palettes)[i], side = 2, las = 1, line = 0.1, cex = 0.7)
    show_step <- 1 / m
    show_rect <- show_step * 0.8
    for (j in seq_along(pal)) {
      graphics::rect(
        xleft = (j - 1) * show_step, ybottom = 0.02, xright = (j - 1) * show_step + show_rect, ytop = 0.98,
        col = pal[j], lwd = 0.2, border = "white"
      )
      if (isTRUE(label)) {
        if (!is.null(names(pal)) && all(names(pal) != "")) {
          graphics::text((j - 1) * show_step + show_rect / 2, 0.5, names(pal)[j], srt = 90, adj = 0.5, cex = 0.5)
        } else {
          graphics::text((j - 1) * show_step + show_rect / 2, 0.5, pal[j], srt = 90, adj = 0.5, cex = 0.5)
        }
      }
    }
  }
}


getpalette <- function(n, random = FALSE) {
  palette <- palette_scp(n = n, palette = "set1")

  out <- palette[1:n]
  if (random) {
    out <- sample(out, n)
  }
  return(out)
}

panel_fix_single <- function(x) {
  stopifnot(inherits(x, "ggplot"))
  if (max(unlist(lapply(x$layers, function(l) nrow(l[["data"]]))), na.rm = TRUE) == 0) {
    stop("No data information was found in the ggplot object.")
  }
  plot_data <- as.data.frame(unlist(lapply(x$layers, function(l) l[["data"]])))
  facet <- unlist(x$facet[["params"]][c("facets", "rows", "cols")])
  if (length(facet) == 0) {
    message("No facets detected.")
    return(x)
  }
  facet_levels <- unique(plot_data[[names(facet)]])
  empty_levels <- setdiff(facet_levels, as.character(unique(x$data[[names(facet)]])))
  if (length(empty_levels) > 0) {
    message("Adding panel for missing levels:", paste0(empty_levels, collapse = ","))
    df_empty <- x$data[rep(1, length(empty_levels)), , drop = FALSE]
    df_empty[[names(facet)]] <- empty_levels
    x$data <- rbind(x$data, df_empty)
  }
  return(x)
}

panel_fix <- function(x) {
  if (is.list(x) && !inherits(x, "ggplot")) {
    for (i in seq_along(x)) {
      x[[i]] <- panel_fix_single(x[[i]])
    }
  } else {
    x <- panel_fix_single(x)
  }
  return(x)
}

replicate_n <- function(x, n) {
  if (inherits(x, "list")) {
    x_out <- list()
    for (l in seq_along(x)) {
      x_out[[l]] <- rep(x[l], times = n)
    }
    x_out <- unlist(x_out, recursive = FALSE)
    return(x_out)
  } else {
    return(rep(list(x), times = n))
  }
}

#' Invoke a function with a list of arguments
#'
#' This is a utility function to call a function with arguments provided as a list.
#' Similar to do.call but with more flexibility for parameter handling.
#'
#' @param .fn Function to call, can be a function object or function name as string
#' @param .args List of arguments to pass to the function
#' @param ... Additional arguments passed to the function
#' @return Result of the function call
#' @export
invoke <- function(.fn, .args = list(), ...) {
  if (is.character(.fn)) {
    .fn <- get(.fn)
  }
  args <- modifyList(.args, list(...))
  do.call(.fn, args)
}


SCP_present <- function() {
  envs <- conda_list()
  SCP_envpath <- ""
  if ("SCP_env" %in% envs[["name"]]) {
    SCP_envpath <- envs[envs[["name"]] == "SCP_env", ][["python"]]
  }
  SCP_envpresent <- SCP_envpath != ""
  return(SCP_envpresent)
}

#' Check if UV is installed
#'
#' @return Logical, TRUE if UV is available
#' @export
check_uv <- function() {
  tryCatch({
    result <- system2("uv", "--version", stdout = TRUE, stderr = FALSE)
    return(length(result) > 0 && grepl("uv", result[1], ignore.case = TRUE))
  }, error = function(e) {
    return(FALSE)
  }, warning = function(w) {
    return(FALSE)
  })
}

#' Install UV if not present
#'
#' @param force Force reinstallation even if UV exists
#' @export
install_uv <- function(force = FALSE) {
  if (!check_uv() || force) {
    message("Installing UV...")

    # Platform-specific installation
    if (Sys.info()["sysname"] == "Windows") {
      # Windows installation using PowerShell
      system2("powershell",
              args = c("-ExecutionPolicy", "ByPass", "-c",
                      "irm https://astral.sh/uv/install.ps1 | iex"),
              wait = TRUE)
    } else {
      # Unix-like systems (macOS, Linux)
      system2("sh",
              args = c("-c", "curl -LsSf https://astral.sh/uv/install.sh | sh"),
              wait = TRUE)
    }

    # Check if installation was successful
    if (!check_uv()) {
      stop("UV installation failed. Please install manually from https://docs.astral.sh/uv/")
    }

    message("UV installed successfully!")
  } else {
    message("UV is already installed.")
  }
}

#' Check if UV environment exists
#'
#' @return Logical, TRUE if .venv exists in package directory
#' @export
uv_env_exists <- function() {
  # Check for .venv in the package root
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    # During development, use current directory
    pkg_dir <- getwd()
  }
  venv_path <- file.path(pkg_dir, ".venv")
  return(dir.exists(venv_path))
}

#' Create UV virtual environment
#'
#' @param python_version Python version to use (e.g., "3.10")
#' @export
uv_create_env <- function(python_version = "3.10") {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    # During development, use current directory
    pkg_dir <- getwd()
  }

  # Change to package directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(pkg_dir)

  # Copy .python-version file if it exists in inst directory
  python_version_locations <- c(
    ".python-version",  # Current directory
    file.path("inst", ".python-version"),  # In inst directory
    system.file(".python-version", package = "SCP")  # Installed location
  )

  for (loc in python_version_locations) {
    if (file.exists(loc) && loc != ".python-version") {
      file.copy(loc, ".python-version", overwrite = TRUE)
      break
    }
  }

  # Create virtual environment
  system2("uv", args = c("venv", "--python", python_version), wait = TRUE)

  # Verify creation
  if (!uv_env_exists()) {
    stop("Failed to create UV virtual environment")
  }

  return(invisible(TRUE))
}

#' Remove UV virtual environment
#'
#' @export
uv_remove_env <- function() {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    pkg_dir <- getwd()
  }
  venv_path <- file.path(pkg_dir, ".venv")

  if (dir.exists(venv_path)) {
    unlink(venv_path, recursive = TRUE)
    message("UV environment removed.")
  }
}

#' Sync Python dependencies using UV
#'
#' @param extras Character vector of extra dependency groups to install
#' @export
uv_sync_deps <- function(extras = "all") {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    pkg_dir <- getwd()
  }

  # Change to package directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(pkg_dir)

  # Look for pyproject.toml in multiple locations
  pyproject_locations <- c(
    "pyproject.toml",  # Current directory
    file.path("inst", "pyproject.toml"),  # In inst directory
    system.file("pyproject.toml", package = "SCP")  # Installed location
  )

  pyproject_path <- NULL
  for (loc in pyproject_locations) {
    if (file.exists(loc)) {
      pyproject_path <- loc
      break
    }
  }

  if (is.null(pyproject_path)) {
    stop("pyproject.toml not found. Cannot sync dependencies.")
  }

  # Copy pyproject.toml to current directory if needed
  if (pyproject_path != "pyproject.toml") {
    file.copy(pyproject_path, "pyproject.toml", overwrite = TRUE)
  }

  # Prepare extras argument
  if (length(extras) == 1 && extras == "all") {
    extras_arg <- "--all-extras"
  } else if (length(extras) > 0) {
    extras_arg <- paste0("--extra ", extras, collapse = " ")
  } else {
    extras_arg <- ""
  }

  # Install dependencies
  # First, compile requirements from pyproject.toml
  if (length(extras) == 0 || (length(extras) == 1 && extras == "")) {
    # Install only core dependencies
    compile_cmd <- c("pip", "compile", "pyproject.toml", "-o", "requirements.txt", "--quiet")
  } else if (extras_arg == "--all-extras") {
    # Install with all extras
    compile_cmd <- c("pip", "compile", "pyproject.toml", "--all-extras", "-o", "requirements.txt", "--quiet")
  } else {
    # Install specific extras
    extras_args <- unlist(lapply(extras, function(e) c("--extra", e)))
    compile_cmd <- c("pip", "compile", "pyproject.toml", extras_args, "-o", "requirements.txt", "--quiet")
  }

  # Compile the requirements
  system2("uv", args = compile_cmd, wait = TRUE, stderr = FALSE, stdout = FALSE)

  # Install from compiled requirements
  if (file.exists("requirements.txt")) {
    result <- system2("uv", args = c("pip", "install", "-r", "requirements.txt"), wait = TRUE)
  } else {
    # Fallback: install packages directly
    packages <- c("numpy", "pandas", "scipy", "matplotlib", "seaborn",
                 "scikit-learn", "h5py", "numba", "anndata")
    result <- system2("uv", args = c("pip", "install", packages), wait = TRUE)
  }

  if (result != 0) {
    warning("UV installation may have encountered issues. Checking installation...")
  }

  # Special handling for Apple Silicon
  if (Sys.info()["sysname"] == "Darwin" && grepl("arm64", Sys.info()["machine"], ignore.case = TRUE)) {
    if ("all" %in% extras || "deeplearning" %in% extras) {
      message("Installing Apple Silicon optimizations...")
      system2("uv", args = c("pip", "install", "-e", ".[apple_silicon]"), wait = TRUE)
    }
  }

  return(invisible(TRUE))
}

#' Configure reticulate to use UV environment
#'
#' @export
use_uv_env <- function() {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    pkg_dir <- getwd()
  }

  venv_path <- file.path(pkg_dir, ".venv")
  if (!dir.exists(venv_path)) {
    stop("UV environment not found. Run PrepareEnv(method='uv') first.")
  }

  # Find Python executable in venv
  if (Sys.info()["sysname"] == "Windows") {
    python_path <- file.path(venv_path, "Scripts", "python.exe")
  } else {
    python_path <- file.path(venv_path, "bin", "python")
  }

  if (!file.exists(python_path)) {
    stop("Python executable not found in UV environment")
  }

  # Use the UV environment
  reticulate::use_python(python_path, required = TRUE)

  return(invisible(TRUE))
}

#' Install the SCP python environment
#'
#' Install all python packages in the SCP environment automatically. The environment will be created
#' using UV (if available) or Conda as a fallback.
#'
#' @param method Installation method. Options are "auto", "uv", or "conda". "auto" will prefer UV if available, otherwise use conda.
#' @param pip Whether to use pip for installing packages. This is only relevant when method is "conda".
#' @param user Whether to install packages into a user site library, instead of the default system site library. Note that this argument is only relevant when using the Conda package manager.
#' @param force Force reinstall the SCP python environment.
#' @param update Whether to update packages that already exist within the environment. For UV, this always syncs to the latest allowed versions. For Conda, options are "default", "all", and "selected".
#' @param python_version Python version to use for the SCP environment. Default is "3.10". Note that with UV, this is controlled by .python-version file.
#' @param extras Character vector of extra dependency groups to install (e.g., c("velocity", "trajectory")). Default is "all" for complete installation.
#' @export
#'
PrepareEnv <- function(method = "auto", pip = FALSE, user = FALSE, force = FALSE, update = "all", python_version = "3.10", extras = "all") {
  # Determine which method to use
  if (method == "auto") {
    if (check_uv()) {
      method <- "uv"
      message("UV detected. Using UV for Python environment management.")
    } else {
      method <- "conda"
      message("UV not found. Using conda for Python environment management.")
      message("To use UV (recommended for faster installations), install it with:")
      message("  curl -LsSf https://astral.sh/uv/install.sh | sh")
    }
  }

  # Check if we should use UV
  if (method == "uv") {
    if (!check_uv()) {
      stop("UV is not installed. Please install UV or use method='conda'")
    }

    # UV-based installation
    if (!uv_env_exists() || isTRUE(force)) {
      if (isTRUE(force) && uv_env_exists()) {
        message("Removing existing UV environment...")
        uv_remove_env()
      }

      message("Creating UV environment with Python ", python_version, "...")
      uv_create_env(python_version = python_version)

      message("Installing Python dependencies with UV...")
      uv_sync_deps(extras = extras)

      message("UV environment setup complete!")
    } else {
      if (isTRUE(update)) {
        message("Updating UV environment...")
        uv_sync_deps(extras = extras)
      } else {
        message("UV environment already exists. Use force=TRUE to reinstall or update=TRUE to update packages.")
      }
    }
  } else {
    # Fallback to original conda-based installation
    if (!SCP_present() || isTRUE(force)) {
      if (isTRUE(force)) {
        message("The SCP python environment will be reinstalled with Python ", python_version, ".")
      } else {
        message("The SCP python environment has not been created. Creating it now with Python ", python_version, ".")
      }
      packages <- c(
        "anndata", "pandas", "numpy", "scipy", "matplotlib", "seaborn", "scikit-learn",
        "scikit-misc", "pynndescent", "umap-learn", "pymde", "opentsne", "phate", "scanorama", "bbknn",
        "leidenalg", "louvain", "rapids-singlecell", "scrublet", "scvi-tools", "torch", "h5py", "numba", "pybind11", "xgboost"
      )

      # Check for Apple Silicon and modify scvi-tools installation
      if (Sys.info()["sysname"] == "Darwin" && grepl("arm64", Sys.info()["machine"], ignore.case = TRUE)) {
        message("Apple Silicon detected. Using optimized scvi-tools installation with Metal support.")
        # Remove scvi-tools from the regular package list
        packages <- packages[packages != "scvi-tools"]
      }
      install_py(packages = packages, method = "conda", pip = pip, user = user, force = force, update = update, python_version = python_version)
    } else {
      message("The SCP python environment is already present.")
    }
  }
}

install_py <- function(packages = NULL, method = "auto", pip = FALSE, user = FALSE, force = FALSE, update = "all", python_version = NULL, ...) {
  if (reticulate::is_osx() && reticulate::is_python_conda()) {
    nomkl <- TRUE
  } else {
    nomkl <- FALSE
  }
  packages <- unique(packages)
  condaenv <- "SCP_env"
  reticulate::conda_create(
    envname = condaenv,
    conda = "auto",
    python_version = python_version,
    force = force
  )
  for (package in packages) {
    message("Install package: ", package, "\n", sep = "")
    if (package == "scikit-misc") {
      system(paste0(reticulate::conda_python(envname = condaenv, conda = "auto"), " -m pip install scikit-misc==0.1.4"))
    } else if (package == "anndata==0.8.0") {
      reticulate::conda_install(envname = condaenv, packages = package, pip = TRUE, ...)
    } else {
      reticulate::conda_install(envname = condaenv, packages = package, ...)
    }
  }
  
  # Handle Apple Silicon scvi-tools installation separately
  if (Sys.info()["sysname"] == "Darwin" && grepl("arm64", Sys.info()["machine"], ignore.case = TRUE)) {
    message("Installing scvi-tools with Metal support for Apple Silicon...")
    system(paste0(reticulate::conda_python(envname = condaenv, conda = "auto"), " -m pip install -U 'scvi-tools[metal]'"))
    # Apply macOS compatibility fix
    Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  }
  message("Done.\n", sep = "")
}

try_get <- function(expr, max_tries = 3, error_message = NULL, default = NULL) {
  result <- NULL
  for (i in 1:max_tries) {
    result <- tryCatch(expr, error = function(e) {
      if (i < max_tries) {
        message("Attempt ", i, " failed. Retrying...")
        Sys.sleep(1)
      } else if (!is.null(error_message)) {
        message(error_message)
      }
      return(NULL)
    })
    if (!is.null(result)) {
      return(result)
    }
  }
  return(default)
}

wait_for_file <- function(file_path, interval = 0.1, timeout = 5) {
  elapsed <- 0
  while (!file.exists(file_path) && elapsed < timeout) {
    Sys.sleep(interval)
    elapsed <- elapsed + interval
  }
  return(file.exists(file_path))
}

strwrap2 <- function(x = NULL, width = 50, whitespace_only = TRUE) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(lapply(x, function(i) paste0(strwrap(i, width = width), collapse = "\n")))
  return(x_wrap)
}

tochunks <- function(x, nchunks) {
  split(x, cut(seq_along(x), nchunks, labels = FALSE))
}

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

layerNames <- function(srt, assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  assay <- assay %||% DefaultAssay(srt)
  assay_obj <- srt[[assay]]

  # Handle Seurat V5 vs V4 differences
  if (IsSeurat5(srt)) {
    # For Seurat V5, use Layers() function
    if (exists("Layers", mode = "function")) {
      return(Layers(assay_obj))
    } else {
      # Fallback for V5 if Layers() doesn't exist
      if (methods::hasMethod("names", class(assay_obj@layers))) {
        return(names(assay_obj@layers))
      }
    }
  } else {
    # For Seurat V4, layer names are the slot names
    layer_names <- c()
    if (length(slot(assay_obj, "counts")) > 0) layer_names <- c(layer_names, "counts")
    if (length(slot(assay_obj, "data")) > 0) layer_names <- c(layer_names, "data")
    if (length(slot(assay_obj, "scale.data")) > 0) layer_names <- c(layer_names, "scale.data")
    return(layer_names)
  }

  # Final fallback
  return(c("counts", "data", "scale.data"))
}

#' Get data from Seurat object in a version-agnostic way
#'
#' This function extracts data from a Seurat object, handling both V4 and V5 formats.
#'
#' @param srt A Seurat object
#' @param layer The layer to extract (e.g., "counts", "data", "scale.data")
#' @param assay Name of the assay. If NULL, uses the default assay.
#'
#' @return A matrix of the requested data
#' @export
get_seurat_data <- function(srt, layer = "data", assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  is_v5 <- IsSeurat5(srt)
  assay <- assay %||% DefaultAssay(srt)

  if (is_v5) {
    return(LayerData(srt[[assay]], layer = layer))
  } else {
    return(GetAssayData(srt, layer = layer, assay = assay))
  }
}

#' Set data in Seurat object in a version-agnostic way
#'
#' This function sets data in a Seurat object, handling both V4 and V5 formats.
#'
#' @param srt A Seurat object
#' @param data The data to set
#' @param layer The layer to set (e.g., "counts", "data", "scale.data")
#' @param assay Name of the assay. If NULL, uses the default assay.
#'
#' @return The modified Seurat object
#' @export
set_seurat_data <- function(srt, data, layer = "data", assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  is_v5 <- IsSeurat5(srt)
  assay <- assay %||% DefaultAssay(srt)

  if (is_v5) {
    LayerData(srt[[assay]], layer = layer) <- data
  } else {
    srt <- SetAssayData(srt, layer = layer, new.data = data, assay = assay)
  }

  return(srt)
}

#' Get feature metadata from Seurat object in a version-agnostic way
#'
#' This function extracts feature metadata from a Seurat object, handling both V4 and V5 formats.
#'
#' @param srt A Seurat object
#' @param assay Name of the assay. If NULL, uses the default assay.
#'
#' @return A data frame of feature metadata
#' @export
get_feature_metadata <- function(srt, assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  is_v5 <- IsSeurat5(srt)
  assay <- assay %||% DefaultAssay(srt)
  assay_obj <- srt[[assay]]

  if (is_v5) {
    # For Seurat V5
    if (exists("Features", mode = "function")) {
      # Get feature metadata using V5 functions
      feature_df <- Features(assay_obj, layer = NULL)
      if (is.character(feature_df)) {
        # If it returns feature names, create a data frame
        feature_df <- data.frame(row.names = feature_df)
      }
      return(feature_df)
    } else {
      # Fallback: access meta.features slot if it exists
      if (.hasSlot(assay_obj, "meta.features")) {
        return(slot(assay_obj, "meta.features"))
      } else {
        # Return empty data frame with correct row names
        return(data.frame(row.names = rownames(assay_obj)))
      }
    }
  } else {
    # For Seurat V4
    if (.hasSlot(assay_obj, "meta.features")) {
      return(slot(assay_obj, "meta.features"))
    } else {
      # Return empty data frame with correct row names
      return(data.frame(row.names = rownames(assay_obj)))
    }
  }
}

#' Set feature metadata in Seurat object in a version-agnostic way
#'
#' This function sets feature metadata in a Seurat object, handling both V4 and V5 formats.
#'
#' @param srt A Seurat object
#' @param metadata A data frame of feature metadata
#' @param assay Name of the assay. If NULL, uses the default assay.
#'
#' @return The modified Seurat object
#' @export
set_feature_metadata <- function(srt, metadata, assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  is_v5 <- IsSeurat5(srt)
  assay <- assay %||% DefaultAssay(srt)
  assay_obj <- srt[[assay]]

  if (is_v5) {
    # For Seurat V5, we need to be careful about how we set metadata
    # The recommended way is to add columns to the existing metadata
    existing_meta <- get_feature_metadata(srt, assay)

    # Merge with existing metadata
    if (ncol(existing_meta) > 0) {
      # Keep existing columns and add new ones
      for (col in colnames(metadata)) {
        existing_meta[[col]] <- metadata[[col]][match(rownames(existing_meta), rownames(metadata))]
      }
      metadata <- existing_meta
    }

    # Set the metadata back
    if (.hasSlot(assay_obj, "meta.features")) {
      slot(assay_obj, "meta.features") <- metadata
      srt[[assay]] <- assay_obj
    } else {
      # For V5, we might need to use a different approach
      warning("Unable to set feature metadata for Seurat V5. This might require manual intervention.")
    }
  } else {
    # For Seurat V4
    if (.hasSlot(assay_obj, "meta.features")) {
      slot(assay_obj, "meta.features") <- metadata
      srt[[assay]] <- assay_obj
    }
  }

  return(srt)
}

#' Get variable features from Seurat object in a version-agnostic way
#'
#' This function extracts variable features from a Seurat object, handling both V4 and V5 formats.
#'
#' @param srt A Seurat object
#' @param assay Name of the assay. If NULL, uses the default assay.
#'
#' @return A character vector of variable feature names
#' @export
get_var_features <- function(srt, assay = NULL) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  is_v5 <- IsSeurat5(srt)
  assay <- assay %||% DefaultAssay(srt)
  assay_obj <- srt[[assay]]

  if (is_v5) {
    # For Seurat V5, use VariableFeatures function
    return(VariableFeatures(srt, assay = assay))
  } else {
    # For Seurat V4, access the slot directly
    if (.hasSlot(assay_obj, "var.features")) {
      return(slot(assay_obj, "var.features"))
    } else {
      return(character(0))
    }
  }
}

#' Get Seurat version information
#'
#' Returns detailed version information about the Seurat object and installed Seurat package.
#'
#' @param srt A Seurat object (optional). If provided, includes object-specific version info.
#' @return A list with version information
#' @export
GetSeuratVersion <- function(srt = NULL) {
  result <- list(
    package_version = packageVersion("Seurat"),
    package_major = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1]),
    is_v5_installed = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1]) >= 5
  )
  
  if (!is.null(srt)) {
    if (!inherits(srt, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    # Object-specific version information
    result$object_version <- attr(srt, "version")
    result$object_class <- class(srt)
    
    # Check default assay information
    tryCatch({
      default_assay <- DefaultAssay(srt)
      assay_obj <- srt[[default_assay]]
      result$assay_class <- class(assay_obj)
      result$has_layers <- !is.null(tryCatch(Layers(assay_obj), error = function(e) NULL))
      result$is_object_v5 <- IsSeurat5(srt)
    }, error = function(e) {
      result$assay_class <- "unknown"
      result$has_layers <- FALSE
      result$is_object_v5 <- FALSE
    })
  }
  
  return(result)
}

#' Check if object is Seurat V5
#'
#' @param srt A Seurat object
#' @return Logical indicating if the object is Seurat V5
#' @export
IsSeurat5 <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Create a comprehensive check for Seurat V5 features
  tryCatch({
    # Method 1: Check object version attribute (most reliable)
    seurat_version <- attr(srt, "version")
    if (!is.null(seurat_version)) {
      # Handle different version formats: "5.0.0", "5", etc.
      version_str <- as.character(seurat_version)
      if (grepl("^[0-9]+", version_str)) {
        major_version <- as.numeric(strsplit(version_str, "\\.")[[1]][1])
        if (!is.na(major_version) && major_version >= 5) {
          return(TRUE)
        }
      }
    }
    
    # Method 2: Check if the assay class is Assay5
    default_assay <- tryCatch(DefaultAssay(srt), error = function(e) NULL)
    if (!is.null(default_assay)) {
      assay <- tryCatch(srt[[default_assay]], error = function(e) NULL)
      if (!is.null(assay) && inherits(assay, "Assay5")) {
        return(TRUE)
      }
    }
    
    # Method 3: Check for presence of layer-based structure (V5 specific)
    if (!is.null(default_assay)) {
      assay <- tryCatch(srt[[default_assay]], error = function(e) NULL)
      if (!is.null(assay)) {
        # Check if Layers function exists and returns valid data
        layers_result <- tryCatch(Layers(assay), error = function(e) NULL)
        if (!is.null(layers_result) && length(layers_result) > 0) {
          # Additional check: V5 objects should have layer-specific structure
          if (is.character(layers_result) || is.list(layers_result)) {
            return(TRUE)
          }
        }
        
        # Method 4: Check for V5-specific slots or methods
        if (.hasSlot(assay, "layers")) {
          layers_slot <- tryCatch(slot(assay, "layers"), error = function(e) NULL)
          if (!is.null(layers_slot) && length(layers_slot) > 0) {
            return(TRUE)
          }
        }
      }
    }
    
    # Method 5: Check package version as final fallback
    seurat_pkg_version <- tryCatch(packageVersion("Seurat"), error = function(e) NULL)
    if (!is.null(seurat_pkg_version)) {
      major_version <- as.numeric(strsplit(as.character(seurat_pkg_version), "\\.")[[1]][1])
      if (!is.na(major_version) && major_version >= 5) {
        # If we have V5 installed but object doesn't have V5 features, 
        # it might be an old object that needs updating
        # Return FALSE here to indicate it's not a V5 object yet
      }
    }
    
    # If none of the above checks pass, it's not a V5 object
    return(FALSE)
    
  }, error = function(e) {
    # If any errors occur during the checks, be conservative and return FALSE
    warning("Error during Seurat version detection: ", e$message, call. = FALSE)
    return(FALSE)
  })
}

#' Ensure Seurat object is V5-compatible
#'
#' This function checks if a Seurat object is V5-compatible and converts it if needed.
#' It requires Seurat V5 to be installed. If the object is already V5, it is returned
#' unchanged.
#'
#' @param srt A Seurat object
#' @param verbose Whether to show verbose output
#' @return A Seurat V5 object
#' 
#' @importFrom Seurat UpdateSeuratObject
#' @importFrom utils packageVersion
#' @export
EnsureSeurat5 <- function(srt, verbose = TRUE) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Check if we have Seurat V5 installed
  tryCatch({
    seurat_version <- packageVersion("Seurat")
    # Extract the major version number properly
    major_version <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1])
    is_v5_installed <- major_version >= 5
    
    if (!is_v5_installed) {
      stop("Seurat V5 is required for this function. Please install Seurat V5 with:\n",
           "  install.packages('Seurat')")
    }
  }, error = function(e) {
    stop("Failed to check Seurat version: ", e$message)
  })
  
  # Check if the object is already V5
  if (IsSeurat5(srt)) {
    if (verbose) message("Object is already Seurat V5")
    return(srt)
  }
  
  # Convert the object to V5
  if (verbose) message("Converting Seurat object from V4 to V5...")
  
  # UpdateSeuratObject is required for the conversion
  tryCatch({
    srt_v5 <- UpdateSeuratObject(srt)
    if (verbose) message("Conversion successful")
    return(srt_v5)
  }, error = function(e) {
    stop("Failed to convert Seurat object to V5: ", e$message, 
         "\nPlease use Seurat::UpdateSeuratObject() directly to resolve any issues.")
  })
}

#' Remove environment helper functions
#'
#' @keywords internal
format_size <- function(size_bytes) {
  if (size_bytes < 1024) {
    return(paste0(size_bytes, " B"))
  } else if (size_bytes < 1024^2) {
    return(paste0(round(size_bytes / 1024, 1), " KB"))
  } else if (size_bytes < 1024^3) {
    return(paste0(round(size_bytes / 1024^2, 1), " MB"))
  } else {
    return(paste0(round(size_bytes / 1024^3, 1), " GB"))
  }
}

#' @keywords internal
dir_size <- function(path) {
  if (!dir.exists(path)) {
    return(0)
  }

  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  sum(file.info(files)$size, na.rm = TRUE)
}

#' Remove SCP Python environment
#'
#' This function removes the SCP Python environment from your system.
#' It will prompt for confirmation before deletion unless `prompt = FALSE`.
#'
#' @param envname The name of the conda environment to remove. Default is "SCP_env".
#' @param prompt Whether to prompt for confirmation before removing. Default is TRUE.
#' @param conda The path to conda executable. Default is "auto" which automatically finds conda.
#'
#' @return NULL (invisibly). Messages will indicate success or failure.
#'
#' @examples
#' \dontrun{
#' # Remove with confirmation prompt
#' RemoveEnv()
#'
#' # Remove without prompt
#' RemoveEnv(prompt = FALSE)
#' }
#'
#' @export
RemoveEnv <- function(envname = "SCP_env", prompt = TRUE, conda = "auto") {
  # Check if environment exists
  envs <- tryCatch(
    {
      reticulate::conda_list(conda = conda)
    },
    error = function(e) {
      message("Could not find conda installation. Is conda installed?")
      return(data.frame(name = character(), python = character()))
    }
  )

  if (!envname %in% envs$name) {
    message("Environment '", envname, "' not found.")
    return(invisible(NULL))
  }

  # Get conda path and environment directory
  conda_binary <- reticulate::conda_binary(conda)
  envs_dir <- dirname(dirname(envs$python[envs$name == envname]))
  if (basename(envs_dir) != "envs") {
    # Try to find the envs directory
    conda_info <- system2(conda_binary, c("info", "--json"), stdout = TRUE)
    conda_info <- jsonlite::fromJSON(paste(conda_info, collapse = ""))
    envs_dir <- conda_info$envs_dirs[1]
  }

  env_path <- paste0(envs_dir, "/", envname)

  # Get environment information
  env_info <- tryCatch(
    {
      python_path <- conda_python(conda = conda, envname = envname)
      python_version <- reticulate:::python_version(python_path)
      list(path = env_path, python = python_path, version = python_version)
    },
    error = function(e) {
      list(path = env_path, python = NA, version = NA)
    }
  )

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
      title = "Are you sure you want to remove this environment?"
    ) == 1
  }

  if (proceed) {
    message("Removing environment '", envname, "'...")

    # Try using conda remove first (cleaner)
    success <- tryCatch(
      {
        result <- system2(conda, c("env", "remove", "-n", envname, "--yes"), stdout = TRUE, stderr = TRUE)
        # Check if the command was successful
        !any(grepl("error", result, ignore.case = TRUE))
      },
      error = function(e) FALSE
    )

    # If conda remove failed, try directly removing the directory
    if (!success) {
      message("Conda removal failed, attempting direct directory removal...")
      success <- tryCatch(
        {
          unlink(env_path, recursive = TRUE)
          !file.exists(env_path)
        },
        error = function(e) FALSE
      )
    }

    if (success) {
      message("[v] Successfully removed environment '", envname, "'")
      # Unset Python environment variables to avoid pointing to deleted environment
      Sys.unsetenv("RETICULATE_PYTHON")
      Sys.unsetenv("RETICULATE_PYTHON_ENV")

      # Clear reticulate's python configuration cache
      if (exists(".globals", envir = asNamespace("reticulate"))) {
        .globals <- get(".globals", envir = asNamespace("reticulate"))
        .globals$py_config <- NULL
      }
    } else {
      message("[x] Failed to remove environment '", envname, "'")
      message("You may need to manually delete: ", env_path)
    }
  } else {
    message("Environment removal cancelled.")
  }

  return(invisible(NULL))
}

#' Test Python Compatibility
#'
#' Test compatibility with different Python versions for SCP environment setup.
#' This function tests the specified Python versions on different platforms.
#'
#' @param versions Character vector of Python versions to test. Default is c("3.10", "3.11", "3.12").
#' @param test_packages Logical. Whether to test package installation. Default is TRUE.
#' @export
#'
TestPythonCompatibility <- function(versions = c("3.10", "3.11", "3.12"), test_packages = TRUE) {
  results <- list()
  
  for (version in versions) {
    cat("Testing Python", version, "compatibility...\n")
    
    # Test basic version support
    version_result <- list(
      version = version,
      platform = Sys.info()[["sysname"]],
      conda_available = FALSE,
      packages_installable = FALSE
    )
    
    # Test conda availability
    tryCatch({
      conda_path <- reticulate::conda_binary()
      if (!is.null(conda_path) && file.exists(conda_path)) {
        version_result$conda_available <- TRUE
      }
    }, error = function(e) {
      version_result$conda_available <- FALSE
    })
    
    # Test package installation (if requested and conda is available)
    if (test_packages && version_result$conda_available) {
      tryCatch({
        test_env <- paste0("SCP_test_", gsub("\\.", "", version))
        
        # Try to create a test environment
        reticulate::conda_create(
          envname = test_env,
          python_version = version,
          packages = "numpy",  # Test with a simple package
          forge = TRUE
        )
        
        version_result$packages_installable <- TRUE
        
        # Clean up test environment
        reticulate::conda_remove(envname = test_env)
        
      }, error = function(e) {
        version_result$packages_installable <- FALSE
        cat("  Warning: Package installation test failed for Python", version, "\n")
      })
    }
    
    results[[version]] <- version_result
    
    cat("  Platform:", version_result$platform, "\n")
    cat("  Conda available:", version_result$conda_available, "\n")
    if (test_packages) {
      cat("  Packages installable:", version_result$packages_installable, "\n")
    }
    cat("\n")
  }
  
  return(invisible(results))
}
