#' @importFrom BiocParallel bpparam
#' @importFrom Matrix t
#' @importFrom reticulate py_module_available use_python import_from_path py_available
#' @importFrom stats as.dist as.formula median model.matrix prcomp sd var aggregate kmeans na.omit setNames quantile weighted.mean
#' @importFrom utils askYesNo head installed.packages menu modifyList packageVersion read.csv read.table setTxtProgressBar tail txtProgressBar download.file
#' @importFrom grDevices col2rgb colorRampPalette dev.cur dev.interactive dev.new dev.off devAskNewPage palette pdf recordPlot rgb
#' @importFrom SeuratObject DefaultAssay CreateAssayObject CreateDimReducObject
#' @importFrom SeuratObject JoinLayers Layers
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

# Removed check_Python() function - exported but never used
# Removed ensure_scp_python() function - internal helper that was never called
# Removed EnsureEnv() function - unnecessary wrapper
# Use use_uv_env() directly when Python environment is needed

#' Unified Color Palette System
#'
#' Consolidated color palette functions for SCPNext package
#'
#' @param type Character, palette type: "default", "category", "dimplot", or "custom"
#' @param n Integer, number of colors to return
#' @param custom_colors Character vector, custom color palette
#' @return Character vector of color codes
#' @export
GetPalette <- function(type = "default", n = NULL, custom_colors = NULL) {
  
  # Define base palettes
  palettes <- list(
    default = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
      "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
      "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
      "#B3B3B3", "#8DA0CB", "#FC8D62", "#66C2A5", "#E6F5C9", "#FFF2AE", "#F4CAE4", "#F1E2CC"
    ),
    category = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999",
      "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
      "#E5C494", "#E78AC3", "#A6D854", "#FFD92F", "#FC8D62", "#66C2A5", "#BC80BD", "#CCEBC5"
    ),
    dimplot = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
      "#E78AC3", "#A6D854", "#FFD92F", "#FC8D62", "#66C2A5", "#8DA0CB", "#E5C494",
      "#B3B3B3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
      "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#F4CAE4", "#F1E2CC", "#CCCCCC"
    )
  )
  
  # Validate type
  if (!type %in% c("default", "category", "dimplot", "custom")) {
    stop("'type' must be one of: 'default', 'category', 'dimplot', 'custom'")
  }
  
  # Get base palette
  if (type == "custom") {
    if (is.null(custom_colors)) {
      stop("'custom_colors' must be provided when type = 'custom'")
    }
    base_palette <- custom_colors
  } else {
    base_palette <- palettes[[type]]
  }
  
  # Adjust length if needed
  if (!is.null(n)) {
    if (n <= length(base_palette)) {
      palette <- base_palette[1:n]
    } else {
      palette <- grDevices::colorRampPalette(base_palette)(n)
    }
  } else {
    palette <- base_palette
  }
  
  return(palette)
}

# Backward compatibility: Keep old function names as thin wrappers
#' @rdname GetPalette
#' @export
palette_default <- function(n = NULL) {
  GetPalette(type = "default", n = n)
}

#' @rdname GetPalette
#' @export
palette_category <- function(n = NULL) {
  GetPalette(type = "category", n = n)
}

#' @rdname GetPalette
#' @export
palette_dimplot <- function(n = NULL) {
  GetPalette(type = "dimplot", n = n)
}

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

#' Convert color to hex format
#'
#' Replacement for gplots::col2hex using base R functions
#' @param colors Character vector of color names or specifications
#' @return Character vector of hex color codes
#' @keywords internal
col2hex <- function(colors) {
  rgb_matrix <- col2rgb(colors)
  rgb(rgb_matrix[1, ], rgb_matrix[2, ], rgb_matrix[3, ], maxColorValue = 255)
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

#' Require packages with helpful error messages
#'
#' @param packages Character vector of package names to require
#' @param context Optional context string for error message
#' @param quietly Logical, whether to suppress package startup messages
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

#' Invoke a function with arguments
#'
#' Wrapper around do.call that handles character function names
#'
#' @param .fn Function or character name of function to call
#' @param .args List of arguments to pass to function
#' @param ... Additional arguments (merged with .args)
#' @keywords internal
invoke <- function(.fn, .args = list(), ...) {
  if (is.character(.fn)) {
    .fn <- get(.fn)
  }
  args <- modifyList(.args, list(...))
  do.call(.fn, args)
}

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

#' Capitalize strings
#'
#' Capitalizes the first letter of each word in a string
#'
#' @param x Character vector to capitalize
#' @param force_tolower Logical, whether to convert to lowercase first
#' @return Character vector with capitalized strings
#' @export
capitalize <- function(x, force_tolower = FALSE) {
  if (force_tolower) {
    x <- tolower(x)
  }
  gsub("(?<!^)(?=[A-Z])", " ", gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", tolower(x), perl = TRUE), perl = TRUE)
}

layer <- function(object, name = NULL) {
  if (!is.null(name)) {
    # For Seurat v5 Assay5 objects, handle data/counts/scale.data as layers
    if (inherits(object, "Assay5") && name %in% c("counts", "data", "scale.data")) {
      # Use LayerData for actual data layers in v5
      # Check if layer exists
      if (name %in% SeuratObject::Layers(object)) {
        return(SeuratObject::LayerData(object, layer = name))
      } else {
        return(NULL)
      }
    }
    # For everything else (meta.features, var.features, etc.), use slot
    return(slot(object, name))
  }
  return(SeuratObject::Layers(object))
}

`layer<-` <- function(object, name = NULL, value) {
  if (!is.null(name)) {
    # For Seurat v5 Assay5 objects, handle data/counts/scale.data as layers
    if (inherits(object, "Assay5") && name %in% c("counts", "data", "scale.data")) {
      # Use LayerData<- for actual data layers in v5
      SeuratObject::LayerData(object, layer = name) <- value
      return(object)
    }
    # For everything else, use slot
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
      "palette_default" = GetPalette("default"),
      "palette_category" = GetPalette("category"),
      "palette_dimplot" = GetPalette("dimplot"),
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

# Removed quiet_ggplot2_scales() function - no longer needed as root cause of scale warnings
# has been fixed by properly configuring manual scales with drop=FALSE and explicit breaks

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


# Removed SCP_present() function - no longer needed with UV-only approach
# Removed invoke() function - unused wrapper around do.call()

# ============================================================================
# UV Environment Helper Functions
# ============================================================================

#' Get SCP package directory
#'
#' @return Character, path to SCP package directory
#' @keywords internal
get_scp_pkg_dir <- function() {
  pkg_dir <- system.file("", package = "SCPNext")
  if (pkg_dir == "") pkg_dir <- getwd()
  # Remove trailing slash to avoid double slashes in file.path()
  pkg_dir <- sub("/$", "", pkg_dir)
  pkg_dir
}

#' Find pyproject.toml file
#'
#' @return Character path to pyproject.toml or NULL if not found
#' @keywords internal
find_pyproject_toml <- function() {
  locations <- c(
    "pyproject.toml",
    file.path("inst", "pyproject.toml"),
    system.file("pyproject.toml", package = "SCPNext")
  )
  for (loc in locations) {
    if (file.exists(loc)) return(loc)
  }
  NULL
}

#' Get Python executable path in UV environment
#'
#' @param venv_path Path to venv directory. If NULL, uses default SCP .venv
#' @return Character path to Python executable
#' @keywords internal
get_uv_python_path <- function(venv_path = NULL) {
  if (is.null(venv_path)) {
    venv_path <- file.path(get_scp_pkg_dir(), ".venv")
  }
  if (Sys.info()["sysname"] == "Windows") {
    file.path(venv_path, "Scripts", "python.exe")
  } else {
    file.path(venv_path, "bin", "python")
  }
}

#' Extract major version number from version string
#'
#' @param version_string Version string (e.g., "5.0.1" or package_version object)
#' @return Numeric major version or NA if cannot parse
#' @keywords internal
extract_major_version <- function(version_string) {
  version_str <- as.character(version_string)
  if (grepl("^[0-9]+", version_str)) {
    major <- as.numeric(strsplit(version_str, "\\.")[[1]][1])
    if (!is.na(major)) return(major)
  }
  NA_real_
}

# ============================================================================
# UV Environment Management Functions
# ============================================================================

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
  venv_path <- file.path(get_scp_pkg_dir(), ".venv")
  dir.exists(venv_path)
}

#' Create UV virtual environment
#'
#' @param python_version Python version to use (e.g., "3.10")
#' @export
uv_create_env <- function(python_version = "3.10") {
  pkg_dir <- get_scp_pkg_dir()

  # Change to package directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(pkg_dir)

  # Copy .python-version file if it exists in inst directory
  python_version_locations <- c(
    ".python-version",  # Current directory
    file.path("inst", ".python-version"),  # In inst directory
    system.file(".python-version", package = "SCPNext")  # Installed location
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

  # Fix missing libpython symlink on macOS (common UV issue)
  if (Sys.info()["sysname"] == "Darwin") {
    venv_path <- file.path(pkg_dir, ".venv")

    # Read pyvenv.cfg to find the actual Python installation
    pyvenv_cfg <- file.path(venv_path, "pyvenv.cfg")
    if (file.exists(pyvenv_cfg)) {
      cfg_lines <- readLines(pyvenv_cfg)
      home_line <- grep("^home\\s*=", cfg_lines, value = TRUE)

      if (length(home_line) > 0) {
        # Extract the home directory path
        python_home <- trimws(sub("^home\\s*=\\s*", "", home_line[1]))

        # Determine Python version from the venv
        python_bin <- file.path(venv_path, "bin", "python")
        if (file.exists(python_bin)) {
          # Check what library version is expected
          otool_output <- tryCatch({
            system2("otool", args = c("-L", python_bin), stdout = TRUE, stderr = FALSE)
          }, error = function(e) NULL)

          if (!is.null(otool_output)) {
            # Look for libpython reference
            libpython_line <- grep("libpython.*\\.dylib", otool_output, value = TRUE)
            if (length(libpython_line) > 0) {
              # Extract the library name (e.g., libpython3.10.dylib)
              lib_name <- sub(".*/(libpython[0-9.]+\\.dylib).*", "\\1", libpython_line[1])

              # Source: actual Python installation lib directory
              source_lib <- file.path(dirname(python_home), "lib", lib_name)
              # Target: venv lib directory
              target_lib <- file.path(venv_path, "lib", lib_name)

              # Create symlink if source exists and target doesn't
              if (file.exists(source_lib) && !file.exists(target_lib)) {
                file.symlink(source_lib, target_lib)
                message("Created symlink for ", lib_name)
              }
            }
          }
        }
      }
    }
  }

  return(invisible(TRUE))
}

# uv_remove_env() removed - use RemoveEnv(prompt = FALSE) instead

#' Sync Python dependencies using UV
#'
#' @param extras Character vector of extra dependency groups to install
#' @export
uv_sync_deps <- function(extras = "all") {
  pkg_dir <- get_scp_pkg_dir()

  # Change to package directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(pkg_dir)

  # Find pyproject.toml
  pyproject_path <- find_pyproject_toml()
  if (is.null(pyproject_path)) {
    stop("pyproject.toml not found. Cannot sync dependencies.")
  }

  # Copy pyproject.toml to current directory if needed
  if (pyproject_path != "pyproject.toml") {
    file.copy(pyproject_path, "pyproject.toml", overwrite = TRUE)
  }

  # Install dependencies using UV's native pyproject.toml support
  message("Installing Python dependencies from pyproject.toml...")

  # Determine which extras to install
  if (length(extras) == 0 || "" %in% extras || is.null(extras)) {
    # Install only core dependencies (no extras)
    result <- system2("uv", args = c("pip", "install", "-e", "."), wait = TRUE)
  } else if ("all" %in% extras) {
    # Install all extras
    result <- system2("uv", args = c("pip", "install", "-e", ".[all]"), wait = TRUE)
  } else {
    # Install specific extras
    extras_str <- paste0("[", paste(extras, collapse = ","), "]")
    result <- system2("uv", args = c("pip", "install", "-e", paste0(".", extras_str)), wait = TRUE)
  }

  if (result != 0) {
    warning("Installation from pyproject.toml failed. Trying fallback method...")
    # Fallback: install core packages directly by name
    core_packages <- c(
      "numpy", "pandas", "scipy", "matplotlib", "seaborn",
      "scikit-learn", "h5py", "numba", "anndata",
      "scanpy", "igraph"
    )
    result <- system2("uv", args = c("pip", "install", core_packages), wait = TRUE)

    if (result != 0) {
      stop("Failed to install Python dependencies. Check UV installation and pyproject.toml.")
    }
  }

  message("Python dependencies installed successfully.")
  return(invisible(TRUE))
}

#' Install Python dependencies to existing UV environment
#'
#' This function installs Python dependencies to an existing UV environment.
#' You can install feature groups (extras) from pyproject.toml and/or individual packages.
#'
#' @param extras Character vector of extra dependency groups to install from pyproject.toml.
#'   Options include: "velocity", "trajectory", "deeplearning", "singlecell", "ml", "spatial", "all".
#'   Can be NULL to skip installing extras.
#' @param packages Character vector of individual Python package names to install.
#'   Can include version specifications (e.g., "numpy>=1.20.0").
#'   Can be NULL to skip installing individual packages.
#'
#' @return Invisible TRUE if successful
#' @export
#'
#' @examples
#' \dontrun{
#' # Install feature groups only
#' uv_install(extras = "velocity")
#' uv_install(extras = c("velocity", "trajectory"))
#'
#' # Install individual packages only
#' uv_install(packages = "requests")
#' uv_install(packages = c("requests>=2.28.0", "beautifulsoup4"))
#'
#' # Install both extras and packages
#' uv_install(extras = "velocity", packages = "custom_package")
#' }
uv_install <- function(extras = NULL, packages = NULL) {
  # Validate inputs
  if (is.null(extras) && is.null(packages)) {
    stop("At least one of 'extras' or 'packages' must be provided.")
  }

  if (!check_uv()) {
    stop("UV is not installed. Please install UV first with install_uv()")
  }

  if (!uv_env_exists()) {
    stop("UV environment not found. Please run PrepareEnv() first.")
  }

  pkg_dir <- get_scp_pkg_dir()

  # Change to package directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(pkg_dir)

  # Install extras if provided
  if (!is.null(extras)) {
    # Find pyproject.toml
    pyproject_path <- find_pyproject_toml()
    if (is.null(pyproject_path)) {
      stop("pyproject.toml not found. Cannot install extras.")
    }

    # Copy pyproject.toml to current directory if needed
    if (pyproject_path != "pyproject.toml") {
      file.copy(pyproject_path, "pyproject.toml", overwrite = TRUE)
    }

    # Install the extras
    message("Installing Python extras: ", paste(extras, collapse = ", "))

    if (length(extras) == 1 && extras == "all") {
      extras_str <- "[all]"
    } else {
      extras_str <- paste0("[", paste(extras, collapse = ","), "]")
    }

    result <- system2("uv", args = c("pip", "install", "-e", paste0(".", extras_str)), wait = TRUE)

    if (result != 0) {
      stop("Failed to install extras. Check UV installation and pyproject.toml.")
    }

    message("Extras installed successfully.")
  }

  # Install individual packages if provided
  if (!is.null(packages)) {
    message("Installing Python packages: ", paste(packages, collapse = ", "))

    result <- system2("uv", args = c("pip", "install", packages), wait = TRUE)

    if (result != 0) {
      stop("Failed to install packages: ", paste(packages, collapse = ", "))
    }

    message("Packages installed successfully.")
  }

  return(invisible(TRUE))
}

#' Configure reticulate to use UV environment
#'
#' @export
use_uv_env <- function() {
  venv_path <- file.path(get_scp_pkg_dir(), ".venv")
  if (!dir.exists(venv_path)) {
    stop("UV environment not found. Run PrepareEnv() first.")
  }

  # Get Python executable path
  python_path <- get_uv_python_path(venv_path)

  if (!file.exists(python_path)) {
    stop("Python executable not found in UV environment")
  }

  # Suppress deprecation warnings and set environment variables
  # 1. Suppress pkg_resources deprecation warning from reticulate
  Sys.setenv(PYTHONWARNINGS = "ignore::DeprecationWarning:pkg_resources")

  # 2. Fix OpenMP deprecation warning by using new API
  Sys.setenv(OMP_NUM_THREADS = "1")  # Avoid nested parallelism issues
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")  # macOS OpenMP compatibility

  # Fix missing libpython symlink on macOS if needed (for existing environments)
  if (Sys.info()["sysname"] == "Darwin") {
    # Check if Python binary can run
    test_run <- tryCatch({
      system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
    }, error = function(e) NULL, warning = function(w) NULL)

    # If Python binary fails to run (likely missing library), create symlink
    if (is.null(test_run) || length(test_run) == 0) {
      pyvenv_cfg <- file.path(venv_path, "pyvenv.cfg")
      if (file.exists(pyvenv_cfg)) {
        cfg_lines <- readLines(pyvenv_cfg)
        home_line <- grep("^home\\s*=", cfg_lines, value = TRUE)

        if (length(home_line) > 0) {
          python_home <- trimws(sub("^home\\s*=\\s*", "", home_line[1]))

          otool_output <- tryCatch({
            system2("otool", args = c("-L", python_path), stdout = TRUE, stderr = FALSE)
          }, error = function(e) NULL)

          if (!is.null(otool_output)) {
            libpython_line <- grep("libpython.*\\.dylib", otool_output, value = TRUE)
            if (length(libpython_line) > 0) {
              lib_name <- sub(".*/(libpython[0-9.]+\\.dylib).*", "\\1", libpython_line[1])
              source_lib <- file.path(dirname(python_home), "lib", lib_name)
              target_lib <- file.path(venv_path, "lib", lib_name)

              if (file.exists(source_lib) && !file.exists(target_lib)) {
                file.symlink(source_lib, target_lib)
                message("Fixed missing Python library symlink")
              }
            }
          }
        }
      }
    }
  }

  # Set environment variables to fix numba caching issues in venv
  # Numba has trouble with file locators in some venv setups
  Sys.setenv(NUMBA_CACHE_DIR = tempdir())
  Sys.setenv(NUMBA_DISABLE_JIT = "0")  # Keep JIT enabled but use temp cache

  # Use the UV environment
  reticulate::use_python(python_path, required = TRUE)

  # Force Python initialization to ensure it's ready
  tryCatch({
    reticulate::py_config()
  }, error = function(e) {
    stop("Failed to initialize Python: ", e$message,
         "\nTry rebuilding the environment with PrepareEnv(force = TRUE)",
         call. = FALSE)
  })

  invisible(TRUE)
}

#' Install the SCP python environment
#'
#' Install all python packages in the SCP environment automatically using UV.
#' UV is a fast Python package installer that is 10-100x faster than pip/conda.
#'
#' @param force Force reinstall the SCP python environment.
#' @param update Whether to update packages that already exist within the environment.
#' @param python_version Python version to use for the SCP environment. Default is "3.10".
#' @param extras Character vector of extra dependency groups to install (e.g., c("velocity", "trajectory")).
#'   Default is NULL for minimal installation (core packages only).
#'   Use "all" for complete installation, or specify individual groups: "velocity", "trajectory", "deeplearning", "singlecell", "ml".
#'   You can install extras later with uv_install().
#' @export
#'
PrepareEnv <- function(force = FALSE, update = FALSE, python_version = "3.10", extras = NULL) {
  # Check if UV is installed
  if (!check_uv()) {
    stop("UV is not installed. Please install UV first.\n",
         "Installation instructions:\n",
         "  macOS/Linux: curl -LsSf https://astral.sh/uv/install.sh | sh\n",
         "  Windows: powershell -ExecutionPolicy ByPass -c \"irm https://astral.sh/uv/install.ps1 | iex\"\n",
         "Or use: install_uv()")
  }

  # UV-based installation
  if (!uv_env_exists() || isTRUE(force)) {
    if (isTRUE(force) && uv_env_exists()) {
      message("Removing existing UV environment...")
      RemoveEnv(prompt = FALSE)
    }

    message("Creating UV environment with Python ", python_version, "...")
    uv_create_env(python_version = python_version)

    message("Installing Python dependencies with UV...")
    if (is.null(extras) || length(extras) == 0) {
      message("Installing minimal (core) Python packages only.")
      message("To add extras later, use: uv_install(extras = 'velocity') or similar.")
    }
    uv_sync_deps(extras = extras)

    message("UV environment setup complete!")
    message("To use this environment, run: use_uv_env()")
  } else {
    if (isTRUE(update)) {
      message("Updating UV environment...")
      uv_sync_deps(extras = extras)
    } else {
      message("UV environment already exists. Use force=TRUE to reinstall or update=TRUE to update packages.")
      message("To add extras: uv_install(extras = 'velocity') or similar.")
    }
  }
}

# Removed install_py() function - no longer needed with UV-only approach
# Removed try_get() function - unused retry wrapper

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
#' For V5 objects with multiple layers (e.g., split by sample), it automatically joins
#' layers to ensure all cells are included, unless join_layers is set to FALSE.
#'
#' @param srt A Seurat object
#' @param layer The layer to extract (e.g., "counts", "data", "scale.data")
#' @param assay Name of the assay. If NULL, uses the default assay.
#' @param join_layers Logical. For V5 objects with multiple layers, should layers be joined?
#'   Default is TRUE to ensure all cells are included. Set to FALSE for per-sample operations.
#'
#' @return A matrix of the requested data
#' @importFrom SeuratObject JoinLayers Layers
#' @export
get_seurat_data <- function(srt, layer = "data", assay = NULL, join_layers = TRUE) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  assay <- assay %||% DefaultAssay(srt)

  # Check SeuratObject package version to determine parameter name
  sobj_version <- tryCatch(packageVersion("SeuratObject"), error = function(e) NULL)
  use_layer_param <- !is.null(sobj_version) && sobj_version >= "5.0.0"

  # For V5 objects with multiple layers, join them if requested
  is_v5 <- IsSeurat5(srt)
  if (is_v5 && isTRUE(join_layers)) {
    assay_obj <- srt[[assay]]

    # Search for layers matching the requested layer name
    matching_layers <- tryCatch({
      SeuratObject::Layers(assay_obj, search = layer)
    }, error = function(e) NULL)

    # If multiple layers found, join them
    if (!is.null(matching_layers) && length(matching_layers) > 1) {
      message("Found ", length(matching_layers), " layers for '", layer, "': ",
              paste(matching_layers, collapse = ", "))
      message("Joining layers to include all cells. Set join_layers=FALSE to disable.")

      # Join only the specific layers we need
      srt[[assay]] <- tryCatch({
        SeuratObject::JoinLayers(assay_obj, layers = matching_layers)
      }, error = function(e) {
        warning("Failed to join layers: ", e$message,
                ". Returning first layer only.", call. = FALSE)
        assay_obj
      })
    }
  }

  # Extract data using appropriate parameter name based on package version
  if (use_layer_param) {
    # SeuratObject >= 5.0.0: use 'layer' parameter
    return(GetAssayData(srt, layer = layer, assay = assay))
  } else {
    # SeuratObject < 5.0.0: use 'slot' parameter
    return(GetAssayData(srt, slot = layer, assay = assay))
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

  assay <- assay %||% DefaultAssay(srt)

  # Check SeuratObject package version to determine parameter name
  sobj_version <- tryCatch(packageVersion("SeuratObject"), error = function(e) NULL)
  use_layer_param <- !is.null(sobj_version) && sobj_version >= "5.0.0"

  # Set data using appropriate parameter name based on package version
  if (use_layer_param) {
    # SeuratObject >= 5.0.0: use 'layer' parameter
    srt <- SetAssayData(srt, layer = layer, new.data = data, assay = assay)
  } else {
    # SeuratObject < 5.0.0: use 'slot' parameter
    srt <- SetAssayData(srt, slot = layer, new.data = data, assay = assay)
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
    obj_version <- attr(srt, "version")
    if (!is.null(obj_version)) {
      major_version <- extract_major_version(obj_version)
      if (!is.na(major_version) && major_version >= 5) {
        return(TRUE)
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
    pkg_version <- tryCatch(packageVersion("Seurat"), error = function(e) NULL)
    if (!is.null(pkg_version)) {
      major_version <- extract_major_version(pkg_version)
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
#' This function removes the UV-based SCP Python environment (.venv) from your system.
#' It will prompt for confirmation before deletion unless `prompt = FALSE`.
#'
#' @param prompt Whether to prompt for confirmation before removing. Default is TRUE.
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
RemoveEnv <- function(prompt = TRUE) {
  # Check if UV environment exists
  if (!uv_env_exists()) {
    message("UV environment (.venv) not found.")
    return(invisible(NULL))
  }

  venv_path <- file.path(get_scp_pkg_dir(), ".venv")

  # Get environment information
  env_size <- dir_size(venv_path)

  # Find Python executable
  python_path <- get_uv_python_path(venv_path)

  python_version <- NA
  if (file.exists(python_path)) {
    python_version <- tryCatch({
      system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
    }, error = function(e) NA, warning = function(w) NA)
  }

  # Display environment information and prompt for confirmation
  message("About to remove SCP UV Python environment:")
  message("- Environment path: ", venv_path)
  # Safely check python_version - handle vectors and NAs
  if (length(python_version) > 0 && !is.na(python_version[1])) {
    message("- Python version: ", python_version[1])
  }
  message("- Disk space to be freed: ", format_size(env_size))

  proceed <- TRUE
  if (prompt && interactive()) {
    proceed <- utils::menu(c("Yes", "No"),
      title = "Are you sure you want to remove this environment?"
    ) == 1
  }

  if (proceed) {
    message("Removing UV environment...")

    success <- tryCatch({
      unlink(venv_path, recursive = TRUE)
      !file.exists(venv_path)
    }, error = function(e) FALSE)

    if (success) {
      message("Successfully removed UV environment")
      # Unset Python environment variables to avoid pointing to deleted environment
      Sys.unsetenv("RETICULATE_PYTHON")
      Sys.unsetenv("RETICULATE_PYTHON_ENV")

      # Clear reticulate's python configuration cache
      if (exists(".globals", envir = asNamespace("reticulate"))) {
        .globals <- get(".globals", envir = asNamespace("reticulate"))
        .globals$py_config <- NULL
      }
    } else {
      message("Failed to remove environment")
      message("You may need to manually delete: ", venv_path)
    }
  } else {
    message("Environment removal cancelled.")
  }

  return(invisible(NULL))
}

# Removed ListEnv(), VerifyEnv(), and TestPythonCompatibility() functions
# These diagnostic functions added unnecessary complexity without providing
# essential functionality. Users can check environment with uv_env_exists()
# and use_uv_env() directly.
