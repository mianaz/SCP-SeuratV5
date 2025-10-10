
#' Check and report the type of data
#'
#' This function checks the type of data and returns a string indicating the type of data. 
#' It checks for the presence of infinite values, negative values, and whether the values are floats or integers.
#'
#' @param srt An object of class 'Seurat'.
#' @param data The input data. If not provided, it will be extracted from the the 'srt' object.
#' @param layer The layer in the 'srt' object from which to extract the data. Default is "data".
#' @param assay The assay to extract the data from. If not provided, the default assay will be used.
#'
#' @return A string indicating the type of data. Possible values are: "raw_counts", "log_normalized_counts", "raw_normalized_counts", or "unknown".
#'
#' @importFrom Seurat DefaultAssay GetAssayData

#' @importFrom SeuratObject LayerData
#' @export
check_DataType <- function(srt, data = NULL, layer = "data", assay = NULL) {
  if (!is.null(srt)) {
    assay <- assay %||% DefaultAssay(srt)
    
    if (is.null(data)) {
      data <- get_seurat_data(srt, layer = layer, assay = assay)
    }
  }
  
  if (is.null(data)) {
    stop("No data provided and unable to extract data from Seurat object")
  }
  
  # Handle empty layers (common in V5 before normalization)
  if (length(data) == 0 || (is.matrix(data) && all(dim(data) == 0))) {
    # For V5 objects, if data layer is empty but counts exist, assume raw_counts
    if (!is.null(srt) && layer == "data") {
      counts_data <- tryCatch({
        get_seurat_data(srt, layer = "counts", assay = assay)
      }, error = function(e) NULL)
      
      if (!is.null(counts_data) && length(counts_data) > 0) {
        return("raw_counts")  # Data layer empty but counts exist = needs normalization
      }
    }
    return("unknown")
  }
  
  if (inherits(data, "dgCMatrix")) {
    # Handle sparse matrix case
    if (length(data@x) == 0) {
      return("unknown")  # Empty sparse matrix
    }
    isfinite <- all(is.finite(range(data@x, na.rm = TRUE)))
    isfloat <- any(data@x %% 1 != 0, na.rm = TRUE)
    islog <- is.finite(expm1(x = max(data@x, na.rm = TRUE)))
  } else {
    # Handle dense matrix case
    isfinite <- all(is.finite(range(data, na.rm = TRUE)))
    isfloat <- any(data[, sample(seq_len(ncol(data)), min(ncol(data), 1000))] %% 1 != 0, na.rm = TRUE)
    islog <- is.finite(expm1(x = max(data, na.rm = TRUE)))
  }
  
  isnegative <- any(data < 0)

  if (!isTRUE(isfinite)) {
    warning("Infinite values detected!", immediate. = TRUE)
    return("unknown")
  } else if (isTRUE(isnegative)) {
    warning("Negative values detected!", immediate. = TRUE)
    return("unknown")
  } else {
    if (!isfloat) {
      return("raw_counts")
    } else if (isfloat && islog) {
      return("log_normalized_counts")
    } else if (isfloat && !islog) {
      if (isFALSE(isnegative)) {
        return("raw_normalized_counts")
      } else {
        return("unknown")
      }
    }
  }
}


#' Check and preprocess a list of seurat objects
#'
#' This function checks and preprocesses a list of seurat objects. It performs various checks on the input, including verification of input types, assay type consistency, feature name consistency, and batch column consistency. It also performs data normalization and variable feature finding based on the specified parameters. Finally, it prepares the data for integration analysis based on the highly variable features.
#'
#' @param srtList A list of Seurat objects to be checked and preprocessed.
#' @param batch A character string specifying the batch variable name.
#' @param assay The name of the assay to be used for downstream analysis.
#' @param do_normalization A logical value indicating whether data normalization should be performed.
#' @param normalization_method The normalization method to be used. Possible values are "LogNormalize", "SCT", and "TFIDF". Default is "LogNormalize".
#' @param do_HVF_finding A logical value indicating whether highly variable feature (HVF) finding should be performed. Default is TRUE.
#' @param HVF_source The source of highly variable features. Possible values are "global" and "separate". Default is "separate".
#' @param HVF_method The method for selecting highly variable features. Default is "vst".
#' @param nHVF The number of highly variable features to select. Default is 2000.
#' @param HVF_min_intersection The feature needs to be present in batches for a minimum number of times in order to be considered as highly variable. The default value is 1.
#' @param HVF A vector of highly variable features. Default is NULL.
#' @param vars_to_regress A vector of variable names to include as additional regression variables. Default is NULL.
#' @param seed An integer specifying the random seed for reproducibility. Default is 11.
#'
#' @return A list containing the preprocessed seurat objects, the highly variable features, the assay name, and the type of assay (e.g., "RNA" or "Chromatin").
#'
#' @importFrom Seurat SplitObject GetAssayData Assays NormalizeData FindVariableFeatures SCTransform SCTResults SelectIntegrationFeatures PrepSCTIntegration DefaultAssay DefaultAssay<- VariableFeatures VariableFeatures<- UpdateSeuratObject
#' @importFrom SeuratObject LayerData
#' @importFrom Signac RunTFIDF
#' @importFrom Matrix rowSums
#' @importFrom utils head packageVersion
#' @export
#'
check_srtList <- function(srtList, batch, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                          nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          vars_to_regress = NULL, seed = 11, check_v5 = TRUE) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  set.seed(seed)

  if (!inherits(srtList, "list") || any(sapply(srtList, function(x) !inherits(x, "Seurat")))) {
    stop("'srtList' is not a list of Seurat object.")
  }
  
  # Check Seurat version compatibility - V5 safe version
  if (check_v5) {
    tryCatch({
      seurat_version <- packageVersion("Seurat")
      # Use string parsing to avoid numeric conversion issues
      if (as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1]) >= 5) {
        # Check if any objects are V4 and need conversion
        v4_objects <- which(!sapply(srtList, IsSeurat5))
        
        if (length(v4_objects) > 0) {
          message("Found ", length(v4_objects), " Seurat v4 objects. Updating to v5...")
          for (i in v4_objects) {
            message("Converting object ", i, " of ", length(v4_objects), " to Seurat V5...")
            srtList[[i]] <- UpdateSeuratObject(srtList[[i]])
          }
          message("All objects updated to Seurat V5")
        }
      } else {
        warning("You are using Seurat v", seurat_version, ". SCP now requires Seurat v5 for full compatibility.",
                " Some functions may not work as expected.")
      }
    }, error = function(e) {
      warning("Failed to check Seurat version compatibility: ", e$message)
    })
  }
  if (!normalization_method %in% c("LogNormalize", "SCT", "TFIDF")) {
    stop("'normalization_method' must be one of: 'LogNormalize', 'SCT', 'TFIDF'")
  }
  if (normalization_method %in% c("SCT")) {
    check_R("glmGamPoi")
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global', 'separate'")
  }
  if (any(sapply(srtList, ncol) < 2)) {
    stop(paste0(
      "Seurat objects in srtList contain less than 2 cells. srtList index: ",
      paste0(which(sapply(srtList, ncol) < 2), collapse = ",")
    ))
  }

  if (is.null(assay)) {
    default_assay <- unique(sapply(srtList, DefaultAssay))
    if (length(default_assay) != 1) {
      stop("The default assay name of the Seurat object in the srtlist is inconsistent.")
    } else {
      assay <- default_assay
    }
  }

  assay_type <- unique(sapply(srtList, function(srt) class(srt[[assay]])))
  if (length(assay_type) != 1) {
    stop("The assay type of the Seurat object in the srtlist is inconsistent.")
  } else {
    if (assay_type == "Assay") {
      type <- "RNA"
    } else if (assay_type == "ChromatinAssay") {
      type <- "Chromatin"
    } else {
      type <- "Unknown"
    }
  }

  features_list <- lapply(srtList, function(srt) {
    sort(rownames(srt[[assay]]))
  })
  if (length(unique(features_list)) != 1) {
    if (type == "Chromatin") {
      warning("The peaks in assay ", assay, " is different between batches. Creating a common set of peaks and may take a long time...")
      srtMerge <- Reduce(merge, srtList)
      srtList <- SplitObject(object = srtMerge, split.by = batch)
    }
    cf <- Reduce(intersect, lapply(srtList, function(srt) rownames(srt[[assay]])))
    warning("'srtList' have different feature names! Will subset the common features(", length(cf), ") for downstream analysis!", immediate. = TRUE)
    for (i in seq_along(srtList)) {
      # V5-compatible feature subsetting
      if (IsSeurat5(srtList[[i]])) {
        # For V5 objects, subset the entire Seurat object instead of just the assay
        srtList[[i]] <- subset(srtList[[i]], features = cf)
      } else {
        # For V4 objects, use original approach
        srtList[[i]][[assay]] <- subset(srtList[[i]][[assay]], features = cf)
      }
    }
  }

  celllist <- unlist(lapply(srtList, colnames))
  if (length(celllist) != length(unique(celllist))) {
    stop("'srtList' have duplicated cell names!")
  }

  if (length(batch) != 1 && length(batch) != length(srtList)) {
    stop("'batch' must be a character to specify the batch column in the meta.data or a vector of the same length of the srtList!")
  }
  if (length(batch) == length(srtList)) {
    srtList_tmp <- list()
    for (bat in unique(batch)) {
      srtList_tmp[[bat]] <- Reduce(merge, srtList[batch == bat])
    }
    srtList <- srtList_tmp
  } else {
    # Handle empty batch string - skip batch processing
    if (batch != "" && !all(sapply(srtList, function(x) {
      batch %in% colnames(x@meta.data)
    }))) {
      stop(paste0("batch column('", batch, "') was not found in one or more object of the srtList!"))
    }
    if (batch != "") {
      for (i in seq_along(srtList)) {
        u <- unique(srtList[[i]][[batch, drop = TRUE]])
        if (length(u) > 1) {
          x <- SplitObject(srtList[[i]], split.by = batch)
          srtList[[i]] <- character(0)
          srtList <- c(srtList, x)
        }
      }
    }
    srtList <- srtList[sapply(srtList, length) > 0]
    if (batch != "") {
      srtList_batch <- sapply(srtList, function(x) unique(x[[batch, drop = TRUE]]))
      batch_to_merge <- names(which(table(srtList_batch) > 1))
      if (length(batch_to_merge) > 0) {
        for (b in batch_to_merge) {
          index <- which(srtList_batch == b)
          srtList_tmp <- Reduce(merge, srtList[index])
          for (i in index) {
            srtList[[i]] <- character(0)
          }
          srtList <- c(srtList, srtList_tmp)
        }
      }
      srtList <- srtList[sapply(srtList, length) > 0]
    }
  }

  for (i in seq_along(srtList)) {
    if (!assay %in% Assays(srtList[[i]])) {
      stop(paste0("srtList ", i, " does not contain '", assay, "' assay."))
    }
    DefaultAssay(srtList[[i]]) <- assay
    if (isTRUE(do_normalization)) {
      if (normalization_method == "LogNormalize") {
        cat("Perform NormalizeData(LogNormalize) on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- NormalizeData(object = srtList[[i]], assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
      }
      if (normalization_method == "TFIDF") {
        cat("Perform RunTFIDF on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- RunTFIDF(object = srtList[[i]], assay = assay, verbose = FALSE)
      }
    } else if (is.null(do_normalization)) {
      # V5-safe check_DataType call
      tryCatch({
        status <- check_DataType(srtList[[i]], layer = "data", assay = assay)
        if (status == "log_normalized_counts") {
          cat("Data ", i, "/", length(srtList), " of the srtList has been log-normalized.\n", sep = "")
        }
        if (status %in% c("raw_counts", "raw_normalized_counts")) {
          if (normalization_method == "LogNormalize") {
            cat("Data ", i, "/", length(srtList), " of the srtList is ", status, ". Perform NormalizeData(LogNormalize) on the data ...\n", sep = "")
            srtList[[i]] <- NormalizeData(object = srtList[[i]], assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
          }
          if (normalization_method == "TFIDF") {
            cat("Data ", i, "/", length(srtList), " of the srtList is ", status, ". Perform RunTFIDF on the data ...\n", sep = "")
            srtList[[i]] <- RunTFIDF(object = srtList[[i]], assay = assay, verbose = FALSE)
          }
        }
        if (status == "unknown") {
          warning("Can not determine whether data ", i, " is log-normalized...", immediate. = TRUE)
        }
      }, error = function(e) {
        warning("Failed to check data type for object ", i, ": ", e$message, ". Assuming normalization is needed.")
        # Default to normalization if we can't check the data type
        if (normalization_method == "LogNormalize") {
          cat("Defaulting to NormalizeData(LogNormalize) for data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
          srtList[[i]] <- NormalizeData(object = srtList[[i]], assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
        }
      })
    }
    if (is.null(HVF)) {
      if (isTRUE(do_HVF_finding) || is.null(do_HVF_finding) || length(VariableFeatures(srtList[[i]], assay = assay)) == 0) {
        # if (type == "RNA") {
        cat("Perform FindVariableFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- FindVariableFeatures(srtList[[i]], assay = assay, nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE)
        # }
        # if (type == "Chromatin") {
        #   cat("Perform FindTopFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        #   srtList[[i]] <- FindTopFeatures(srtList[[i]], assay = assay, min.cutoff = HVF_min_cutoff, verbose = FALSE)
        # }
      }
    }

    if (normalization_method %in% c("SCT") && type == "RNA") {
      if (isTRUE(do_normalization) || isTRUE(do_HVF_finding) || !"SCT" %in% Assays(srtList[[i]])) {
        cat("Perform SCTransform on the data", i, "of the srtList...\n")
        srtList[[i]] <- SCTransform(
          object = srtList[[i]],
          variable.features.n = nHVF,
          vars.to.regress = vars_to_regress,
          assay = assay,
          method = "glmGamPoi",
          new.assay.name = "SCT",
          verbose = FALSE
        )
      } else {
        DefaultAssay(srtList[[i]]) <- "SCT"
      }
      if (!"residual_variance" %in% colnames(srtList[[i]]@assays$SCT@meta.features)) {
        if (length(srtList[[i]]@assays$SCT@SCTModel.list) > 1) {
          index <- which(sapply(srtList[[i]]@assays$SCT@SCTModel.list, function(x) nrow(x@cell.attributes) == ncol(srtList[[i]])))
        } else {
          index <- 1
        }
        model <- srtList[[i]]@assays$SCT@SCTModel.list[[index]]
        feature.attr <- SCTResults(object = model, layer = "feature.attributes")
      } else {
        feature.attr <- srtList[[i]]@assays$SCT@meta.features
      }
      nfeatures <- min(nHVF, nrow(x = feature.attr))
      top.features <- rownames(x = feature.attr)[head(order(feature.attr$residual_variance, decreasing = TRUE), n = nfeatures)]
      VariableFeatures(srtList[[i]], assay = DefaultAssay(srtList[[i]])) <- top.features
      srtList[[i]]@assays$SCT@meta.features <- feature.attr
    }
  }

  if (is.null(HVF)) {
    if (HVF_source == "global") {
      cat("Use the global HVF from merged dataset...\n")
      srtMerge <- Reduce(merge, srtList)
      # if (type == "RNA") {
      srtMerge <- FindVariableFeatures(srtMerge, assay = DefaultAssay(srtMerge), nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE)
      # }
      # if (type == "Chromatin") {
      #   srtMerge <- FindTopFeatures(srtMerge, assay = DefaultAssay(srtMerge), min.cutoff = HVF_min_cutoff, verbose = FALSE)
      # }
      HVF <- VariableFeatures(srtMerge)
    }
    if (HVF_source == "separate") {
      cat("Use the separate HVF from srtList...\n")
      # if (type == "RNA") {
      HVF <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
      HVF_sort <- sort(table(unlist(lapply(srtList, VariableFeatures))), decreasing = TRUE)
      HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      HVF <- intersect(HVF, names(HVF_filter))
      # }
      # if (type == "Chromatin") {
      #   nHVF <- min(sapply(srtList, function(srt) length(VariableFeatures(srt))))
      #   HVF_sort <- sort(table(unlist(lapply(srtList, VariableFeatures))), decreasing = TRUE)
      #   HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      #   HVF <- names(head(HVF_filter, nHVF))
      # }
      if (length(HVF) == 0) {
        stop("No HVF available.")
      }
    }
  } else {
    cf <- Reduce(intersect, lapply(srtList, function(srt) {
      rownames(get_seurat_data(srt, layer = "counts", assay = DefaultAssay(srt)))
    }))
    HVF <- HVF[HVF %in% cf]
  }
  message("Number of available HVF: ", length(HVF))
  
  # Use the utility function get_seurat_data for data access

  hvf_sum <- lapply(srtList, function(srt) {
    count_data <- get_seurat_data(srt, layer = "counts")
    colSums(count_data[HVF, , drop = FALSE])
  })
  cell_all <- unlist(unname(hvf_sum))
  cell_abnormal <- names(cell_all)[cell_all == 0]
  if (length(cell_abnormal) > 0) {
    warning("Some cells do not express any of the highly variable features: ", paste(cell_abnormal, collapse = ","), immediate. = TRUE)
  }

  if (normalization_method == "SCT" && type == "RNA") {
    srtList <- PrepSCTIntegration(object.list = srtList, anchor.features = HVF, assay = "SCT", verbose = FALSE)
  }
  cat(paste0("[", Sys.time(), "]", " Finished checking.\n"))

  return(list(
    srtList = srtList,
    HVF = HVF,
    assay = assay,
    type = type
  ))
}

#' Check and preprocess a merged seurat object
#'
#' This function checks and preprocesses a merged seurat object.
#' @seealso \code{\link{check_srtList}}
#'
#' @inheritParams check_srtList
#' @param srtMerge A merged Seurat object that includes the batch information.
#'
#' @inheritParams Integration_SCP
#' @importFrom Seurat GetAssayData SplitObject SetAssayData VariableFeatures VariableFeatures<- UpdateSeuratObject
#' @importFrom SeuratObject LayerData
#' @export
check_srtMerge <- function(srtMerge, batch = NULL, assay = NULL,
                           do_normalization = NULL, normalization_method = "LogNormalize",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                           nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                           vars_to_regress = NULL, seed = 11, check_v5 = TRUE) {
  if (!inherits(srtMerge, "Seurat")) {
    stop("'srtMerge' is not a Seurat object.")
  }
  
  # Check Seurat version compatibility
  if (check_v5) {
    tryCatch({
      seurat_version <- packageVersion("Seurat")
      is_v5_installed <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1]) >= 5
      
      if (is_v5_installed) {
        # Check if object is V4 and needs conversion
        if (!IsSeurat5(srtMerge)) {
          message("Seurat object is V4. Updating to V5...")
          srtMerge <- UpdateSeuratObject(srtMerge)
          message("Object updated to Seurat V5")
        }
      } else {
        warning("You are using Seurat v", seurat_version, ". SCP now requires Seurat v5 for full compatibility.",
                " Some functions may not work as expected.")
      }
    }, error = function(e) {
      warning("Failed to check Seurat version compatibility: ", e$message)
    })
  }
  if (length(batch) != 1) {
    stop("'batch' must be provided to specify the batch column in the meta.data")
  }
  if (!batch %in% colnames(srtMerge[[]])) {
    stop(paste0("No batch column('", batch, "') found in the meta.data"))
  }
  if (!is.factor(srtMerge[[batch, drop = TRUE]])) {
    srtMerge[[batch, drop = TRUE]] <- factor(srtMerge[[batch, drop = TRUE]],
      levels = unique(srtMerge[[batch, drop = TRUE]])
    )
  }
  assay <- assay %||% DefaultAssay(srtMerge)
  srtMerge_raw <- srtMerge

  cat(paste0("[", Sys.time(), "]", " Spliting srtMerge into srtList by column ", batch, "... ...\n"))
  srtList <- SplitObject(object = srtMerge_raw, split.by = batch)

  checked <- check_srtList(
    srtList = srtList, batch = batch, assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srtList <- checked[["srtList"]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  srtMerge <- Reduce(merge, srtList)

  srtMerge <- SrtAppend(
    srt_raw = srtMerge, srt_append = srtMerge_raw, pattern = "",
    layers = "reductions", overwrite = TRUE, verbose = FALSE
  )
  if (normalization_method == "SCT" && type == "RNA") {
    DefaultAssay(srtMerge) <- "SCT"
  } else {
    DefaultAssay(srtMerge) <- assay
  }
  VariableFeatures(srtMerge) <- HVF

  return(list(
    srtMerge = srtMerge,
    srtList = srtList,
    HVF = HVF,
    assay = assay,
    type = type
  ))
}

#' Attempt to recover raw counts from the normalized matrix.
#'
#' @param srt A Seurat object.
#' @param assay Name of assay to recover counts.
#' @param trans The transformation function to applied when data is presumed to be log-normalized.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#'  However, due to decimal point preservation during normalization, the calculated nCount is usually a floating point number close to the integer.
#'  The tolerance is its difference from the integer. Default is 0.1
#' @param sf Set the scaling factor manually.
#' @param verbose Whether to show messages.
#'
#' @examples
#' data("pancreas_sub")
#' raw_counts <- pancreas_sub@assays$RNA$counts
#'
#' # Normalized the data
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#'
#' # Now replace counts with the log-normalized data matrix
#' pancreas_sub@assays$RNA$counts <- pancreas_sub@assays$RNA$data
#'
#' # Recover the counts and compare with the raw counts matrix
#' pancreas_sub <- RecoverCounts(pancreas_sub)
#' identical(raw_counts, pancreas_sub@assays$RNA$counts)
#' @importFrom Seurat GetAssayData SetAssayData
#' @importFrom SeuratObject LayerData<-
#' @importFrom SeuratObject LayerData
#' @importFrom SeuratObject as.sparse
#' @export
RecoverCounts <- function(srt, assay = NULL, trans = c("expm1", "exp", "none"), min_count = c(1, 2, 3), tolerance = 0.1, sf = NULL, verbose = TRUE) {
  assay <- assay %||% DefaultAssay(srt)
  counts <- get_seurat_data(srt, layer = "counts", assay = assay)
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as.sparse(counts[1:nrow(counts), , drop = FALSE])
  }
  status <- check_DataType(data = counts)
  if (status == "raw_counts") {
    if (isTRUE(verbose)) {
      message("The data is already raw counts.")
    }
    return(srt)
  }
  if (status == "log_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be log-normalized.")
    }
    trans <- match.arg(trans)
    if (trans %in% c("expm1", "exp")) {
      if (isTRUE(verbose)) {
        message("Perform ", trans, " on the raw data.")
      }
      counts <- do.call(trans, list(counts))
    }
  }
  if (status == "raw_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be normalized without log transformation.")
    }
  }
  if (is.null(sf)) {
    sf <- unique(round(colSums(counts)))
    if (isTRUE(verbose)) {
      message("The presumed scale factor: ", paste0(head(sf, 10), collapse = ", "))
    }
  }
  if (length(sf) == 1) {
    counts <- counts / sf
    elements <- split(counts@x, rep(1:ncol(counts), diff(counts@p)))
    min_norm <- sapply(elements, min)
    nCount <- NULL
    for (m in min_count) {
      if (is.null(nCount)) {
        presumed_nCount <- m / min_norm
        diff_value <- abs(presumed_nCount - round(presumed_nCount))
        if (max(diff_value, na.rm = TRUE) < tolerance) {
          nCount <- round(presumed_nCount)
        }
      }
    }
    if (is.null(nCount)) {
      warning("The presumed nCount of some cells is not valid: ", paste0(head(colnames(counts)[diff_value < tolerance], 10), collapse = ","), ", ...", immediate. = TRUE)
      return(srt)
    }
    counts@x <- round(counts@x * rep(nCount, diff(counts@p)))
    srt <- set_seurat_data(srt, data = counts, layer = "counts", assay = assay)
    srt[[paste0("nCount_", assay)]] <- nCount
  } else {
    warning("Scale factor is not unique. No changes to be made.", immediate. = TRUE)
  }
  return(srt)
}

#' Rename features for the Seurat object
#'
#' @param srt A Seurat object.
#' @param newnames A vector with the same length of features in Seurat object, or characters named with old features.
#' @param assays Assays to rename.
#'
#' @examples
#' data("panc8_sub")
#' head(rownames(panc8_sub))
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
#' panc8_rename <- RenameFeatures(panc8_sub, newnames = genenames)
#' head(rownames(panc8_rename))
#'
#' @importFrom Seurat Assays GetAssay
#' @export
RenameFeatures <- function(srt, newnames = NULL, assays = NULL) {
  assays <- assays[assays %in% Assays(srt)] %||% Assays(srt)
  if (is.null(names(newnames))) {
    if (!identical(length(newnames), nrow(srt))) {
      stop("'newnames' must be named or the length of features in the srt.")
    }
    if (length(unique(sapply(pancreas_sub@assays[assays], nrow))) > 1) {
      stop("Assays in the srt object have different number of features. Please use a named vectors.")
    }
    names(newnames) <- rownames(srt[[assays[1]]])
  }
  
  # Check if using Seurat V5
  is_v5 <- IsSeurat5(srt)
  
  for (assay in assays) {
    message("Rename features for the assay: ", assay)
    assay_obj <- GetAssay(srt, assay)
    
    if (is_v5) {
      # Seurat V5 approach
      if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' is required for Seurat V5 support")
      }
      
      # Handle feature metadata renaming
      meta_features <- SeuratObject::FetchData(assay_obj, vars = NULL, layer = "meta.features")
      if (nrow(meta_features) > 0) {
        index <- which(rownames(meta_features) %in% names(newnames))
        rownames(meta_features)[index] <- newnames[rownames(meta_features)[index]]
        assay_obj <- SeuratObject::AddMetaData(assay_obj, meta_features, layer = "meta.features")
      }
      
      # Handle layer data renaming for different layers
      for (layer_name in c("counts", "data", "scale.data")) {
        if (layer_name %in% SeuratObject::Layers(assay_obj)) {
          layer_data <- SeuratObject::LayerData(assay_obj, layer = layer_name)
          if (nrow(layer_data) > 0) {
            index <- which(rownames(layer_data) %in% names(newnames))
            rownames(layer_data)[index] <- newnames[rownames(layer_data)[index]]
            SeuratObject::LayerData(assay_obj, layer = layer_name) <- layer_data
          }
        }
      }
      
      # Handle variable features renaming
      var_features <- assay_obj[["var.features"]]
      if (length(var_features) > 0) {
        index <- which(var_features %in% names(newnames))
        var_features[index] <- newnames[var_features[index]]
        assay_obj[["var.features"]] <- var_features
      }
    } else {
      # Seurat V4 approach - use direct layer access
      for (d in c("meta.features", "scale.data", "counts", "data")) {
        if (length(dim(layer(assay_obj, d))) > 0) {
          index <- which(rownames(layer(assay_obj, d)) %in% names(newnames))
          rownames(layer(assay_obj, d))[index] <- newnames[rownames(layer(assay_obj, d))[index]]
        }
      }
      if (length(layer(assay_obj, "var.features")) > 0) {
        index <- which(layer(assay_obj, "var.features") %in% names(newnames))
        layer(assay_obj, "var.features")[index] <- newnames[layer(assay_obj, "var.features")[index]]
      }
    }
    
    srt[[assay]] <- assay_obj
  }
  
  return(srt)
}

#' Rename clusters for the Seurat object
#'
#' @param srt A Seurat object.
#' @param group.by The old group used to rename cells.
#' @param nameslist A named list of new cluster value.
#' @param name The name of the new cluster stored in the Seurat object.
#' @param keep_levels If the old group is a factor, keep the order of the levels.
#'
#' @examples
#' data("pancreas_sub")
#' levels(pancreas_sub@meta.data[["SubCellType"]])
#'
#' # Rename all clusters
#' pancreas_sub <- RenameClusters(pancreas_sub, group.by = "SubCellType", nameslist = letters[1:8])
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Rename specified clusters
#' pancreas_sub <- RenameClusters(pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list("a" = "Alpha", "b" = "Beta")
#' )
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Merge and rename clusters
#' pancreas_sub <- RenameClusters(pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list("EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")),
#'   name = "Merged", keep_levels = TRUE
#' )
#' CellDimPlot(pancreas_sub, "Merged")
#'
#' @importFrom stats setNames
#' @export
RenameClusters <- function(srt, group.by, nameslist = list(), name = "newclusters", keep_levels = FALSE) {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (missing(group.by)) {
    stop("group.by must be provided")
  }
  
  # Use [[]] operator instead of direct @meta.data access
  if (!group.by %in% colnames(srt[[]])) {
    stop(paste0(group.by, " is not in the meta.data of srt object."))
  }
  
  if (length(nameslist) > 0 && is.null(names(nameslist))) {
    names(nameslist) <- levels(srt[[group.by, drop=TRUE]])
  }
  
  if (is.list(nameslist) && length(nameslist) > 0) {
    names_assign <- setNames(rep(names(nameslist), sapply(nameslist, length)), nm = unlist(nameslist))
  } else {
    if (is.null(names(nameslist))) {
      if (!is.factor(srt[[group.by, drop=TRUE]])) {
        stop("'nameslist' must be named when srt[[group.by]] is not a factor")
      }
      if (!identical(length(nameslist), length(unique(srt[[group.by, drop=TRUE]])))) {
        stop("'nameslist' must be named or the length of ", length(unique(srt[[group.by, drop=TRUE]])))
      }
      names(nameslist) <- levels(srt[[group.by, drop=TRUE]])
    }
    names_assign <- nameslist
  }
  
  if (all(!names(names_assign) %in% srt[[group.by, drop=TRUE]])) {
    stop("No group name mapped.")
  }
  
  if (is.factor(srt[[group.by, drop=TRUE]])) {
    levels <- levels(srt[[group.by, drop=TRUE]])
  } else {
    levels <- NULL
  }
  
  # Convert to metadata dataframe to make modifications
  metadata <- srt[[]]
  index <- which(as.character(metadata[[group.by]]) %in% names(names_assign))
  metadata[[name]] <- as.character(metadata[[group.by]])
  metadata[[name]][index] <- names_assign[metadata[[name]][index]]
  
  if (!is.null(levels)) {
    levels[levels %in% names(names_assign)] <- names_assign[levels[levels %in% names(names_assign)]]
    if (isFALSE(keep_levels)) {
      levels <- unique(c(names_assign, levels))
    } else {
      levels <- unique(levels)
    }
    metadata[[name]] <- factor(metadata[[name]], levels = levels)
  }
  
  # Add metadata column back to the Seurat object
  srt[[name]] <- metadata[[name]]
  
  return(srt)
}

#' Reorder idents by the gene expression
#'
#' @param srt A Seurat object.
#' @param features Features used to reorder idents.
#' @param reorder_by Reorder groups instead of idents.
#' @param layer Specific layer to get data from.
#' @param assay Specific assay to get data from.
#' @param log Whether log1p transformation needs to be applied. Default is \code{TRUE}.
#' @param distance_metric Metric to compute distance. Default is "euclidean".
#'
#' @importFrom Seurat VariableFeatures DefaultAssay DefaultAssay<- AverageExpression AggregateExpression Idents<-
#' @importFrom SeuratObject as.sparse
#' @importFrom stats hclust reorder as.dendrogram as.dist
#' @importFrom Matrix t colMeans
#' @importFrom proxyC simil dist
#' @export
SrtReorder <- function(srt, features = NULL, reorder_by = NULL, layer = "data", assay = NULL, log = TRUE,
                       distance_metric = "euclidean") {
  if (!inherits(srt, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  
  # Add ident column to metadata
  if (is.null(reorder_by)) {
    srt$ident <- Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  
  if (length(unique(srt[[reorder_by, drop = TRUE]])) == 1) {
    warning("Only one cluster found.", immediate. = TRUE)
    return(srt)
  }
  
  simil_method <- c(
    "cosine", "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman",
    "simple matching", "faith"
  )
  dist_method <- c(
    "euclidean", "chisquared", "kullback", "manhattan", "maximum", "canberra",
    "minkowski", "hamming"
  )
  if (!distance_metric %in% c(simil_method, dist_method, "pearson", "spearman")) {
    stop(distance_metric, " method is invalid.")
  }

  # Use AggregateExpression which should work in both V4 and V5
  # Note: layer parameter is passed via ... for compatibility
  data.avg <- tryCatch({
    AggregateExpression(
      object = srt, 
      features = features, 
      assays = assay, 
      group.by = "ident", 
      verbose = FALSE
    )[[1]][features, , drop = FALSE]
  }, error = function(e) {
    # Fallback: use AverageExpression if AggregateExpression fails
    AverageExpression(
      object = srt, 
      features = features, 
      assays = assay, 
      group.by = "ident", 
      verbose = FALSE
    )[[1]][features, , drop = FALSE]
  })
  if (isTRUE(log)) {
    data.avg <- log1p(data.avg)
  }
  mat <- t(x = data.avg[features, , drop = FALSE])
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as.sparse(mat[1:nrow(mat), , drop = FALSE])
  }

  if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 - simil(as.sparse(mat[1:nrow(mat), , drop = FALSE]), method = distance_metric)
  } else if (distance_metric %in% dist_method) {
    d <- dist(as.sparse(mat[1:nrow(mat), , drop = FALSE]), method = distance_metric)
  }
  data.dist <- as.dist(d)
  hc <- hclust(d = data.dist)
  dd <- as.dendrogram(hc)
  dd_ordered <- reorder(dd, wts = colMeans(data.avg[features, , drop = FALSE]), agglo.FUN = mean)
  ident_new <- unname(setNames(object = seq_along(labels(dd_ordered)), nm = labels(dd_ordered))[as.character(srt$ident)])
  ident_new <- factor(ident_new, levels = seq_along(labels(dd_ordered)))
  Idents(srt) <- srt$ident <- ident_new
  return(srt)
}

#' Append a Seurat object to another
#'
#' @param srt_raw A Seurat object to be appended.
#' @param srt_append New Seurat object to append.
#' @param layers layers names.
#' @param pattern A character string containing a regular expression. All data with matching names will be considered for appending.
#' @param overwrite Whether to overwrite.
#' @param verbose Show messages.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom methods slot hasMethod
#' @export
SrtAppend <- function(srt_raw, srt_append,
                      layers = layerNames(srt_append), pattern = NULL, overwrite = FALSE,
                      verbose = TRUE) {
  if (!inherits(srt_raw, "Seurat") || !inherits(srt_append, "Seurat")) {
    stop("'srt_raw' or 'srt_append' is not a Seurat object.")
  }
  
  # Check if objects are Seurat V5
  is_v5_raw <- IsSeurat5(srt_raw)
  is_v5_append <- IsSeurat5(srt_append)
  
  if (is_v5_raw != is_v5_append) {
    warning("Seurat objects have different versions (V4 vs V5). Consider updating both to V5 with UpdateSeuratObject() for best compatibility.")
  }
  
  pattern <- pattern %||% ""
  for (layer_nm in layerNames(srt_append)) {
    if (!layer_nm %in% layers) {
      if (isTRUE(verbose)) {
        message("layer ", layer_nm, " is not appended.")
      }
      next
    }
    
    if (identical(layer_nm, "active.ident") && isTRUE(overwrite)) {
      layer(srt_raw, name = "active.ident") <- layer(srt_append, name = "active.ident")
      next
    }
    
    for (info in names(layer(srt_append, name = layer_nm))) {
      if (is.null(info)) {
        if (length(layer(srt_append, name = layer_nm)) > 0 && isTRUE(overwrite)) {
          layer(srt_raw, name = layer_nm) <- layer(srt_append, name = layer_nm)
        }
        next
      }
      
      if (!grepl(pattern = pattern, x = info)) {
        if (isTRUE(verbose)) {
          message(info, " in layer ", layer_nm, " is not appended.")
        }
        next
      }
      
      if (!info %in% names(layer(srt_raw, name = layer_nm)) || isTRUE(overwrite)) {
        if (layer_nm %in% c("assays", "graphs", "neighbors", "reductions", "images")) {
          if (identical(layer_nm, "graphs")) {
            # Use bracket accessor instead of direct slot access
            if ("graphs" %in% layerNames(srt_raw)) {
              slot_or_layer <- layer(srt_raw, name = "graphs")
              slot_or_layer[[info]] <- srt_append[[info]]
              layer(srt_raw, name = "graphs") <- slot_or_layer
            } else {
              # For V4 compatibility
              if (!is_v5_raw) {
                srt_raw[["graphs"]][[info]] <- srt_append[[info]]
              }
            }
          } else if (identical(layer_nm, "assays")) {
            if (!info %in% Assays(srt_raw)) {
              # Add the entire assay
              srt_raw[[info]] <- srt_append[[info]]
            } else {
              # Update assay components
              # For V5, use LayerData
              if (is_v5_raw) {
                # For V5 we use layers
                if ("counts" %in% Layers(srt_append[[info]])) {
                  layer(srt_raw[[info]], "counts") <- layer(srt_append[[info]], "counts")
                }
                if ("data" %in% Layers(srt_append[[info]])) {
                  layer(srt_raw[[info]], "data") <- layer(srt_append[[info]], "data")
                }
                # Copy other assay attributes
                srt_raw[[info]] <- SetAssayData(srt_raw[[info]], 
                                                slot = "meta.features", 
                                                new.data = cbind(GetAssayData(srt_raw, slot = "meta.features", assay = info),
                                                                GetAssayData(srt_append, slot = "meta.features", assay = info)[
                                                                  rownames(GetAssayData(srt_raw, slot = "meta.features", assay = info)),
                                                                  setdiff(colnames(GetAssayData(srt_append, slot = "meta.features", assay = info)), 
                                                                          colnames(GetAssayData(srt_raw, slot = "meta.features", assay = info)))
                                                                ]))
                # Copy variable features                                                
                VariableFeatures(srt_raw[[info]]) <- VariableFeatures(srt_append[[info]])
              } else {
                # For V4 we can use direct slot access for backward compatibility
                srt_raw[[info]]$counts <- srt_append[[info]]$counts
                srt_raw[[info]]$data <- srt_append[[info]]$data
                VariableFeatures(srt_raw[[info]]) <- VariableFeatures(srt_append[[info]])
                # Update meta.features
                srt_raw[[info]] <- SetAssayData(srt_raw[[info]], 
                                                slot = "meta.features", 
                                                new.data = cbind(GetAssayData(srt_raw, slot = "meta.features", assay = info),
                                                               GetAssayData(srt_append, slot = "meta.features", assay = info)[
                                                                 rownames(GetAssayData(srt_raw, slot = "meta.features", assay = info)),
                                                                 setdiff(colnames(GetAssayData(srt_append, slot = "meta.features", assay = info)), 
                                                                         colnames(GetAssayData(srt_raw, slot = "meta.features", assay = info)))
                                                               ]))
              }
            }
          } else {
            # For other data types, use standard accessor
            srt_raw[[info]] <- srt_append[[info]]
          }
        } else if (identical(layer_nm, "meta.data")) {
          # Update metadata using the standard accessor
          # First remove any existing column with the same name
          if (info %in% colnames(srt_raw[[]])) {
            srt_raw[[info]] <- NULL
          }
          # Add the column from srt_append
          srt_raw[[info]] <- srt_append[[info]][colnames(srt_raw)]
        } else {
          # For other layer types, use the layer accessor
          layer(srt_raw, name = layer_nm)[[info]] <- layer(srt_append, name = layer_nm)[[info]]
        }
      }
    }
  }
  return(srt_raw)
}

#' Run dimensionality reduction
#'
#' @param srt A Seurat object.
#' @param prefix The prefix used to name the result.
#' @param features Use features expression data to run linear or nonlinear dimensionality reduction.
#' @param assay Specific assay to get data from.
#' @param layer Specific layer to get data from.
#' @param linear_reduction Method of linear dimensionality reduction. Options are "pca", "ica", "nmf", "mds", "glmpca".
#' @param linear_reduction_dims Total number of dimensions to compute and store for \code{linear_reduction}.
#' @param linear_reduction_params Other parameters passed to the \code{linear_reduction} method.
#' @param force_linear_reduction Whether force to do linear dimensionality reduction.
#' @param nonlinear_reduction Method of nonlinear dimensionality reduction. Options are "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"
#' @param nonlinear_reduction_dims Total number of dimensions to compute and store for \code{nonlinear_reduction}.
#' @param reduction_use Which dimensional reduction to use as input for \code{nonlinear_reduction}.
#' @param reduction_dims Which dimensions to use as input for \code{nonlinear_reduction}, used only if \code{features} is \code{NULL}.
#' @param neighbor_use Name of neighbor to use for the \code{nonlinear_reduction}.
#' @param graph_use Name of graph to use for the \code{nonlinear_reduction}.
#' @param nonlinear_reduction_params  Other parameters passed to the \code{nonlinear_reduction} method.
#' @param force_nonlinear_reduction Whether force to do nonlinear dimensionality reduction.
#' @param verbose Show messages.
#' @param seed Set a seed.
#'
#' @importFrom Seurat Embeddings RunPCA RunICA RunTSNE Reductions DefaultAssay DefaultAssay<- Key Key<-
#' @importFrom Signac RunSVD
#' @export
RunDimReduction <- function(srt, prefix = "", features = NULL, assay = NULL, layer = "data",
                            linear_reduction = NULL, linear_reduction_dims = 50,
                            linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = NULL, nonlinear_reduction_dims = 2,
                            reduction_use = NULL, reduction_dims = NULL,
                            graph_use = NULL, neighbor_use = NULL,
                            nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            verbose = TRUE, seed = 11) {
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    type <- "Chromatin"
  } else {
    type <- "RNA"
  }
  linear_reduction_dims <- min(linear_reduction_dims, nrow(srt[[assay]]) - 1, ncol(srt[[assay]]) - 1, na.rm = TRUE)
  nonlinear_reduction_dims <- min(nonlinear_reduction_dims, nrow(srt[[assay]]) - 1, ncol(srt[[assay]]) - 1, na.rm = TRUE)
  if (!is.null(linear_reduction)) {
    if (any(!linear_reduction %in% c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt))) || length(linear_reduction) > 1) {
      stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
    }
  }
  if (!is.null(nonlinear_reduction)) {
    if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr", Reductions(srt))) || length(nonlinear_reduction) > 1) {
      stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
    }
    if (is.null(features) && is.null(reduction_use) && is.null(neighbor_use) && is.null(graph_use)) {
      stop("'features', 'reduction_use', 'neighbor_use', or 'graph_use' must be provided when running non-linear dimensionality reduction.")
    }
    if (nonlinear_reduction %in% c("fr")) {
      if (!is.null(graph_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Graphs(", graph_use, ") as input")
      } else if (!is.null(neighbor_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Neighbors(", neighbor_use, ") as input")
      } else if (!is.null(features)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Features(length:", length(features), ") as input")
      } else if (!is.null(reduction_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Reduction(", reduction_use, ", dims:", min(reduction_dims), "-", max(reduction_dims), ") as input")
      }
    } else {
      if (!is.null(features)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Features(length:", length(features), ") as input")
      } else if (!is.null(reduction_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Reduction(", reduction_use, ", dims:", min(reduction_dims), "-", max(reduction_dims), ") as input")
      } else if (!is.null(neighbor_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Neighbors(", neighbor_use, ") as input")
      } else if (!is.null(graph_use)) {
        message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Graphs(", graph_use, ") as input")
      }
    }
  }
  if (!is.null(linear_reduction)) {
    if (!isTRUE(force_linear_reduction)) {
      if (linear_reduction %in% Reductions(srt)) {
        if (srt[[linear_reduction]]@assay.used == assay) {
          message("linear_reduction(", linear_reduction, ") is already existed. Skip calculation.")
          reduc <- srt[[linear_reduction]]
          Key(reduc) <- paste0(prefix, linear_reduction, "_")
          srt[[paste0(prefix, linear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
          return(srt)
        } else {
          message("assay.used is ", srt[[linear_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the linear reduction(pca)")
          linear_reduction <- "pca"
        }
      }
    }
    if (is.null(features) || length(features) == 0) {
      message("No features provided. Use variable features.")
      if (length(VariableFeatures(srt, assay = assay)) == 0) {
        srt <- FindVariableFeatures(srt, assay = assay, verbose = FALSE)
      }
      features <- VariableFeatures(srt, assay = assay)
    }
    fun_use <- switch(linear_reduction,
      "pca" = "RunPCA",
      "svd" = "RunSVD",
      "ica" = "RunICA",
      "nmf" = "RunNMF",
      "mds" = "RunMDS",
      "glmpca" = "RunGLMPCA"
    )
    key_use <- switch(linear_reduction,
      "pca" = "PC_",
      "svd" = "LSI_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    components_nm <- switch(linear_reduction,
      "pca" = "npcs",
      "svd" = "n",
      "ica" = "nics",
      "nmf" = "nbes",
      "mds" = "nmds",
      "glmpca" = "L"
    )
    params <- list(
      object = srt, assay = assay, layer = layer,
      features = features, components_nm = linear_reduction_dims,
      reduction.name = paste0(prefix, linear_reduction),
      reduction.key = paste0(prefix, key_use),
      verbose = verbose, seed.use = seed
    )
    # Parameter filtering is now handled in the switch statement below
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(linear_reduction_params)) {
      params[[nm]] <- linear_reduction_params[[nm]]
    }
    
    # Use direct function calls instead of invoke to avoid parameter conflicts
    srt <- switch(linear_reduction,
      "pca" = {
        # For RunPCA, remove layer parameter as it's not supported
        pca_params <- params
        pca_params[["layer"]] <- NULL
        do.call(RunPCA, pca_params)
      },
      "svd" = {
        params[["assay"]] <- assay  # SVD might need assay parameter
        do.call(RunSVD, params)
      },
      "ica" = {
        params[["assay"]] <- assay  # ICA might need assay parameter
        do.call(RunICA, params)
      },
      "nmf" = {
        params[["assay"]] <- assay
        params[["layer"]] <- layer
        do.call(RunNMF, params)
      },
      "mds" = {
        params[["assay"]] <- assay
        params[["layer"]] <- layer
        do.call(RunMDS, params)
      },
      "glmpca" = {
        params[["assay"]] <- assay
        params[["layer"]] <- "counts"  # GLMPCA always uses counts
        do.call(RunGLMPCA, params)
      }
    )

    if (is.null(rownames(srt[[paste0(prefix, linear_reduction)]]))) {
      rownames(srt[[paste0(prefix, linear_reduction)]]@cell.embeddings) <- colnames(srt)
    }
    if (linear_reduction == "pca") {
      pca.out <- srt[[paste0(prefix, linear_reduction)]]
      center <- rowMeans(get_seurat_data(srt, layer = "scale.data", assay = assay)[features, , drop = FALSE])
      model <- list(sdev = pca.out@stdev, rotation = pca.out@feature.loadings, center = center, scale = FALSE, x = pca.out@cell.embeddings)
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["model"]] <- model
    }
    if (linear_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:linear_reduction_dims
    } else {
      dim_est <- tryCatch(expr = {
        min(
          intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(prefix, linear_reduction)), k = 20)[["dim.est"]],
          ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction)))
        )
      }, error = function(e) {
        message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
        return(NA)
      })
      if (!is.na(dim_est)) {
        dims_estimate <- seq_len(max(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 10), ceiling(dim_est)))
      } else {
        dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 30))
      }
    }
    srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["dims_estimate"]] <- dims_estimate
    srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
  } else if (!is.null(nonlinear_reduction)) {
    if (!isTRUE(force_nonlinear_reduction)) {
      if (nonlinear_reduction %in% Reductions(srt)) {
        if (srt[[nonlinear_reduction]]@assay.used == assay) {
          message("nonlinear_reduction(", nonlinear_reduction, ") is already existed. Skip calculation.")
          reduc <- srt[[nonlinear_reduction]]
          Key(reduc) <- paste0(prefix, nonlinear_reduction, "_")
          srt[[paste0(prefix, nonlinear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction)
          return(srt)
        } else {
          message("assay.used is ", srt[[nonlinear_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the nonlinear reduction(umap)")
          nonlinear_reduction <- "umap"
        }
      }
    }
    # if (!is.null(neighbor_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'neighbor_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    # if (!is.null(graph_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'graph_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    fun_use <- switch(nonlinear_reduction,
      "umap" = "RunUMAP2",
      "umap-naive" = "RunUMAP2",
      "tsne" = "RunTSNE",
      "dm" = "RunDM",
      "phate" = "RunPHATE",
      "pacmap" = "RunPaCMAP",
      "trimap" = "RunTriMap",
      "largevis" = "RunLargeVis",
      "fr" = "RunFR"
    )
    components_nm <- switch(nonlinear_reduction,
      "umap" = "n.components",
      "umap-naive" = "n.components",
      "tsne" = "dim.embed",
      "dm" = "ndcs",
      "phate" = "n_components",
      "pacmap" = "n_components",
      "trimap" = "n_components",
      "largevis" = "n_components",
      "fr" = "ndim"
    )
    other_params <- switch(nonlinear_reduction,
      "umap" = list(umap.method = "uwot", return.model = TRUE),
      "umap-naive" = list(umap.method = "naive", return.model = TRUE),
      "tsne" = list(tsne.method = "Rtsne", num_threads = 0, check_duplicates = FALSE),
      "dm" = list(),
      "phate" = list(),
      "pacmap" = list(),
      "trimap" = list(),
      "largevis" = list(),
      "fr" = list()
    )
    nonlinear_reduction_sim <- toupper(gsub(pattern = "-.*", replacement = "", x = nonlinear_reduction))
    params <- list(
      object = srt, assay = assay, layer = layer, components_nm = nonlinear_reduction_dims,
      features = features, reduction = reduction_use, dims = reduction_dims,
      reduction.name = paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D"),
      reduction.key = paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D_"),
      verbose = verbose, seed.use = seed
    )
    if (!is.null(neighbor_use)) {
      params[["neighbor"]] <- neighbor_use
    }
    if (!is.null(graph_use)) {
      params[["graph"]] <- graph_use
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(other_params)) {
      params[[nm]] <- other_params[[nm]]
    }
    for (nm in names(nonlinear_reduction_params)) {
      params[[nm]] <- nonlinear_reduction_params[[nm]]
    }
    
    # Use direct function calls instead of invoke to avoid parameter conflicts
    srt <- switch(nonlinear_reduction,
      "umap" = {
        # For RunUMAP2, ensure proper parameter mapping
        do.call(RunUMAP2, params)
      },
      "umap-naive" = {
        do.call(RunUMAP2, params)
      },
      "tsne" = {
        do.call(RunTSNE, params)
      },
      "dm" = {
        do.call(RunDM, params)
      },
      "phate" = {
        do.call(RunPHATE, params)
      },
      "pacmap" = {
        do.call(RunPaCMAP, params)
      },
      "trimap" = {
        do.call(RunTriMap, params)
      },
      "largevis" = {
        do.call(RunLargeVis, params)
      },
      "fr" = {
        do.call(RunFR, params)
      }
    )

    srt@reductions[[paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D")]]@misc[["reduction_dims"]] <- reduction_dims
    srt@reductions[[paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D")]]@misc[["reduction_use"]] <- reduction_use
    srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction_sim)
  }
  return(srt)
}

#' Find the default reduction name in a Seurat object.
#'
#' @param srt A Seurat object.
#' @param pattern Character string containing a regular expression to search for.
#' @param min_dim Minimum dimension threshold.
#' @param max_distance Maximum distance allowed for a match.
#'
#' @examples
#' data("pancreas_sub")
#' names(pancreas_sub@reductions)
#' DefaultReduction(pancreas_sub)
#'
#' # Searches for matches to "pca"
#' DefaultReduction(pancreas_sub, pattern = "pca")
#'
#' # Searches for approximate matches to "pc"
#' DefaultReduction(pancreas_sub, pattern = "pc")
#'
#' @return Default reduction name.
#'
#' @export
DefaultReduction <- function(srt, pattern = NULL, min_dim = 2, max_distance = 0.1) {
  if (length(srt@reductions) == 0) {
    stop("Unable to find any reductions.")
  }
  pattern_default <- c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr", "pca", "svd", "ica", "nmf", "mds", "glmpca")
  pattern_dim <- c("2D", "3D")
  reduc_all <- names(srt@reductions)
  reduc_all <- reduc_all[unlist(lapply(reduc_all, function(x) {
    dim(srt@reductions[[x]]@cell.embeddings)[2] >= min_dim
  }))]
  if (length(reduc_all) == 0) {
    stop("No dimensional reduction found in the srt object.")
  }
  if (length(reduc_all) == 1) {
    return(reduc_all)
  }
  if (is.null(pattern)) {
    if (("Default_reduction" %in% names(srt@misc))) {
      pattern <- srt@misc[["Default_reduction"]]
    } else {
      pattern <- pattern_default
    }
  }

  pattern <- c(pattern, paste0(pattern, min_dim, "D"))
  if (any(pattern %in% reduc_all)) {
    return(pattern[pattern %in% reduc_all][1])
  }
  index <- c(unlist(sapply(pattern, function(pat) {
    grep(pattern = pat, x = reduc_all, ignore.case = TRUE)
  })))
  if (length(index) > 0) {
    default_reduc <- reduc_all[index]
  } else {
    index <- c(unlist(sapply(pattern, function(pat) {
      agrep(pattern = pat, x = reduc_all, max.distance = max_distance, ignore.case = TRUE)
    })))
    if (length(index) > 0) {
      default_reduc <- reduc_all[index]
    } else {
      default_reduc <- reduc_all
    }
  }
  if (length(default_reduc) > 1) {
    default_reduc <- default_reduc[unlist(sapply(c(pattern_default, pattern_dim), function(pat) {
      grep(pattern = pat, x = default_reduc, ignore.case = TRUE)
    }))]
    default_reduc <- default_reduc[which.min(sapply(default_reduc, function(x) dim(srt@reductions[[x]])[2]))]
  }
  return(default_reduc)
}

#' Uncorrected_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData SetAssayData VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData
#' @export
Uncorrected_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                                  do_normalization = NULL, normalization_method = "LogNormalize",
                                  do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                                  do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                                  linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                                  nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                  neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                                  seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Uncorrected) on the data...\n"))

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtMerge, layer = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    if (normalization_method != "SCT") {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srtMerge <- ScaleData(object = srtMerge, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srtMerge,
    prefix = "Uncorrected", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("Uncorrected", linear_reduction)]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtMerge <- tryCatch(
    {
      srtMerge <- FindNeighbors(
        object = srtMerge, reduction = paste0("Uncorrected", linear_reduction), dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("Uncorrected_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtMerge <- FindClusters(object = srtMerge, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "Uncorrected_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtMerge <- SrtReorder(srtMerge, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtMerge[["seurat_clusters"]] <- NULL
      srtMerge[[paste0("Uncorrected", linear_reduction, "clusters")]] <- Idents(srtMerge)
      srtMerge
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtMerge)
    }
  )

  srtMerge <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtMerge <- RunDimReduction(
            srtMerge,
            prefix = "Uncorrected",
            reduction_use = paste0("Uncorrected", linear_reduction), reduction_dims = linear_reduction_dims_use,
            graph_use = "Uncorrected_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtMerge
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtMerge)
    }
  )

  DefaultAssay(srtMerge) <- assay
  VariableFeatures(srtMerge) <- srtMerge@misc[["Uncorrected_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtMerge, pattern = paste0(assay, "|Uncorrected|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtMerge)
  }
}

#' Seurat_integrate
#'
#' @inheritParams Integration_SCP
#' @param FindIntegrationAnchors_params A list of parameters for the Seurat::FindIntegrationAnchors function, default is an empty list.
#' @param IntegrateData_params A list of parameters for the Seurat::IntegrateData function, default is an empty list.
#' @param IntegrateEmbeddings_params A list of parameters for the Seurat::IntegrateEmbeddings function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData FindIntegrationAnchors IntegrateData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents
#' @importFrom SeuratObject LayerData
#' @importFrom Signac RunTFIDF
#' @export
Seurat_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                             linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                             FindIntegrationAnchors_params = list(), IntegrateData_params = list(), IntegrateEmbeddings_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "svd", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    if (normalization_method == "TFIDF") {
      srtMerge <- Reduce(merge, srtList)
      VariableFeatures(srtMerge) <- HVF
    }
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srtList, ncol)) < 50) {
    warning("The cell count in some batches is lower than 50, which may not be suitable for the current integration method.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srtMerge)
    }
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'rlsi' integration workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      FindIntegrationAnchors_params[["dims"]] <- 2:min(linear_reduction_dims, 30)
    }
    srtMerge <- RunTFIDF(object = srtMerge, assay = DefaultAssay(srtMerge), verbose = FALSE)
    srtMerge <- RunDimReduction(
      srtMerge,
      prefix = "", features = HVF, assay = DefaultAssay(srtMerge),
      linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    srtMerge[["lsi"]] <- srtMerge[["svd"]]
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (svd) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt,
        prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srtList[[i]] <- srt
    }
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    cat(paste0("[", Sys.time(), "]", " Use 'rpca' integration workflow...\n"))
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srt, layer = "scale.data", assay = DefaultAssay(srt)))))) {
        cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
        srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
      }
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (pca) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt,
        prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "pca", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srtList[[i]] <- srt
    }
  }

  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  if (normalization_method %in% c("LogNormalize", "SCT")) {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = normalization_method,
      anchor.features = HVF,
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)

    cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
    params2 <- list(
      anchorset = srt_anchors,
      new.assay.name = "Seuratcorrected",
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateData, .args = params2)

    DefaultAssay(srtIntegrated) <- "Seuratcorrected"
    VariableFeatures(srtIntegrated[["Seuratcorrected"]]) <- HVF

    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtIntegrated, layer = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srtIntegrated <- ScaleData(object = srtIntegrated, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }

    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
    srtIntegrated <- RunDimReduction(
      srtIntegrated,
      prefix = "Seurat", features = HVF, assay = DefaultAssay(srtIntegrated),
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("Seurat", linear_reduction)]]@misc[["dims_estimate"]] %||% 1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)

    cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srtMerge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateEmbeddings, .args = params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:max(srtIntegrated@reductions[[paste0("Seurat", linear_reduction)]]@misc[["dims_estimate"]]) %||% 2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = paste0("Seurat", linear_reduction), dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("Seurat_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "Seurat_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Seurat",
            reduction_use = paste0("Seurat", linear_reduction), reduction_dims = linear_reduction_dims_use,
            graph_use = "Seurat_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Seurat|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' scVI_integrate
#'
#' @inheritParams Integration_SCP
#' @param scVI_dims_use A vector specifying the dimensions returned by scVI that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param model A string indicating the scVI model to be used. Options are "SCVI" and "PEAKVI". Default is "SCVI".
#' @param SCVI_params A list of parameters for the SCVI model, default is an empty list.
#' @param PEAKVI_params A list of parameters for the PEAKVI model, default is an empty list.
#' @param num_threads An integer setting the number of threads for scVI, default is 8.
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData<-
#' @importFrom SeuratObject LayerData
#' @importFrom reticulate import
#' @export
scVI_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                           do_normalization = NULL, normalization_method = "LogNormalize",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                           scVI_dims_use = NULL,
                           nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                           neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                           model = "SCVI", SCVI_params = list(), PEAKVI_params = list(), num_threads = 8, seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (.Platform$OS.type == "windows" && !exist_Python_pkgs(packages = "scvi-tools")) {
    suppressWarnings(system2(command = conda_python(), args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver", stdout = TRUE))
  }

  check_Python("scvi-tools")
  scvi <- import("scvi")
  scipy <- import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(num_threads)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  adata <- srt_to_adata(srtMerge, features = HVF, assay_X = DefaultAssay(srtMerge), assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  if (model == "SCVI") {
    scvi$model$SCVI$setup_anndata(adata, batch_key = batch)
    params <- list(
      adata = adata
    )
    for (nm in names(SCVI_params)) {
      params[[nm]] <- SCVI_params[[nm]]
    }
    model <- invoke(.fn = scvi$model$SCVI, .args = params)
    model$train()
    srtIntegrated <- srtMerge
    srtMerge <- NULL
    corrected <- t(as_matrix(model$get_normalized_expression()))
    srtIntegrated[["scVIcorrected"]] <- CreateAssayObject(counts = corrected)
    DefaultAssay(srtIntegrated) <- "scVIcorrected"
    VariableFeatures(srtIntegrated[["scVIcorrected"]]) <- HVF
  } else if (model == "PEAKVI") {
    message("Assay is ChromatinAssay. Using PeakVI workflow.")
    scvi$model$PEAKVI$setup_anndata(adata, batch_key = batch)
    params <- list(
      adata = adata
    )
    for (nm in names(PEAKVI_params)) {
      params[[nm]] <- PEAKVI_params[[nm]]
    }
    model <- invoke(.fn = scvi$model$PEAKVI, .args = params)
    model$train()
    srtIntegrated <- srtMerge
    srtMerge <- NULL
  }

  latent <- as_matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
  srtIntegrated[["scVI"]] <- CreateDimReducObject(embeddings = latent, key = "scVI_", assay = DefaultAssay(srtIntegrated))
  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srtIntegrated[["scVI"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "scVI", dims = scVI_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("scVI_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "scVI_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["scVIclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "scVI",
            reduction_use = "scVI", reduction_dims = scVI_dims_use,
            graph_use = "scVI_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["scVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|scVI|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' SCANVI_integrate
#'
#' @inheritParams Integration_SCP
#' @param SCANVI_dims_use A vector specifying the dimensions returned by scANVI that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param labels_key A character string indicating the column name in the metadata that stores cell labels for semi-supervised training. This argument is required.
#' @param unlabeled_category A string specifying the category in \code{labels_key} that represents unlabeled cells. Default is "Unknown".
#' @param SCVI_params A list of parameters passed to the \code{SCVI} model used for pretraining, default is an empty list.
#' @param SCANVI_params A list of parameters passed to \code{SCANVI\$from_scvi_model}, default is an empty list.
#' @param num_threads An integer setting the number of threads for scvi-tools, default is 8.
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData<-
#' @importFrom SeuratObject LayerData
#' @importFrom reticulate import
#' @export
SCANVI_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             SCANVI_dims_use = NULL,
                             labels_key = NULL, unlabeled_category = "Unknown",
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                             SCVI_params = list(), SCANVI_params = list(), num_threads = 8, seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (is.null(labels_key) || length(labels_key) != 1 || !is.character(labels_key)) {
    stop("'labels_key' must be a single character string specifying the metadata column with labels.")
  }

  if (.Platform$OS.type == "windows" && !exist_Python_pkgs(packages = "scvi-tools")) {
    suppressWarnings(system2(command = conda_python(), args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver", stdout = TRUE))
  }

  check_Python("scvi-tools")
  scvi <- import("scvi")
  scipy <- import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(num_threads)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  if (!labels_key %in% colnames(srtMerge@meta.data)) {
    stop(paste0("labels_key '", labels_key, "' was not found in the merged object's metadata."))
  }
  label_values <- srtMerge@meta.data[[labels_key]]
  if (!unlabeled_category %in% unique(label_values)) {
    stop(paste0("unlabeled_category '", unlabeled_category, "' was not found in column '", labels_key, "'."))
  }

  adata <- srt_to_adata(srtMerge, features = HVF, assay_X = DefaultAssay(srtMerge), assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  scvi$model$SCVI$setup_anndata(adata, batch_key = batch, labels_key = labels_key)
  scvi_params <- list(
    adata = adata
  )
  for (nm in names(SCVI_params)) {
    scvi_params[[nm]] <- SCVI_params[[nm]]
  }
  scvi_model <- invoke(.fn = scvi$model$SCVI, .args = scvi_params)
  scvi_model$train()

  scanvi_args <- list(
    scvi_model = scvi_model,
    unlabeled_category = unlabeled_category
  )
  for (nm in names(SCANVI_params)) {
    scanvi_args[[nm]] <- SCANVI_params[[nm]]
  }
  model <- invoke(.fn = scvi$model$SCANVI$from_scvi_model, .args = scanvi_args)
  model$train()

  srtIntegrated <- srtMerge
  srtMerge <- NULL

  corrected <- t(as_matrix(model$get_normalized_expression()))
  srtIntegrated[["SCANVIcorrected"]] <- CreateAssayObject(counts = corrected)
  DefaultAssay(srtIntegrated) <- "SCANVIcorrected"
  VariableFeatures(srtIntegrated[["SCANVIcorrected"]]) <- HVF

  latent <- as_matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("SCANVI_", seq_len(ncol(latent)))
  srtIntegrated[["SCANVI"]] <- CreateDimReducObject(embeddings = latent, key = "SCANVI_", assay = DefaultAssay(srtIntegrated))
  if (is.null(SCANVI_dims_use)) {
    SCANVI_dims_use <- 1:ncol(srtIntegrated[["SCANVI"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "SCANVI", dims = SCANVI_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("SCANVI_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "SCANVI_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["SCANVIclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "SCANVI",
            reduction_use = "SCANVI", reduction_dims = SCANVI_dims_use,
            graph_use = "SCANVI_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["SCANVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|SCANVI|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' MNN_integrate
#'
#' @inheritParams Integration_SCP
#' @param mnnCorrect_params A list of parameters for the batchelor::mnnCorrect function, default is an empty list.
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData<-
#' @importFrom SeuratObject LayerData
#' @export
MNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                          linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                          mnnCorrect_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  sceList <- lapply(srtList, function(srt) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = get_seurat_data(srt, layer = "data", assay = DefaultAssay(srt))[HVF, , drop = FALSE]))
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- as_matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(MNN) on the data...\n"))
  params <- list(
    sceList,
    cos.norm.out = FALSE
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- invoke(.fn = batchelor::mnnCorrect, .args = params)

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["MNNcorrected"]] <- CreateAssayObject(counts = out@assays@data$corrected)
  VariableFeatures(srtIntegrated[["MNNcorrected"]]) <- HVF
  DefaultAssay(srtIntegrated) <- "MNNcorrected"

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtIntegrated, layer = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(object = srtIntegrated, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srtIntegrated,
    prefix = "MNN", features = HVF, assay = DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("MNN", linear_reduction)]]@misc[["dims_estimate"]] %||% 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = paste0("MNN", linear_reduction), dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("MNN_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "MNN_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["MNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "MNN",
            reduction_use = paste0("MNN", linear_reduction), reduction_dims = linear_reduction_dims_use,
            graph_use = "MNN_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["MNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|MNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' fastMNN_integrate
#'
#' @inheritParams Integration_SCP
#' @param fastMNN_dims_use A vector specifying the dimensions returned by fastMNN that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param fastMNN_params A list of parameters for the batchelor::fastMNN function, default is an empty list.
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
fastMNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                              do_normalization = NULL, normalization_method = "LogNormalize",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                              fastMNN_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                              fastMNN_params = list(), seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  sceList <- lapply(srtList, function(srt) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = get_seurat_data(srt, layer = "data", assay = DefaultAssay(srt))[HVF, , drop = FALSE]))
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- as_matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(fastMNN) on the data...\n"))
  params <- list(
    sceList
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke(.fn = batchelor::fastMNN, .args = params)

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["fastMNNcorrected"]] <- CreateAssayObject(counts = as_matrix(out@assays@data$reconstructed))
  DefaultAssay(srtIntegrated) <- "fastMNNcorrected"
  VariableFeatures(srtIntegrated[["fastMNNcorrected"]]) <- HVF
  reduction <- out@int_colData$reducedDims$corrected
  colnames(reduction) <- paste0("fastMNN_", seq_len(ncol(reduction)))
  srtIntegrated[["fastMNN"]] <- CreateDimReducObject(embeddings = reduction, key = "fastMNN_", assay = "fastMNNcorrected")

  if (is.null(fastMNN_dims_use)) {
    fastMNN_dims_use <- 1:ncol(srtIntegrated[["fastMNN"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "fastMNN", dims = fastMNN_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("fastMNN", "_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "fastMNN_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["fastMNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "fastMNN",
            reduction_use = "fastMNN", reduction_dims = fastMNN_dims_use,
            graph_use = "fastMNN_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["fastMNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|fastMNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Harmony_integrate
#'
#' @inheritParams Integration_SCP
#' @param Harmony_dims_use A vector specifying the dimensions returned by RunHarmony that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param RunHarmony_params A list of parameters for the harmony::RunHarmony function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
Harmony_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                              do_normalization = NULL, normalization_method = "LogNormalize",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                              linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                              Harmony_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                              RunHarmony_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("harmony@1.1.0")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtMerge, layer = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srtMerge,
    prefix = "Harmony", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("Harmony", linear_reduction)]]@misc[["dims_estimate"]] %||% 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Harmony) on the data...\n"))
  message("Harmony integration using Reduction(", paste0("Harmony", linear_reduction), ", dims:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
  params <- list(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = DefaultAssay(srtMerge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  if (nrow(get_seurat_data(srtMerge, layer = "scale.data", assay = DefaultAssay(srtMerge))) == 0) {
    params[["project.dim"]] <- FALSE
  }
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srtIntegrated <- invoke(.fn = RunHarmony2, .args = params)

  if (is.null(Harmony_dims_use)) {
    Harmony_dims_use <- 1:ncol(srtIntegrated[["Harmony"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "Harmony", dims = Harmony_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("Harmony", "_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "Harmony_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Harmonyclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Harmony",
            reduction_use = "Harmony", reduction_dims = Harmony_dims_use,
            graph_use = "Harmony_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            force_nonlinear_reduction = force_nonlinear_reduction,
            nonlinear_reduction_params = nonlinear_reduction_params,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Harmony_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Harmony|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Scanorama_integrate
#'
#' @inheritParams Integration_SCP
#' @param Scanorama_dims_use  A vector specifying the dimensions returned by Scanorama that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param return_corrected Logical indicating whether to return the corrected data. Default is FALSE.
#' @param Scanorama_params A list of parameters for the scanorama.correct function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @importFrom SeuratObject LayerData
#' @importFrom Matrix t
#' @importFrom reticulate import
#' @importFrom stats sd
#' @export
Scanorama_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                                do_normalization = NULL, normalization_method = "LogNormalize",
                                do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                                do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                Scanorama_dims_use = NULL,
                                nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                                return_corrected = FALSE, Scanorama_params = list(), seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_Python("scanorama")
  scanorama <- import("scanorama")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  srtIntegrated <- Reduce(merge, srtList)

  cat(paste0("[", Sys.time(), "]", " Perform integration(Scanorama) on the data...\n"))
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srtList)) {
    assaylist[[i]] <- t(as_matrix(get_seurat_data(srtList[[i]], layer = "data", assay = DefaultAssay(srtList[[i]]))[HVF, , drop = FALSE])) # test
    genelist[[i]] <- HVF
  }
  if (isTRUE(return_corrected)) {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      return_dimred = TRUE,
      return_dense = TRUE,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    corrected <- invoke(.fn = scanorama$correct, .args = params)

    cor_value <- t(do.call(rbind, corrected[[2]]))
    rownames(cor_value) <- corrected[[3]]
    colnames(cor_value) <- unlist(sapply(assaylist, rownames))
    srtIntegrated[["Scanoramacorrected"]] <- CreateAssayObject(data = cor_value)
    VariableFeatures(srtIntegrated[["Scanoramacorrected"]]) <- HVF

    dim_reduction <- do.call(rbind, corrected[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0("Scanorama_", seq_len(ncol(dim_reduction)))
  } else {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    integrated <- invoke(.fn = scanorama$integrate, .args = params)

    dim_reduction <- do.call(rbind, integrated[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0("Scanorama_", seq_len(ncol(dim_reduction)))
  }
  srtIntegrated[["Scanorama"]] <- CreateDimReducObject(embeddings = dim_reduction, key = "Scanorama_", assay = DefaultAssay(srtIntegrated))

  if (is.null(Scanorama_dims_use)) {
    Scanorama_dims_use <- 1:ncol(srtIntegrated[["Scanorama"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "Scanorama", dims = Scanorama_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("Scanorama_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "Scanorama_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Scanoramaclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Scanorama",
            reduction_use = "Scanorama", reduction_dims = Scanorama_dims_use,
            graph_use = "Scanorama_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Scanorama_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Scanorama|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' BBKNN_integrate
#'
#' @inheritParams Integration_SCP
#' @param bbknn_params A list of parameters for the bbknn.matrix.bbknn function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- as.Graph Embeddings FindClusters Idents VariableFeatures VariableFeatures<- as.sparse
#' @importFrom SeuratObject LayerData
#' @importFrom Matrix t
#' @importFrom reticulate import
#' @export
BBKNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            cluster_algorithm = "louvain", cluster_resolution = 0.6,
                            bbknn_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_Python("bbknn")
  bbknn <- import("bbknn")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtMerge, layer = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srtMerge,
    prefix = "BBKNN", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("BBKNN", linear_reduction)]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(BBKNN) on the data...\n"))
  message("BBKNN integration using Reduction(", paste0("BBKNN", linear_reduction), ", dims:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
  emb <- Embeddings(srtMerge, reduction = paste0("BBKNN", linear_reduction))[, linear_reduction_dims_use, drop = FALSE]
  params <- list(
    pca = emb,
    batch_list = srtMerge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke(.fn = bbknn$matrix$bbknn, .args = params)
  n.neighbors <- bem[[3]]$n_neighbors
  srtIntegrated <- srtMerge

  bbknn_graph <- as.sparse(bem[[2]][1:nrow(bem[[2]]), , drop = FALSE])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- t(as.sparse(bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]))
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["BBKNN_dist"]] <- bbknn_dist

  val <- split(bbknn_dist@x, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  pos <- split(bbknn_dist@i + 1, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  idx <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = F), n.neighbors - 1)]
    length(out) <- n.neighbors - 1
    return(out)
  }, x = val, y = pos))
  idx[is.na(idx)] <- sample(seq_len(nrow(idx)), size = sum(is.na(idx)), replace = TRUE)
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = F), n.neighbors - 1)]
    length(out) <- n.neighbors - 1
    out[is.na(out)] <- 0
    return(out)
  }, x = val, y = val))
  dist <- cbind(0, dist)
  srtIntegrated[["BBKNN_neighbors"]] <- new(Class = "Neighbor", nn.idx = idx, nn.dist = dist, alg.info = list(), cell.names = rownames(emb))
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

  srtIntegrated <- tryCatch(
    {
      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "BBKNN", resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["BBKNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors
        }
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "BBKNN", neighbor_use = "BBKNN_neighbors",
            graph_use = "BBKNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["BBKNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|BBKNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' CSS_integrate
#'
#' @inheritParams Integration_SCP
#' @param CSS_dims_use A vector specifying the dimensions returned by CSS that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param CSS_params A list of parameters for the simspec::cluster_sim_spectrum function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData
#' @export
CSS_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                          linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          CSS_dims_use = NULL,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                          CSS_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R(c("quadbiolab/simspec", "qlcMatrix"))
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtMerge, layer = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srtMerge,
    prefix = "CSS", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("CSS", linear_reduction)]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(CSS) on the data...\n"))
  message("CSS integration using Reduction(", paste0("CSS", linear_reduction), ", dims:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
  params <- list(
    object = srtMerge,
    use_dr = paste0("CSS", linear_reduction),
    dims_use = linear_reduction_dims_use,
    var_genes = HVF,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  for (nm in names(CSS_params)) {
    params[[nm]] <- CSS_params[[nm]]
  }
  srtIntegrated <- invoke(.fn = get("cluster_sim_spectrum", envir = getNamespace("simspec")), .args = params)

  if (any(is.na(srtIntegrated@reductions[["CSS"]]@cell.embeddings))) {
    stop("NA detected in the CSS embeddings. You can try to use a lower resolution value in the CSS_param.")
  }
  if (is.null(CSS_dims_use)) {
    CSS_dims_use <- 1:ncol(srtIntegrated[["CSS"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "CSS", dims = CSS_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("CSS", "_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "CSS_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["CSSclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "CSS",
            reduction_use = "CSS", reduction_dims = CSS_dims_use,
            graph_use = "CSS_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["CSS_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|CSS|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' LIGER_integrate
#'
#' @inheritParams Integration_SCP
#' @param LIGER_dims_use A vector specifying the dimensions returned by LIGER that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param optimizeALS_params A list of parameters for the rliger::optimizeALS function, default is an empty list.
#' @param quantilenorm_params A list of parameters for the rliger::quantile_norm function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData
#' @export
LIGER_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            LIGER_dims_use = NULL,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                            optimizeALS_params = list(), quantilenorm_params = list(), seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("rliger")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srtList, ncol)) < 30) {
    warning("The cell count in some batches is lower than 30, which may not be suitable for the current integration method.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srtMerge)
    }
  }

  scale.data <- list()
  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srt, layer = "scale.data", assay = DefaultAssay(srt)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, do.center = FALSE, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    scale.data[[i]] <- t(x = get_seurat_data(srt, layer = "scale.data", assay = DefaultAssay(srt)))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(LIGER) on the data...\n"))
  params1 <- list(
    object = scale.data,
    k = 20,
    verbose = FALSE
  )
  for (nm in names(optimizeALS_params)) {
    params1[[nm]] <- optimizeALS_params[[nm]]
  }
  out1 <- invoke(.fn = rliger::optimizeALS, .args = params1)
  cat("\n")
  colnames(x = out1$W) <- colnames(scale.data[[1]])
  reduction1 <- do.call(what = "rbind", args = out1$H)
  colnames(reduction1) <- paste0("riNMF_", seq_len(ncol(reduction1)))
  loadings1 <- t(x = out1$W)
  rownames(loadings1) <- colnames(scale.data[[1]])
  colnames(loadings1) <- paste0("riNMF_", seq_len(ncol(loadings1)))
  srtMerge[["iNMF_raw"]] <- CreateDimReducObject(
    embeddings = reduction1,
    loadings = loadings1,
    assay = DefaultAssay(srtMerge),
    key = "riNMF_"
  )

  embeddings <- sapply(
    X = SplitObject(object = srtMerge, split.by = batch),
    FUN = function(x) {
      return(Embeddings(object = x[["iNMF_raw"]]))
    }, simplify = FALSE, USE.NAMES = TRUE
  )
  num.samples <- vapply(X = embeddings, FUN = nrow, FUN.VALUE = integer(length = 1L))
  ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  params2 <- list(
    object = embeddings,
    ref_dataset = ref_dataset
  )
  for (nm in names(quantilenorm_params)) {
    params2[[nm]] <- quantilenorm_params[[nm]]
  }
  out2 <- invoke(.fn = rliger::quantile_norm, .args = params2)
  srtMerge[["LIGER"]] <- CreateDimReducObject(
    embeddings = out2$H.norm,
    assay = DefaultAssay(srtMerge),
    key = "LIGER_"
  )
  srtIntegrated <- srtMerge
  srtMerge <- NULL
  if (is.null(LIGER_dims_use)) {
    LIGER_dims_use <- 1:ncol(srtIntegrated[["LIGER"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = "LIGER", dims = LIGER_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("LIGER", "_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "LIGER_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["LIGERclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "LIGER",
            reduction_use = "LIGER", reduction_dims = LIGER_dims_use,
            graph_use = "LIGER_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["LIGER_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|LIGER|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Conos_integrate
#'
#' @inheritParams Integration_SCP
#' @param buildGraph_params A list of parameters for the buildGraph function, default is an empty list.
#' @param num_threads  An integer setting the number of threads for Conos, default is 2.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData
#' @importFrom igraph as_adjacency_matrix
#' @export
Conos_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            cluster_algorithm = "louvain", cluster_resolution = 0.6,
                            buildGraph_params = list(), num_threads = 2, seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("conos")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srtList, ncol)) < 30) {
    warning("The cell count in some batches is lower than 30, which may not be suitable for the current integration method.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srtMerge)
    }
  }

  srtIntegrated <- srtMerge
  srtMerge <- NULL

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srt, layer = "scale.data", assay = DefaultAssay(srt)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data ", i, " ...\n"))
    srt <- RunDimReduction(
      srt,
      prefix = "Conos", features = HVF, assay = DefaultAssay(srt),
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    srt[["pca"]] <- srt[[paste0("Conos", linear_reduction)]]
    srtList[[i]] <- srt
  }
  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  if (is.null(linear_reduction_dims_use)) {
    maxdims <- max(unlist(sapply(srtList, function(srt) max(srt@reductions[[paste0("Conos", linear_reduction)]]@misc[["dims_estimate"]]))))
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Conos) on the data...\n"))
  message("Conos integration using Reduction(", linear_reduction, ", dims_max:", maxdims, ") as input")
  srtList_con <- conos::Conos$new(srtList, n.cores = num_threads)
  params <- list(
    ncomps = maxdims,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invoke(.fn = srtList_con[["buildGraph"]], .args = params)
  conos_graph <- as_adjacency_matrix(srtList_con$graph, type = "both", attr = "weight", names = TRUE, sparse = TRUE)
  conos_graph <- as.Graph(conos_graph)
  conos_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

  srtIntegrated <- tryCatch(
    {
      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "Conos", resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Conosclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]
        }
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Conos", graph_use = "Conos",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Conos_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Conos|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Combat_integrate
#'
#' @inheritParams Integration_SCP
#' @param ComBat_params A list of parameters for the sva::ComBat function, default is an empty list.
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @importFrom SeuratObject LayerData
#' @export
ComBat_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                             linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                             ComBat_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("sva")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Combat) on the data...\n"))
  dat <- get_seurat_data(srtMerge, layer = "data", assay = DefaultAssay(srtMerge))[HVF, , drop = FALSE]
  batch <- srtMerge[[batch, drop = TRUE]]
  params <- list(
    dat = dat,
    batch = batch
  )
  for (nm in names(ComBat_params)) {
    params[[nm]] <- ComBat_params[[nm]]
  }
  corrected <- suppressMessages(invoke(.fn = sva::ComBat, .args = params))

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["ComBatcorrected"]] <- CreateAssayObject(data = corrected)
  DefaultAssay(srtIntegrated) <- "ComBatcorrected"
  VariableFeatures(srtIntegrated[["ComBatcorrected"]]) <- HVF

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srtIntegrated, layer = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(srtIntegrated, split.by = if (isTRUE(scale_within_batch)) batch else NULL, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srtIntegrated,
    prefix = "ComBat", features = HVF, assay = DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("ComBat", linear_reduction)]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = paste0("ComBat", linear_reduction), dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        force.recalc = TRUE, graph.name = paste0("ComBat_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "ComBat_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["ComBatclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "ComBat",
            reduction_use = paste0("ComBat", linear_reduction), reduction_dims = linear_reduction_dims_use,
            graph_use = "ComBat_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["ComBat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|ComBat|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Standard SCP
#'
#' This function performs a standard single-cell analysis workflow.
#'
#' @param srt A Seurat object.
#' @param prefix A prefix to add to the names of intermediate objects created by the function (default is "Standard").
#' @param assay The name of the assay to use for the analysis. If NULL, the default assay of the Seurat object will be used.
#' @param do_normalization A logical value indicating whether to perform normalization. If NULL, normalization will be performed if the specified assay does not have scaled data.
#' @param normalization_method The method to use for normalization. Options are "LogNormalize", "SCT", or "TFIDF" (default is "LogNormalize").
#' @param do_HVF_finding A logical value indicating whether to perform high variable feature finding. If TRUE, the function will force to find the highly variable features (HVF) using the specified HVF method.
#' @param HVF_method The method to use for finding highly variable features. Options are "vst", "mvp" or "disp" (default is "vst").
#' @param nHVF The number of highly variable features to select. If NULL, all highly variable features will be used.
#' @param HVF A vector of feature names to use as highly variable features. If NULL, the function will use the highly variable features identified by the HVF method.
#' @param do_scaling A logical value indicating whether to perform scaling. If TRUE, the function will force to scale the data using the ScaleData function.
#' @param vars_to_regress A vector of feature names to use as regressors in the scaling step. If NULL, no regressors will be used.
#' @param regression_model The regression model to use for scaling. Options are "linear", "poisson", or "negativebinomial" (default is "linear").
#' @param linear_reduction The linear dimensionality reduction method to use. Options are "pca", "svd", "ica", "nmf", "mds", or "glmpca" (default is "pca").
#' @param linear_reduction_dims The number of dimensions to keep after linear dimensionality reduction (default is 50).
#' @param linear_reduction_dims_use The dimensions to use for downstream analysis. If NULL, all dimensions will be used.
#' @param linear_reduction_params A list of parameters to pass to the linear dimensionality reduction method.
#' @param force_linear_reduction A logical value indicating whether to force linear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' @param nonlinear_reduction The nonlinear dimensionality reduction method to use. Options are "umap","umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", or "fr" (default is "umap").
#' @param nonlinear_reduction_dims The number of dimensions to keep after nonlinear dimensionality reduction. If a vector is provided, different numbers of dimensions can be specified for each method (default is c(2, 3)).
#' @param nonlinear_reduction_params A list of parameters to pass to the nonlinear dimensionality reduction method.
#' @param force_nonlinear_reduction A logical value indicating whether to force nonlinear dimensionality reduction even if the specified reduction is already present in the Seurat object.
#' @param neighbor_metric The distance metric to use for finding neighbors. Options are "euclidean", "cosine", "manhattan", or "hamming" (default is "euclidean").
#' @param neighbor_k The number of nearest neighbors to use for finding neighbors (default is 20).
#' @param cluster_algorithm The clustering algorithm to use. Options are "louvain", "slm", or "leiden" (default is "louvain").
#' @param cluster_resolution The resolution parameter to use for clustering. Larger values result in fewer clusters (default is 0.6).
#' @param seed The random seed to use for reproducibility (default is 11).
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link{Integration_SCP}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType")
#'
#' # Use a combination of different linear or non-linear dimension reduction methods
#' linear_reductions <- c("pca", "ica", "nmf", "mds", "glmpca")
#' pancreas_sub <- Standard_SCP(
#'   pancreas_sub,
#'   linear_reduction = linear_reductions,
#'   nonlinear_reduction = "umap"
#' )
#' plist1 <- lapply(linear_reductions, function(lr) {
#'   CellDimPlot(pancreas_sub,
#'     group.by = "SubCellType",
#'     reduction = paste0("Standard", lr, "UMAP2D"),
#'     xlab = "", ylab = "", title = lr,
#'     legend.position = "none",
#'     theme_use = "theme_blank"
#'   )
#' })
#' patchwork::wrap_plots(plotlist = plist1)
#'
#' nonlinear_reductions <- c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr")
#' pancreas_sub <- Standard_SCP(
#'   pancreas_sub,
#'   linear_reduction = "pca",
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' plist2 <- lapply(nonlinear_reductions, function(nr) {
#'   CellDimPlot(pancreas_sub,
#'     group.by = "SubCellType",
#'     reduction = paste0("Standardpca", toupper(nr), "2D"),
#'     xlab = "", ylab = "", title = nr,
#'     legend.position = "none",
#'     theme_use = "theme_blank"
#'   )
#' })
#' patchwork::wrap_plots(plotlist = plist2)
#'
#' @importFrom Seurat Assays GetAssayData NormalizeData SCTransform SCTResults ScaleData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject LayerData
#' @importFrom Matrix rowSums
#' @export
Standard_SCP <- function(srt, prefix = "Standard", assay = NULL,
                         do_normalization = NULL, normalization_method = "LogNormalize",
                         do_HVF_finding = TRUE, HVF_method = "vst", nHVF = 2000, HVF = NULL,
                         do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                         linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                         nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                         neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                         seed = 11) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  if (any(!linear_reduction %in% c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt)))) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  time_start <- Sys.time()
  set.seed(seed)

  cat(paste0("[", time_start, "] ", "Start Standard_SCP\n"))

  checked <- check_srtList(
    srtList = list(srt), batch = "", assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = "separate", HVF_method = HVF_method, nHVF = nHVF, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srt <- checked[["srtList"]][[1]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(get_seurat_data(srt, layer = "scale.data", assay = DefaultAssay(srt)))))) {
    if (normalization_method != "SCT") {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
  }

  for (lr in linear_reduction) {
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", lr, ") on the data...\n"))
    srt <- RunDimReduction(
      srt,
      prefix = prefix, features = HVF, assay = DefaultAssay(srt),
      linear_reduction = lr, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use_current <- srt@reductions[[paste0(prefix, lr)]]@misc[["dims_estimate"]]
      if (normalization_method == "TFIDF") {
        linear_reduction_dims_use_current <- 2:max(linear_reduction_dims_use_current)
      }
    } else {
      linear_reduction_dims_use_current <- linear_reduction_dims_use
    }

    srt <- tryCatch(
      {
        srt <- FindNeighbors(
          object = srt, reduction = paste0(prefix, lr), dims = linear_reduction_dims_use_current,
          annoy.metric = neighbor_metric, k.param = neighbor_k,
          force.recalc = TRUE, graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")), verbose = FALSE
        )

        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srt <- FindClusters(object = srt, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0(prefix, lr, "_SNN"), verbose = FALSE)
        cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
        srt <- SrtReorder(srt, features = HVF, reorder_by = "seurat_clusters", layer = "data")
        srt[["seurat_clusters"]] <- NULL
        srt[[paste0(prefix, lr, "clusters")]] <- Idents(srt)
        srt
      },
      error = function(error) {
        message(error)
        message("Error when performing FindClusters. Skip this step...")
        return(srt)
      }
    )

    srt <- tryCatch(
      {
        for (nr in nonlinear_reduction) {
          cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
          for (n in nonlinear_reduction_dims) {
            srt <- RunDimReduction(
              srt,
              prefix = paste0(prefix, lr),
              reduction_use = paste0(prefix, lr), reduction_dims = linear_reduction_dims_use_current,
              graph_use = paste0(prefix, lr, "_SNN"),
              nonlinear_reduction = nr, nonlinear_reduction_dims = n,
              nonlinear_reduction_params = nonlinear_reduction_params,
              force_nonlinear_reduction = force_nonlinear_reduction,
              verbose = FALSE, seed = seed
            )
          }
        }
        srt
      },
      error = function(error) {
        message(error)
        message("Error when performing nonlinear dimension reduction. Skip this step...")
        return(srt)
      }
    )
  }

  if (paste0(prefix, linear_reduction[1], "clusters") %in% colnames(srt@meta.data)) {
    srt[[paste0(prefix, "clusters")]] <- srt[[paste0(prefix, linear_reduction[1], "clusters")]]
  }
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      if (paste0(prefix, linear_reduction[1], toupper(nr), n, "D") %in% names(srt@reductions)) {
        reduc <- srt@reductions[[paste0(prefix, linear_reduction[1], toupper(nr), n, "D")]]
        srt@reductions[[paste0(prefix, toupper(nr), n, "D")]] <- reduc
      }
    }
    srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nr))
  }

  DefaultAssay(srt) <- assay
  VariableFeatures(srt) <- srt@misc[["Standard_HVF"]] <- HVF

  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "Standard_SCP done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

  return(srt)
}

#' Integration_SCP
#'
#' Integrate single-cell RNA-seq data using various integration methods.
#'
#' @inheritParams check_srtList
#' @inheritParams check_srtMerge
#' @inheritParams Standard_SCP
#' @param scale_within_batch  Whether to scale data within each batch. Only valid when the \code{integration_method} is one of \code{"Uncorrected"}, \code{"Seurat"}, \code{"MNN"}, \code{"Harmony"}, \code{"BBKNN"}, \code{"CSS"}, \code{"ComBat"}.
#' @param integration_method  A character string specifying the integration method to use.
#'   Supported methods are: \code{"Uncorrected"}, \code{"Seurat"}, \code{"scVI"}, \code{"SCANVI"}, \code{"MNN"}, \code{"fastMNN"}, \code{"Harmony"},
#'   \code{"Scanorama"}, \code{"BBKNN"}, \code{"CSS"}, \code{"LIGER"}, \code{"Conos"}, \code{"ComBat"}. Default is \code{"Uncorrected"}.
#' @param append Logical, if TRUE, the integrated data will be appended to the original Seurat object (srtMerge).
#' @param ... Additional arguments to be passed to the integration method function.
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link{Seurat_integrate}} \code{\link{scVI_integrate}} \code{\link{SCANVI_integrate}} \code{\link{MNN_integrate}} \code{\link{fastMNN_integrate}} \code{\link{Harmony_integrate}} \code{\link{Scanorama_integrate}} \code{\link{BBKNN_integrate}} \code{\link{CSS_integrate}} \code{\link{LIGER_integrate}} \code{\link{Conos_integrate}} \code{\link{ComBat_integrate}} \code{\link{Standard_SCP}}
#'
#' @examples
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5, scale_within_batch = TRUE
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat",
#'   FindIntegrationAnchors_params = list(reduction = "rpca")
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' \dontrun{
#' integration_methods <- c(
#'   "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony",
#'   "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat"
#' )
#' for (method in integration_methods) {
#'   panc8_sub <- Integration_SCP(
#'     srtMerge = panc8_sub, batch = "tech",
#'     integration_method = method,
#'     linear_reduction_dims_use = 1:50,
#'     nonlinear_reduction = "umap"
#'   )
#'   print(CellDimPlot(panc8_sub,
#'     group.by = c("tech", "celltype"),
#'     reduction = paste0(method, "UMAP2D"),
#'     xlab = "", ylab = "", title = method,
#'     legend.position = "none", theme_use = "theme_blank"
#'   ))
#' }
#'
#' nonlinear_reductions <- c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr")
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat",
#'   linear_reduction_dims_use = 1:50,
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' for (nr in nonlinear_reductions) {
#'   print(CellDimPlot(panc8_sub,
#'     group.by = c("tech", "celltype"),
#'     reduction = paste0("Seurat", nr, "2D"),
#'     xlab = "", ylab = "", title = nr,
#'     legend.position = "none", theme_use = "theme_blank"
#'   ))
#' }
#' }
#'
#' @export
Integration_SCP <- function(srtMerge = NULL, batch, append = TRUE, srtList = NULL, assay = NULL,
                            integration_method = "Uncorrected",
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear", scale_within_batch = FALSE,
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            neighbor_metric = "euclidean", neighbor_k = 20L, cluster_algorithm = "louvain", cluster_resolution = 0.6,
                            seed = 11, ...) {
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("Neither 'srtList' nor 'srtMerge' was found.")
  }
  if (length(integration_method) == 1 && integration_method %in% c("Uncorrected", "Seurat", "scVI", "SCANVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat")) {
    # Convert the arguments of the function call to a list and remove the function itself
    args <- as.list(match.call())[-1]

    # Create a new environment, the parent of which is the environment that called the foo function
    new_env <- new.env(parent = parent.frame())

    # Evaluate the arguments in the new environment to get the correct values
    args <- lapply(args, function(x) eval(x, envir = new_env))

    # Keep srtMerge and srtList as type of 'symbol' when use `do.call` function
    # args[!names(args) %in% c("srtMerge", "srtList")] <- lapply(args[!names(args) %in% c("srtMerge", "srtList")], function(x) eval(x, envir = new_env))

    # print("================ args ================ ")
    # print(args)

    # Get the function's formal arguments and their default values
    formals <- mget(names(formals()))
    formals <- formals[names(formals) != "..."]

    # print("================ formals ================ ")
    # print(formals)

    # Merge the formal arguments with the actual arguments, so that all arguments are included
    args <- modifyList(formals, args)

    time_start <- Sys.time()
    cat(paste0("[", time_start, "] ", paste0("Start ", integration_method, "_integrate"), "\n"))
    srtIntegrated <- invoke(
      .fn = paste0(integration_method, "_integrate"),
      .args = args[names(args) %in% formalArgs(paste0(integration_method, "_integrate"))]
    )
    time_end <- Sys.time()
    cat(paste0("[", time_end, "] ", paste0(integration_method, "_integrate done\n")))
    cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

    return(srtIntegrated)
  } else {
    stop(paste(integration_method, "is not a supported integration method!"))
  }
}
