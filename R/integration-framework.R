# Integration Framework - Common Helper Functions
# These functions provide shared functionality for all integration methods
# Extracted from duplicated code across 13 integration functions

#' Validate Reduction Parameters (Internal)
#'
#' Common validation logic for linear/nonlinear reduction and clustering parameters.
#' Handles leiden algorithm setup if needed.
#'
#' @param linear_reduction Character. Linear reduction method(s).
#' @param nonlinear_reduction Character. Nonlinear reduction method(s).
#' @param cluster_algorithm Character. Clustering algorithm to use.
#' @param linear_reduction_dims Integer. Number of dimensions for linear reduction.
#' @param linear_reduction_dims_use Integer vector. Specific dimensions to use (optional).
#' @param srtMerge Seurat object (optional, for checking available reductions).
#'
#' @return List with validated parameters including cluster_algorithm_index
#' @keywords internal
.validate_reduction_params <- function(linear_reduction = "pca",
                                      nonlinear_reduction = "umap",
                                      cluster_algorithm = "louvain",
                                      linear_reduction_dims = 50,
                                      linear_reduction_dims_use = NULL,
                                      srtMerge = NULL) {

  # Validate linear reduction
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

  # Adjust linear_reduction_dims if needed
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }

  # Validate nonlinear reduction
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'.")
  }

  # Validate cluster algorithm
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }

  # Setup leiden if needed
  if (cluster_algorithm == "leiden") {
    # Leiden algorithm requires Python environment
    if (uv_env_exists()) {
      tryCatch(use_uv_env(), error = function(e) {
        stop("Failed to configure Python environment for leiden algorithm: ", e$message)
      })
    } else {
      stop("Leiden algorithm requires Python environment. Run PrepareEnv() first.")
    }
    use_uv_env()
    if (!reticulate::py_module_available("leidenalg")) {
      stop("Python module 'leidenalg' is required. Install with: uv_install(packages = 'leidenalg')")
    }
  }

  # Get cluster algorithm index for Seurat
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  return(list(
    linear_reduction = linear_reduction,
    nonlinear_reduction = nonlinear_reduction,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_dims_use = linear_reduction_dims_use
  ))
}


#' Post-Integration Clustering (Internal)
#'
#' Common clustering workflow: FindNeighbors → FindClusters → SrtReorder.
#' Used after integration to cluster cells based on the integrated reduction.
#'
#' @param srt Seurat object after integration
#' @param prefix Character. Prefix for naming (e.g., "Uncorrected", "Seurat", "Harmony")
#' @param reduction_name Character. Name of the reduction to use for clustering
#' @param dims Integer vector. Dimensions to use from the reduction
#' @param neighbor_metric Character. Distance metric for FindNeighbors
#' @param neighbor_k Integer. Number of neighbors
#' @param cluster_algorithm_index Integer. Algorithm index (1=louvain, 2=louvain_refined, 3=slm, 4=leiden)
#' @param cluster_resolution Numeric. Resolution for clustering
#' @param HVF Character vector. Highly variable features for SrtReorder
#' @param cluster_algorithm Character. Name of algorithm (for messages)
#' @param linear_reduction Character. Name of linear reduction method
#'
#' @return Seurat object with clustering results
#' @keywords internal
.post_integration_clustering <- function(srt, prefix, reduction_name, dims,
                                        neighbor_metric = "euclidean", neighbor_k = 20L,
                                        cluster_algorithm_index = 1, cluster_resolution = 0.6,
                                        HVF = NULL, cluster_algorithm = "louvain",
                                        linear_reduction = "pca") {

  srt <- tryCatch(
    {
      # Find neighbors
      srt <- FindNeighbors(
        object = srt, reduction = reduction_name, dims = dims,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        graph.name = paste0(prefix, "_", c("KNN", "SNN")), verbose = FALSE
      )

      # Find clusters
      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srt <- FindClusters(
        object = srt, resolution = cluster_resolution,
        algorithm = cluster_algorithm_index, method = "igraph",
        graph.name = paste0(prefix, "_SNN"), verbose = FALSE
      )

      # Reorder clusters
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srt <- SrtReorder(srt, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srt[["seurat_clusters"]] <- NULL
      srt[[paste0(prefix, linear_reduction, "clusters")]] <- Idents(srt)

      srt
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srt)
    }
  )

  return(srt)
}


#' Post-Integration Nonlinear Reductions (Internal)
#'
#' Common nonlinear reduction workflow: RunDimReduction for each method and dimension.
#' Used after integration to create visualizations (UMAP, tSNE, etc.).
#'
#' @param srt Seurat object after integration and clustering
#' @param prefix Character. Prefix for naming
#' @param reduction_use Character. Name of the reduction to use as input
#' @param reduction_dims Integer vector. Dimensions to use from the input reduction
#' @param graph_use Character. Name of the SNN graph to use
#' @param nonlinear_reduction Character vector. Nonlinear methods to compute
#' @param nonlinear_reduction_dims Integer vector. Dimensions for each nonlinear reduction
#' @param nonlinear_reduction_params List. Additional parameters
#' @param force_nonlinear_reduction Logical. Force recomputation
#' @param seed Integer. Random seed
#'
#' @return Seurat object with nonlinear reductions
#' @keywords internal
.post_integration_reductions <- function(srt, prefix, reduction_use, reduction_dims,
                                        graph_use,
                                        nonlinear_reduction = "umap",
                                        nonlinear_reduction_dims = c(2, 3),
                                        nonlinear_reduction_params = list(),
                                        force_nonlinear_reduction = TRUE,
                                        seed = 11) {

  srt <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srt <- RunDimReduction(
            srt,
            prefix = prefix,
            reduction_use = reduction_use,
            reduction_dims = reduction_dims,
            graph_use = graph_use,
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
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

  return(srt)
}
