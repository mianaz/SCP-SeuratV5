#' Run imputation methods on Seurat object
#'
#' This function performs imputation on a Seurat object using various methods including ALRA, MAGIC, etc.
#' The imputed values are stored as a new assay in the Seurat object.
#'
#' @param srt An object of class Seurat.
#' @param assay Name of assay to use for imputation.
#' @param slot Name of slot to use for imputation (default: "data").
#' @param method Imputation method to use. One of "alra", "magic", "knn_smooth".
#' @param new_assay Name of the new assay to store imputed data (default: paste0(method, "_imputed")).
#' @param features Features to use for imputation. If NULL, uses variable features.
#' @param verbose Print progress messages.
#'
#' @return A Seurat object with imputed values stored in a new assay.
#'
#' @importFrom Seurat CreateAssayObject DefaultAssay VariableFeatures
#' @importFrom SeuratObject LayerData
#' @export
RunImputation <- function(srt, assay = NULL, slot = "data", 
                         method = c("alra", "magic", "knn_smooth"),
                         new_assay = NULL, features = NULL,
                         verbose = TRUE) {
  
  if (!inherits(srt, "Seurat")) {
    stop("'srt' must be a Seurat object")
  }
  
  assay <- assay %||% DefaultAssay(srt)
  method <- match.arg(method)
  new_assay <- new_assay %||% paste0(method, "_imputed")
  
  is_v5 <- IsSeurat5(srt)
  if (is_v5) {
    data_mat <- LayerData(srt[[assay]], layer = slot)
  } else {
    data_mat <- GetAssayData(srt, slot = slot, assay = assay)
  }
  
  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
    if (length(features) == 0) {
      stop("No variable features found. Please run FindVariableFeatures() first or specify features.")
    }
  }
  data_mat <- data_mat[features, ]
  
  if (verbose) message("Running ", method, " imputation...")
  
  imputed_mat <- switch(method,
    "alra" = RunALRA(data_mat, verbose = verbose),
    "magic" = RunMAGIC(data_mat, verbose = verbose),
    "knn_smooth" = RunKNNSmooth(data_mat, verbose = verbose)
  )
  
  if (verbose) message("Creating new assay: ", new_assay)
  
  if (is_v5) {
    new_assay_obj <- CreateAssayObject(counts = NULL)
    new_assay_obj <- SetAssayData(new_assay_obj, slot = "data", 
                                 new.data = imputed_mat)
  } else {
    new_assay_obj <- CreateAssayObject(data = imputed_mat)
  }
  
  srt[[new_assay]] <- new_assay_obj
  
  return(srt)
}

#' Run ALRA imputation
#'
#' This function performs ALRA (Adaptively-thresholded Low Rank Approximation) imputation
#' on a gene expression matrix.
#'
#' @param data_mat Expression matrix with genes as rows and cells as columns
#' @param k Number of dimensions for SVD (default: NULL, automatically determined)
#' @param q Quantile for automatic k selection (default: 0.01)
#' @param verbose Print progress messages
#'
#' @return Imputed expression matrix
#'
#' @importFrom irlba irlba
#' @importFrom Matrix t
RunALRA <- function(data_mat, k = NULL, q = 0.01, verbose = TRUE) {
  require_packages("irlba")
  
  if (verbose) message("Centering and scaling data...")
  data_centered <- scale(t(data_mat), center = TRUE, scale = TRUE)
  
  if (verbose) message("Performing SVD...")
  svd_res <- irlba::irlba(data_centered, nv = k)
  
  if (is.null(k)) {
    k <- choose_k(svd_res$d, q = q)
    if (verbose) message("Automatically selected k = ", k)
  }
  
  if (verbose) message("Computing low-rank approximation...")
  A_norm <- svd_res$u[, 1:k] %*% diag(svd_res$d[1:k]) %*% t(svd_res$v[, 1:k])
  
  if (verbose) message("Completing missing values...")
  imputed_mat <- t(A_norm)
  
  return(imputed_mat)
}

#' Run MAGIC imputation
#'
#' This function performs MAGIC (Markov Affinity-based Graph Imputation of Cells)
#' imputation on a gene expression matrix.
#'
#' @param data_mat Expression matrix with genes as rows and cells as columns
#' @param t Time parameter for diffusion (default: 3)
#' @param k Number of nearest neighbors (default: 10)
#' @param verbose Print progress messages
#'
#' @return Imputed expression matrix
#'
#' @importFrom Matrix t
#' @importFrom FNN get.knn
RunMAGIC <- function(data_mat, t = 3, k = 10, verbose = TRUE) {
  require_packages("FNN")
  
  if (verbose) message("Computing distance matrix...")
  dist_mat <- FNN::get.knn(t(data_mat), k = k)$nn.dist
  
  if (verbose) message("Creating affinity matrix...")
  sigma <- median(dist_mat)
  W <- exp(-dist_mat^2 / (2 * sigma^2))
  
  if (verbose) message("Normalizing transition matrix...")
  P <- W / rowSums(W)
  
  if (verbose) message("Running diffusion...")
  P_t <- matrix_power(P, t)
  
  if (verbose) message("Imputing values...")
  imputed_mat <- t(P_t %*% t(data_mat))
  
  return(imputed_mat)
}

#' Run KNN-based smoothing imputation
#'
#' This function performs k-nearest neighbor based smoothing imputation
#' on a gene expression matrix.
#'
#' @param data_mat Expression matrix with genes as rows and cells as columns
#' @param k Number of nearest neighbors (default: 15)
#' @param verbose Print progress messages
#'
#' @return Imputed expression matrix
#'
#' @importFrom FNN get.knn
#' @importFrom Matrix t
RunKNNSmooth <- function(data_mat, k = 15, verbose = TRUE) {
  require_packages("FNN")
  
  # Find k nearest neighbors
  if (verbose) message("Finding ", k, " nearest neighbors...")
  knn_res <- FNN::get.knn(t(data_mat), k = k)
  
  # Perform smoothing
  if (verbose) message("Performing KNN smoothing...")
  n_cells <- ncol(data_mat)
  imputed_mat <- matrix(0, nrow = nrow(data_mat), ncol = n_cells)
  
  for (i in 1:n_cells) {
    neighbors <- knn_res$nn.index[i, ]
    imputed_mat[, i] <- rowMeans(data_mat[, neighbors, drop = FALSE])
  }
  
  return(imputed_mat)
}

#' Helper function to choose optimal k for ALRA
#'
#' @param d Singular values
#' @param q Quantile threshold
#' @return Optimal k value
choose_k <- function(d, q = 0.01) {
  ratios <- d[-length(d)] / d[-1]
  k <- which.max(ratios > quantile(ratios, 1 - q))[1]
  return(k)
}

#' Helper function for matrix power
#'
#' @param mat Input matrix
#' @param power Power to raise matrix to
#' @return Matrix raised to specified power
matrix_power <- function(mat, power) {
  if (power == 0) return(diag(nrow(mat)))
  if (power == 1) return(mat)
  if (power %% 2 == 0) {
    temp <- matrix_power(mat, power/2)
    return(temp %*% temp)
  }
  return(mat %*% matrix_power(mat, power - 1))
}
