#' Create example Chromium scRNA dataset (downsampled)
#'
#' @param h5_path 10x filtered feature-barcode matrix in HDF5 format.
#' @param out Path to write (.rds or .h5seurat).
#' @param n_cells Number of cells to sample (default 2000).
#' @return Invisible path to the written file.
#' @export
SaveExampleChromium <- function(h5_path, out, n_cells = 2000) {
  if (!file.exists(h5_path)) stop("h5_path not found: ", h5_path)
  mat <- Seurat::Read10X_h5(h5_path)
  srt <- Seurat::CreateSeuratObject(counts = mat)
  if (ncol(srt) > n_cells) {
    set.seed(1)
    srt <- subset(srt, cells = sample(colnames(srt), n_cells))
  }
  if (grepl("\\.rds$", out, ignore.case = TRUE)) {
    saveRDS(srt, file = out)
  } else if (grepl("\\.h5seurat$", out, ignore.case = TRUE)) {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("SeuratDisk is required to write .h5seurat. Install SeuratDisk.")
    }
    SeuratDisk::SaveH5Seurat(srt, filename = out, overwrite = TRUE)
  } else {
    stop("Unsupported output extension. Use .rds or .h5seurat")
  }
  invisible(out)
}

