# Generate downsampled example datasets from local sources

library(SCP)

# 1) Visium HD (binned outputs) example
# Set your path here:
# hd_bin <- "/Users/you/Visium_HD/.../binned_outputs/square_008um"
# SaveExampleSpatial(hd_bin, out = "inst/extdata/visium_hd_subset.rds", n_spots = 5000)

# 2) Chromium 10x example
# h5 <- "/Users/you/Human_Kidney_4k/..._filtered_feature_bc_matrix.h5"
# SaveExampleChromium(h5, out = "inst/extdata/chromium_subset.rds", n_cells = 2000)

