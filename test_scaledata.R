# Minimal test for ScaleData with split layers
devtools::load_all()
data("panc8_sub")

cat("=== Testing ScaleData with split layers ===\n")
cat("Initial object:\n")
print(panc8_sub)

# Split layers
cat("\nSplitting layers by tech...\n")
panc8_sub[["RNA"]] <- split(panc8_sub[["RNA"]], f = panc8_sub$tech)

cat("\nAfter split:\n")
print(panc8_sub)

# Normalize
cat("\nNormalizing...\n")
panc8_sub <- NormalizeData(panc8_sub)

# Find variable features
cat("\nFinding variable features...\n")
panc8_sub <- FindVariableFeatures(panc8_sub, nfeatures = 2000)
hvf <- VariableFeatures(panc8_sub)
cat("Found", length(hvf), "variable features\n")

# Try ScaleData with variable features only
cat("\nTrying ScaleData with features parameter...\n")
system.time({
  panc8_sub <- ScaleData(panc8_sub, features = hvf, verbose = TRUE)
})

cat("\n=== SUCCESS ===\n")
print(panc8_sub)
