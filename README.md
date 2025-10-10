# SCP - Seurat V5 Compatible

An adaptation of the [SCP package](https://github.com/zhanghao-njmu/SCP) for Seurat V5 support.

## Key Features

### Seurat V5 Compatibility
- **V5 Support**: Full compatibility with Seurat V5's layer-based data structure
- **Version Agnostic**: Seamlessly handles both V4 and V5 objects
- **Automatic Conversion**: Detection and conversion of legacy objects

### Core Functions
- `get_seurat_data()` / `set_seurat_data()` - Unified data access across versions
- `IsSeurat5()` / `EnsureSeurat5()` - Version detection and conversion
- `get_feature_metadata()` / `set_feature_metadata()` - Feature metadata management
- `get_var_features()` - Variable feature extraction

## Installation

```r
# Install from GitHub
devtools::install_github("mianaz/SCP-SeuratV5")
```

## Requirements

### Core Dependencies
- R >= 4.1.0
- Seurat >= 5.0.0
- SeuratObject >= 5.0.0

### Python Environment (Optional)
For Python-based features (PAGA, scVelo, etc.):
```r
# Set up Python environment
SCP::PrepareEnv()
```

## Quick Start

```r
library(SCP)

# Load your data (automatic V5 conversion if needed)
srt <- LoadSeuratObject("your_data.rds")

# Run standard workflow
srt <- RunSCP(srt,
  species = "Homo_sapiens",
  db = "GO_BP",
  reduction = "umap"
)

# Visualize results
CellDimPlot(srt, group.by = "cell_type")
```

## Documentation

For detailed documentation and tutorials, see the original [SCP repository](https://github.com/zhanghao-njmu/SCP).

## License

GPL (>= 3)
