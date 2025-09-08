# SCP - Seurat V5 Compatible

This is a minimal adaptation of the [SCP package](https://github.com/zhanghao-njmu/SCP) to support Seurat V5 compatibility.

## Key Changes for Seurat V5

### Core Compatibility Functions
- `get_seurat_data()` - Version-agnostic data access (LayerData for V5, GetAssayData for V4)  
- `set_seurat_data()` - Version-agnostic data writing
- `IsSeurat5()` - Version detection utility
- `get_feature_metadata()` / `set_feature_metadata()` - Feature metadata handling
- `get_var_features()` - Variable features extraction

### Parameter Updates
- Updated `slot` → `layer` parameter throughout workflow functions
- Enhanced `check_srtList()` and `check_DataType()` for V5 object validation
- Automatic V4→V5 object conversion when needed

## Installation

```r
# Install from GitHub
devtools::install_github("your-username/SCP-SeuratV5")
```

## Dependencies

- Seurat >= 5.0.0
- SeuratObject >= 5.0.0
- All original SCP dependencies

## Usage

The package maintains the same API as the original SCP package, with automatic V5/V4 compatibility handling behind the scenes.

For complete documentation and examples, see the original [SCP repository](https://github.com/zhanghao-njmu/SCP).