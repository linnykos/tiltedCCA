# Compute how aligned the common component for a feature is with that feature in the original data modality

The common matrix for modality `input_assay` is extracted, and then
compared against another matrix (depending on `bool_use_denoised` is
set). (For example, if `bool_use_denoised=TRUE`, then this is the common
plus distinct matrix.) Then, a regression is performed, one per feature
(i.e., gene or protein) that regresses latter matrix onto the common
matrix, and the R-squared (one per feature) is returned.

## Usage

``` r
postprocess_modality_alignment(
  input_obj,
  bool_use_denoised,
  input_assay,
  bool_center = T,
  bool_scale = T,
  bool_regression_include_intercept = T,
  min_subsample_cell = NULL,
  seurat_celltype_variable = "celltype",
  seurat_obj = NULL,
  seurat_assay = NULL,
  seurat_slot = "data",
  verbose = 1
)
```

## Arguments

- input_obj:

  a `multiSVD_obj` that was the output of
  [`tiltedCCA::tiltedCCA_decomposition`](https://linnykos.github.io/tiltedCCA/reference/tiltedCCA_decomposition.md)

- bool_use_denoised:

  boolean. If `TRUE`, then the common component is compared against the
  common plus distinct component. If `FALSE`, then the common component
  is compared against the original data matrix in slot `seurat_slot` in
  `seurat_obj[[seurat_assay]]`

- input_assay:

  integer of `1` or `2`, denoting which modality is being analyzed

- bool_center:

  boolean if all the features in the common component are centered prior
  to the comparison

- bool_scale:

  boolean if all the features in the common component are rescaled prior
  to the comparison

- bool_regression_include_intercept:

  boolean if the regression analysis

- min_subsample_cell:

  if not `NULL`, subsample `min_subsample_cell` number of cells of each
  cell type (denoted by in `seurat_obj$seurat_celltype_variable`)

- seurat_celltype_variable:

  a character where `seurat_obj$seurat_celltype_variable` denotes the
  cell type for each cell in `seurat_obj`

- seurat_obj:

  the `Seurat` object that was used to compute `input_obj`, the
  `multiSVD_obj`

- seurat_assay:

  the assay to extract the data matrix, which is relevant
  `bool_use_denoised=FALSE`

- seurat_slot:

  the slot to extract the data matrix, which is relevant
  `bool_use_denoised=FALSE`

- verbose:

  non-negative integer

## Value

a vector of R-squared values for each variable
