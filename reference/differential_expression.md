# Compute the differential expression across each pair of cell types

Compute the differential expression across each pair of cell types

## Usage

``` r
differential_expression(
  seurat_obj,
  assay,
  idents,
  slot = "data",
  test_use = "wilcox",
  verbose = T
)
```

## Arguments

- seurat_obj:

  seurat object

- assay:

  name of assay in the `seurat_obj`

- idents:

  variables name inside `seuart_obj` that contains the cell-type
  information (or clustering information) for each cell

- slot:

  slot of `seurat_obj[[assay]]` that informs with data matrix will be
  used in the DE test

- test_use:

  either `"MAST"` or `"wilcox"` dictates which DE test will be used

- verbose:

  non-negative integer

## Value

a list with elements `"combn_mat"`, `"de_list"` and `"level_vec"`
