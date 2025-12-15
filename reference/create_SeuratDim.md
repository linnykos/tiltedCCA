# Create Seurat dimension reduction objects of Tilted-CCA

Create Seurat dimension reduction objects of Tilted-CCA

## Usage

``` r
create_SeuratDim(
  input_obj,
  what,
  aligned_umap_assay = NULL,
  scale_max_1 = NULL,
  scale_max_2 = NULL,
  seurat_obj = NULL,
  seurat_assay = "RNA",
  suppress_warnings = TRUE,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  a `multiSVD` object that is the output of `tiltedCCA_decomposition`

- what:

  a character, either `"common"`, `"distinct_1"`, or `"distinct_2"` for
  which embedding to construct the visualization for

- aligned_umap_assay:

  either `NULL` or a UMAP assay in `seurat_assay` (in which case the
  resulting UMAP will be rotated to best mimic the relative orientation
  of cells in `seurat_obj[[aligned_umap_assay]]`)

- scale_max_1:

  numeric or `NULL`, to threshold Modality 1 in magnitude prior to
  computing latent dimensions

- scale_max_2:

  numeric or `NULL`, to threshold Modality 2 in magnitude prior to
  computing latent dimensions

- seurat_obj:

  the `Seurat` object that was used to compute `input_obj`, the
  `multiSVD_obj`

- seurat_assay:

  the assay in `seurat_obj` to assign the resulting embedding to

- suppress_warnings:

  boolean to suppress the warning when running
  [`Seurat::RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html)

- verbose:

  non-negative integer

- ...:

  extra arguments to
  [`Seurat::RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html)

## Value

a `DimReduc` `SeuratObject`
