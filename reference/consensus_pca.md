# Computing Consensus PCA

Computing Consensus PCA

## Usage

``` r
consensus_pca(
  mat_1,
  mat_2,
  dims_1,
  dims_2,
  dims_consensus,
  apply_pca = T,
  center_1 = T,
  center_2 = T,
  center_consensus = T,
  normalize_row = F,
  normalize_singular_value = T,
  recenter_1 = F,
  recenter_2 = F,
  recenter_consensus = F,
  rescale_1 = F,
  rescale_2 = F,
  rescale_consensus = F,
  scale_1 = T,
  scale_2 = T,
  scale_consensus = T,
  scale_max_1 = NULL,
  scale_max_2 = NULL,
  scale_max_consensus = NULL,
  svd_1 = NULL,
  svd_2 = NULL,
  tol = 0.001,
  verbose = 0
)
```

## Arguments

- mat_1:

  data matrix of `n` cells and `p1` features for Modality 1

- mat_2:

  data matrix of `n` cells and `p2` features for Modality 2

- dims_1:

  vector of latent dimensions for `mat_1` used for analysis

- dims_2:

  vector of latent dimensions for `mat_2` used for analysis

- dims_consensus:

  vector of latent dimensions for the Consensus PCA

- apply_pca:

  boolean, where PCA is combined set of latent dimensions (from both
  modalities) if `T` or not if `F`

- center_1:

  boolean, to center each the feature in Modality 1 prior to computing
  latent dimensions

- center_2:

  boolean, to center each the feature in Modality 2 prior to computing
  latent dimensions

- center_consensus:

  boolean, to center the combined set of latent dimensions (from both
  modalities)

- normalize_row:

  boolean, to normalize each cell's latent vector after
  dimension-reduction for both modalities

- normalize_singular_value:

  boolean, to normalize each modality by its largest singular value

- recenter_1:

  boolean, to center each latent dimension in Modality 1 after computing
  latent dimensions

- recenter_2:

  boolean, to center each latent dimension in Modality 2 after computing
  latent dimensions

- recenter_consensus:

  boolean, to center the latent dimension after Consensus PCA

- rescale_1:

  boolean, to rescale each latent dimension in Modality 1 after
  computing latent dimensions

- rescale_2:

  boolean, to rescale each latent dimension in Modality 2 after
  computing latent dimensions

- rescale_consensus:

  boolean, to rescale the latent dimension after Consensus PCA

- scale_1:

  boolean, to rescale each the feature in Modality 1 prior to computing
  latent dimensions

- scale_2:

  boolean, to rescale each the feature in Modality 2 prior to computing
  latent dimensions

- scale_consensus:

  boolean, to rescale the combined set of latent dimensions (from both
  modalities)

- scale_max_1:

  numeric or `NULL`, to threshold Modality 1 in magnitude prior to
  computing latent dimensions

- scale_max_2:

  numeric or `NULL`, to threshold Modality 2 in magnitude prior to
  computing latent dimensions

- scale_max_consensus:

  numeric or `NULL`, to threshold the combined set of latent dimensions
  in magnitude prior to computing latent dimensions

- svd_1:

  list of `u`, `d`, `v` for the SVD of Modality 1 if it's already
  computed

- svd_2:

  list of `u`, `d`, `v` for the SVD of Modality 2 if it's already
  computed

- tol:

  small positive number

- verbose:

  non-negative integer

## Value

`consensusPCA` object
