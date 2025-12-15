# Create the object to initialize Tilted-CCA

Initialize the `multiSVD` object by computing each modality's SVD.

## Usage

``` r
create_multiSVD(
  mat_1,
  mat_2,
  dims_1,
  dims_2,
  center_1 = T,
  center_2 = T,
  normalize_row = T,
  normalize_singular_value = T,
  recenter_1 = F,
  recenter_2 = F,
  rescale_1 = F,
  rescale_2 = F,
  scale_1 = T,
  scale_2 = T,
  scale_max_1 = NULL,
  scale_max_2 = NULL,
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

- center_1:

  boolean, to center each the feature in Modality 1 prior to computing
  latent dimensions

- center_2:

  boolean, to center each the feature in Modality 2 prior to computing
  latent dimensions

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

- rescale_1:

  boolean, to rescale each latent dimension in Modality 1 after
  computing latent dimensions

- rescale_2:

  boolean, to rescale each latent dimension in Modality 2 after
  computing latent dimensions

- scale_1:

  boolean, to rescale each the feature in Modality 1 prior to computing
  latent dimensions

- scale_2:

  boolean, to rescale each the feature in Modality 2 prior to computing
  latent dimensions

- scale_max_1:

  numeric or `NULL`, to threshold Modality 1 in magnitude prior to
  computing latent dimensions

- scale_max_2:

  numeric or `NULL`, to threshold Modality 2 in magnitude prior to
  computing latent dimensions

- verbose:

  non-negative integer

## Value

`multiSVD` object
