# Main workhorse of dcca_factor

Given the two matrices (given by `svd_1` and `svd_2`) and the CCA
solution in `cca_res`, compute the common scores. This calls the
functions `.common_decomposition` and `.compute_distinct_score`.

## Usage

``` r
.tiltedCCA_common_score(
  averaging_mat,
  cca_res,
  discretization_gridsize,
  enforce_boundary,
  fix_tilt_perc,
  snn_bool_cosine,
  snn_bool_intersect,
  snn_k,
  snn_min_deg,
  snn_num_neigh,
  svd_1,
  svd_2,
  target_dimred,
  verbose = 0
)
```

## Arguments

- averaging_mat:

  sparse matrix

- cca_res:

  returned object from `.cca`

- discretization_gridsize:

  positive integer for how many values between 0 and 1 (inclusive) to
  search the appropriate amount of tilt over

- enforce_boundary:

  boolean, on whether or not the tilt is required to stay between the
  two canonical score vectors

- fix_tilt_perc:

  boolean or a numeric. If `FALSE`, then the tilt is adaptively
  determined, and if `TRUE`, then the tilt is set to be equal to `0.5`.
  If numeric, the value should be between `0` and `1`, which the tilt
  will be set to.

- snn_bool_cosine:

  boolean

- snn_bool_intersect:

  boolean

- snn_k:

  integer

- snn_min_deg:

  integer

- snn_num_neigh:

  integer

- svd_1:

  SVD of the denoised variant of `mat_1` from `dcca_factor`

- svd_2:

  SVD of the denoised variant of `mat_2` from `dcca_factor`

- target_dimred:

  matrix

- verbose:

  non-negative integer

## Value

list
