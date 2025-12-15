# Helper function with the CCA function

Called by the `.cca` function. Recall when computing CCA, the main
matrix we need compute is, roughly speaking, half-inverse of the first
covariance times the cross-covariance matrix times the half-inverse of
the second covariance. If we had the SVD of the two original matrices,
this is actually equivalent to the product of the left singular vectors.

## Usage

``` r
.compute_cca_aggregate_matrix(svd_1, svd_2, augment)
```

## Arguments

- svd_1:

  SVD of the denoised variant of `mat_1` from `dcca_factor`

- svd_2:

  SVD of the denoised variant of `mat_2` from `dcca_factor`

- augment:

  boolean. If `TRUE`, augment the matrix with either rows or columns
  with 0's so the dimension of the output matrix matches those in
  `svd_1` and `svd_2`

## Value

matrix
