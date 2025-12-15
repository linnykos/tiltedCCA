# Using the CCA solution, compute the score matrices.

This is called by `.tiltedCCA_common_score`. It's called unnormalized
scores since if there are `n` rows (i.e., `nrow(svd_1$u)`), then the
return matrices will be orthogonal matrices where the matrix multiplied
by itself is a diagonal matrix with all values `n`.

## Usage

``` r
.compute_unnormalized_scores(svd_1, svd_2, cca_res)
```

## Arguments

- svd_1:

  SVD of the denoised variant of `mat_1` from `dcca_factor`

- svd_2:

  SVD of the denoised variant of `mat_2` from `dcca_factor`

- cca_res:

  returned object from `.cca`

## Value

list of two matrices
