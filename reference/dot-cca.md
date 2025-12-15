# Perform CCA

This function takes either two matrices or two SVDs. Both `input_1` and
`input_2` must be "of the same type." Calls the
`.compute_cca_aggregate_matrix` function.

## Usage

``` r
.cca(input_1, input_2, dims_1, dims_2, return_scores, tol = 1e-06)
```

## Arguments

- input_1:

  first input

- input_2:

  second input

- dims_1:

  desired latent dimensions of data matrix 1. Only used if `input_1` is
  a matrix, not if it's a list representing the SVD

- dims_2:

  desired latent dimensions of data matrix 2. Only used if `input_2` is
  a matrix, not if it's a list representing the SVD

- return_scores:

  boolean. If `TRUE`, return the scores (i.e., matrices where the rows
  are the cells). If `FALSE`, return the loadings (i.e., matrices where
  the rows are the variables). Either way, one of the output matrices
  will have `rank_1` columns and another will have `rank_2` columns

- tol:

  small numeric

## Value

list
