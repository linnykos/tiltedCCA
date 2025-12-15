# Computing cell enrichments for multiSVD object

Computing cell enrichments for multiSVD object

## Usage

``` r
# S3 method for class 'multiSVD'
postprocess_cell_enrichment(
  input_obj,
  membership_vec,
  max_subsample = min(1000, length(membership_vec)),
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  `multiSVD` object, after using
  [`tiltedCCA_decomposition()`](https://linnykos.github.io/tiltedCCA/reference/tiltedCCA_decomposition.md)

- membership_vec:

  factor vector of the same length as the number of cells in `multiSVD`,
  denoting the cell-types for each cell

- max_subsample:

  maximum of cells to sample for each cell-type. If there are more than
  `max_subsample` cells in a cell-type (dictated by `membership_vec`), a
  random subset of cells will be selected for the sake of this function

- verbose:

  non-negative integer

- ...:

  additional arguments
