# Computing cell enrichments for matrix object

Computing cell enrichments for matrix object

## Usage

``` r
# Default S3 method
postprocess_cell_enrichment(
  input_obj,
  membership_vec,
  num_neigh,
  bool_cosine = T,
  bool_intersect = T,
  max_subsample = min(1000, length(membership_vec)),
  min_deg = 1,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  matrix of `n` cells and `p` variable

- membership_vec:

  factor vector of the same length as the number of cells in `multiSVD`,
  denoting the cell-types for each cell

- num_neigh:

  number of nearest neighbors

- bool_cosine:

  boolean, for using cosine distance if `T` or Euclidean distance if `F`

- bool_intersect:

  boolean, on whether or not to symmetrize (via the AND function) the
  SNN

- max_subsample:

  maximum of cells to sample for each cell-type. If there are more than
  `max_subsample` cells in a cell-type (dictated by `membership_vec`), a
  random subset of cells will be selected for the sake of this function

- min_deg:

  minimum degree of each cell in the SNN

- verbose:

  non-negative integer

- ...:

  additional arguments
