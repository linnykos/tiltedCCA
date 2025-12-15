# Include SNN graphs to multiSVD

Computing and including the shared nearest-neighbor (SNN) graphs for
each modality and the target common manifold for an existing `multiSVD`
object

## Usage

``` r
compute_snns(
  input_obj,
  latent_k,
  num_neigh,
  bool_cosine = T,
  bool_intersect = T,
  min_deg = 1,
  tol = 1e-04,
  verbose = 0
)
```

## Arguments

- input_obj:

  `multiSVD` class, after creation via
  [`create_multiSVD()`](https://linnykos.github.io/tiltedCCA/reference/create_multiSVD.md)
  or
  [`form_metacells()`](https://linnykos.github.io/tiltedCCA/reference/form_metacells.md)

- latent_k:

  number of latent dimensions for the graph Laplacian bases

- num_neigh:

  number of neighbors for each cell, when constructing the SNN graphs

- bool_cosine:

  boolean, for using cosine distance if `T` or Euclidean distance if `F`

- bool_intersect:

  boolean, on whether or not to symmetrize (via the AND function) the
  SNN

- min_deg:

  minimum degree of each cell in the SNN

- tol:

  small positive number

- verbose:

  non-negative integer

## Value

updated `multiSVD` object
