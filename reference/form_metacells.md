# Include meta-cell information to multiSVD

Include the hard clustering information and compute the meta-cells for
an existing `multiSVD` object.

## Usage

``` r
form_metacells(
  input_obj,
  large_clustering_1,
  large_clustering_2,
  num_metacells = NULL,
  min_size = 5,
  verbose = 0
)
```

## Arguments

- input_obj:

  `multiSVD` class, after creation via
  [`create_multiSVD()`](https://linnykos.github.io/tiltedCCA/reference/create_multiSVD.md)

- large_clustering_1:

  factor among the `n` cells or `NULL`, representing the hard clustering
  structure of Modality 1

- large_clustering_2:

  factor among the `n` cells or `NULL`, representing the hard clustering
  structure of Modality 2

- num_metacells:

  number of desired meta-cells

- min_size:

  smallest number of cells for a meta-cell

- verbose:

  non-negative integer

## Value

updated `multiSVD` object
