# Tilted-CCA Decomposition

After computing the appropriate tilts, recover the full decomposition of
both modalities.

## Usage

``` r
tiltedCCA_decomposition(
  input_obj,
  bool_modality_1_full = T,
  bool_modality_2_full = T,
  verbose = 0
)
```

## Arguments

- input_obj:

  `multiSVD` class, after creation via
  [`tiltedCCA()`](https://linnykos.github.io/tiltedCCA/reference/tiltedCCA.md)
  or
  [`fine_tuning()`](https://linnykos.github.io/tiltedCCA/reference/fine_tuning.md)

- bool_modality_1_full:

  boolean, to compute the full-matrix decomposition for Modality 1 if
  `T` or only the low-dimensional representation if `F`

- bool_modality_2_full:

  boolean, to compute the full-matrix decomposition for Modality 2 if
  `T` or only the low-dimensional representation if `F`

- verbose:

  non-negative integer

## Value

updated `multiSVD` object
