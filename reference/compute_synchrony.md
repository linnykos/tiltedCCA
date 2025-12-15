# Compute the synchrony based on Tilted-CCA

Compute the synchrony based on Tilted-CCA

## Usage

``` r
compute_synchrony(input_obj, anchor_modality = 1)
```

## Arguments

- input_obj:

  `multiSVD` class, after applying
  [`tiltedCCA_decomposition()`](https://linnykos.github.io/tiltedCCA/reference/tiltedCCA_decomposition.md)

- anchor_modality:

  numeric of `1` or `2`, denoting which modality the synchrony is
  computed on

## Value

a matrix with `n` cells and 2 columns, which one column (named
`synchrony`) denotes the synchrony score (between 0 and 1, where a
number closer to 1 means the cell is synchronous between both
modalities) and (named `synchrony_rescaled`), which is simply a monotone
transformation of the first column for easier visualization
