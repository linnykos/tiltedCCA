# Tilted-CCA Factorization

Computing the common-distinct decomposition via CCA, by setting each
latent dimension to the same common tilt.

## Usage

``` r
tiltedCCA(
  input_obj,
  discretization_gridsize = 21,
  enforce_boundary = F,
  fix_tilt_perc = F,
  verbose = 0
)
```

## Arguments

- input_obj:

  `multiSVD` class, after creation via
  [`compute_snns()`](https://linnykos.github.io/tiltedCCA/reference/compute_snns.md)

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

- verbose:

  non-negative integer

  For the tilt values (possibly set in `fix_tilt_perc`), values close to
  0 or 1 means the common space resembles the canonical scores of
  `mat_2` or `mat_1` respectively.

## Value

updated `multiSVD` object
