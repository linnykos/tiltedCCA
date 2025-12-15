# Fine tune the common tilts, one for each latent dimension

Tune each latent dimension (representing a pair of canonical score
vectors) for an appropriate tilt of the common component

## Usage

``` r
fine_tuning(
  input_obj,
  max_iter = 5,
  fix_tilt_perc = NA,
  temp_path = NULL,
  tol = 1e-05,
  verbose = 0
)
```

## Arguments

- input_obj:

  `multiSVD` class, after creation via
  [`tiltedCCA()`](https://linnykos.github.io/tiltedCCA/reference/tiltedCCA.md)

- max_iter:

  maximum number of epochs (i.e., cycling through all latent dimensions)

- fix_tilt_perc:

  scalar between 0 and 1, or `NA`, for setting all the tilts in each
  latent dimension to the scalar if not `NA`

- temp_path:

  filepath for saving temporary progress to

- tol:

  small positive number

- verbose:

  non-negative integer

## Value

updated `multiSVD` object
