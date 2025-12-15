# Side-by-side plot of the canonical scores, colored by membership

Side-by-side plot of the canonical scores, colored by membership

## Usage

``` r
plot_scores(
  obj,
  membership_vec,
  col_vec = (scales::hue_pal())(length(levels(membership_vec))),
  xlim = NA,
  decomposition = F
)
```

## Arguments

- obj:

  output of either `generate_data` or `dcca_decomposition`

- membership_vec:

  factor vector

- col_vec:

  vector of colors

- xlim:

  custom `xlim` graphical argument

- decomposition:

  boolean

## Value

shows a plot but returns nothing
