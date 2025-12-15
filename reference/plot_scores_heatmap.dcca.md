# Side-by-side plot of the canonical scores as heatmaps

Side-by-side plot of the canonical scores as heatmaps

## Usage

``` r
plot_scores_heatmap.dcca(
  obj,
  main_vec = c("Common score", "Distinct score 1", "Distinct score 2"),
  membership_vec = NA,
  num_col = 10,
  log_scale = F,
  scaling_power = 1,
  luminosity = F
)
```

## Arguments

- obj:

  output of either `generate_data` or `dcca_decomposition`

- main_vec:

  vector of characters for the title of the plots

- membership_vec:

  factor vector

- num_col:

  positive integers for number of distinct colors

- log_scale:

  boolean

- scaling_power:

  positive numeric

- luminosity:

  boolean

## Value

shows a plot but returns nothing
