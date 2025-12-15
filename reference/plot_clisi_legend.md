# Plot the cLISI legend

Plot the cLISI legend

## Usage

``` r
plot_clisi_legend(
  clisi_obj,
  col_vec = (scales::hue_pal())(nrow(clisi_obj$common_clisi$membership_info)),
  percent_coverage = 1,
  pch = 16,
  cex_point = 1,
  cex_text = 1,
  text_nudge = 0,
  xlim = c(0, 1),
  ...
)
```

## Arguments

- clisi_obj:

  output of `clisi_information`

- col_vec:

  vector of colors

- percent_coverage:

  numeric

- pch:

  `pch` parameter

- cex_point:

  `cex` parameter for points

- cex_text:

  `cex` parameter for the text

- text_nudge:

  x-axis offset for the text

- xlim:

  `xlim` parameter for the plot

- ...:

  additional graphical parameters

## Value

nothing
