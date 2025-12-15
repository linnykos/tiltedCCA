# Heatmap of the data

If `reserve_zero = T`, then reserve zero for white. If
`reserve_zero = F`, then all the greens are negative and all the reds
are positive

## Usage

``` r
plot_heatmat(
  dat,
  luminosity = F,
  asp = nrow(dat)/ncol(dat),
  reserve_zero = T,
  ...
)
```

## Arguments

- dat:

  matrix

- luminosity:

  boolean

- asp:

  numeric

- reserve_zero:

  boolean

- ...:

  additional graphical parameters

## Value

shows a plot but returns nothing
