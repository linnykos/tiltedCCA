# Plot cell enrichment

Plot cell enrichment

## Usage

``` r
plot_cell_enrichment(
  cell_enrichment_res,
  col_palette,
  cex_axis = 1,
  cex_lab = 1,
  lwd_diag = 3,
  lwd_grid = 2,
  mar = c(4, 5, 0.5, 0.5),
  xlab_1 = "Modality 1 Distinct",
  xlab_2 = "Modality 2 Distinct",
  ylab = "Common"
)
```

## Arguments

- cell_enrichment_res:

  the result of
  [`tiltedCCA::postprocess_cell_enrichment`](https://linnykos.github.io/tiltedCCA/reference/postprocess_cell_enrichment.md)

- col_palette:

  vector of colors (as strings), equal to the number of celltypes in
  `cell_enrichment_res`

- cex_axis:

  positive number (graphical parameter)

- cex_lab:

  positive number (graphical parameter)

- lwd_diag:

  positive number (graphical parameter)

- lwd_grid:

  positive number (graphical parameter)

- mar:

  vector of 4 positive numbers (graphical parameter)

- xlab_1:

  character (graphical parameter)

- xlab_2:

  character (graphical parameter)

- ylab:

  character (graphical parameter)

## Value

makes a plot but does not return anything
