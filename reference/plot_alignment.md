# Plot alignment of variables

Plot alignment of variables

## Usage

``` r
plot_alignment(
  rsquare_vec,
  logpval_vec,
  bool_hide_points = F,
  bool_mark_ymedian = F,
  bool_polygon_mean = T,
  bool_truncate_xaxis = T,
  bool_white_bg = T,
  bty = "n",
  cex_axis = 1,
  cex_gene_highlight_inner = 3,
  cex_gene_highlight_outer = 4,
  cex_lab = 1,
  cex_points = 2,
  col_points = grDevices::rgb(0.5, 0.5, 0.5, 0.1),
  col_gene_highlight = 2,
  col_gene_highlight_border = "white",
  col_grid = "gray",
  density = 20,
  gene_names = NULL,
  lty_grid_major = 2,
  lty_grid_minor = 3,
  lty_polygon = 2,
  lwd_axis = 1,
  lwd_axis_ticks = 1,
  lwd_grid_major = 1,
  lwd_grid_minor = 0.5,
  lwd_polygon = 1,
  lwd_polygon_bold = 2,
  mark_median_xthres = 1,
  pch = 16,
  xaxt_num_ticks = 10,
  xaxt_grid_spacing = 2,
  xlab = "Separability (Median -Log10(p-value))",
  xlim = NULL,
  yaxt_by = 0.1,
  yaxt_num_ticks = 2,
  ylab = "Alignment with common space (R^2)",
  ylim = NULL,
  verbose = T,
  ...
)
```

## Arguments

- rsquare_vec:

  output of
  [`tiltedCCA::postprocess_modality_alignment`](https://linnykos.github.io/tiltedCCA/reference/postprocess_modality_alignment.md)

- logpval_vec:

  output of
  [`tiltedCCA::postprocess_depvalue`](https://linnykos.github.io/tiltedCCA/reference/postprocess_depvalue.md)

- bool_hide_points:

  boolean (graphical parameter)

- bool_mark_ymedian:

  boolean (graphical parameter)

- bool_polygon_mean:

  boolean (graphical parameter)

- bool_truncate_xaxis:

  boolean (graphical parameter)

- bool_white_bg:

  boolean (graphical parameter)

- bty:

  character (graphical parameter)

- cex_axis:

  positive number (graphical parameter)

- cex_gene_highlight_inner:

  positive number (graphical parameter)

- cex_gene_highlight_outer:

  positive number (graphical parameter)

- cex_lab:

  positive number (graphical parameter)

- cex_points:

  positive number (graphical parameter)

- col_points:

  character for a color (graphical parameter)

- col_gene_highlight:

  character for a color (graphical parameter)

- col_gene_highlight_border:

  character for a color (graphical parameter)

- col_grid:

  character for a color (graphical parameter)

- density:

  positive number, `NA`, or `NULL` (graphical parameter)

- gene_names:

  vector of charaters for which genes (in `names(rsquare_vec)`) to
  highlight

- lty_grid_major:

  positive number (graphical parameter)

- lty_grid_minor:

  positive number (graphical parameter)

- lty_polygon:

  positive number (graphical parameter)

- lwd_axis:

  positive number (graphical parameter)

- lwd_axis_ticks:

  positive number (graphical parameter)

- lwd_grid_major:

  positive number (graphical parameter)

- lwd_grid_minor:

  positive number (graphical parameter)

- lwd_polygon:

  positive number (graphical parameter)

- lwd_polygon_bold:

  positive number (graphical parameter)

- mark_median_xthres:

  boolean (graphical parameter)

- pch:

  positive integer (graphical parameter)

- xaxt_num_ticks:

  positive integer (graphical parameter)

- xaxt_grid_spacing:

  positive integer (graphical parameter)

- xlab:

  character (graphical parameter)

- xlim:

  vector of 2 numerics, or `NULL` (graphical parameter)

- yaxt_by:

  positive number (graphical parameter)

- yaxt_num_ticks:

  positive integer (graphical parameter)

- ylab:

  character (graphical parameter)

- ylim:

  vector of 2 numerics, or `NULL` (graphical parameter)

- verbose:

  non-negative integer

- ...:

  extra arguments to
  [`graphics::plot`](https://rdrr.io/r/graphics/plot.default.html)

## Value

makes a plot but does not return anything
