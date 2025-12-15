# Find the marker genes the separate one cell type from all other cell types

Based on the output of
[`tiltedCCA::differential_expression`](https://linnykos.github.io/tiltedCCA/reference/differential_expression.md),
for a particular cell type, 1) when comparing against every other cell
type, find all the variables that have an adjusted p-value smaller than
`p_val_thres` and also have a log fold change larger than
`log_thres_name`, and 2) if the variable was selected for more than
`num_quantile_threshold` percentage of the comparisons with other cell
types, we deem it as a marker gene for this cell type.

## Usage

``` r
postprocess_marker_variables(
  de_list,
  log_thres_name = "avg_log2FC",
  log_thres = 2,
  num_quantile_threshold = 0.25,
  p_val_name = "p_val",
  p_val_adj_name = "p_val_adj",
  p_val_thres = 1e-04
)
```

## Arguments

- de_list:

  output of the
  [`tiltedCCA::differential_expression`](https://linnykos.github.io/tiltedCCA/reference/differential_expression.md)
  function

- log_thres_name:

  character of the column in each of `de_list$de_list` for the log fold
  change

- log_thres:

  positive value, for selecting variables with log fold change larger
  than `log_thres`

- num_quantile_threshold:

  value between 0 and 1

- p_val_name:

  character of the column in each of `de_list$de_list` for the p-value

- p_val_adj_name:

  character of the column in each of `de_list$de_list` for the
  multiple-testing adjusted p-value

- p_val_thres:

  small positive value between 0 and 1, for selecting variables with
  p-values smaller than `p_val_thres`

## Value

a list of variable vectors, one vector for each element in
`de_list$level_vec`
