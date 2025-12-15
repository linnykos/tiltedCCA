# Compute the negative log-10 p-values across all cell types

This function does the following: There are varible cell types stored in
`de_list$level_vec`. For each cell type, enumerate the p-value (across
different features, i.e., genes or proteins) when comparing that cell
type against all other cell types. Return the p-value at the quantile of
`quantile_value` among all these comparisons for this particular cell
type. Then, repeat for all cell types, resulting in a p-value for each
cell type and each feature. Lastly, for each feature, take the maximum
of the -log10 p-value (across all the cell types).

## Usage

``` r
postprocess_depvalue(
  de_list,
  exclude_celltypes = NULL,
  maximum_output = NULL,
  quantile_value = 0.75,
  verbose = 1
)
```

## Arguments

- de_list:

  output of the
  [`tiltedCCA::differential_expression`](https://linnykos.github.io/tiltedCCA/reference/differential_expression.md)
  function

- exclude_celltypes:

  a vector of strings containing cell types to be excluded in this
  function's calculations. All the cell types should be in
  `de_list$level_vec`

- maximum_output:

  any output (i.e., negative log10 pvalue) larger than `maximum_output`
  are set to be `maximum_output`

- quantile_value:

  for a particular variable, when comparing one cell type to all other
  cell types, output the p-value at the quantile of `quantile_value` .
  (prior to taking the maximum of the -log10 across all cell types)

- verbose:

  non-negative integer

## Value

a vector of the -log10 p-value for each variable in `de_list`
