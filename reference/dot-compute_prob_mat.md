# Compute probability matrix

Compute probability matrix

## Usage

``` r
.compute_prob_mat(B_mat, membership_vec)
```

## Arguments

- B_mat:

  symmetric connectivity matrix

- membership_vec:

  vector containing values `1` through `ncol(B_mat)`

## Value

symmetric matrix of dimension `length(membership_vec)`
