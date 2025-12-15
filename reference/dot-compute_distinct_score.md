# Compute the distinct scores

Given `score_1` and `score_2`, and having already computed
`common_score`, compute the distinct scores. This is more-or-less a
simple subtraction, but we need to handle situations where we might need
to "fill-in extra dimensions"

## Usage

``` r
.compute_distinct_score(score_1, score_2, common_score)
```

## Arguments

- score_1:

  matrix

- score_2:

  matrix

- common_score:

  matrix

## Value

list of two matrices
