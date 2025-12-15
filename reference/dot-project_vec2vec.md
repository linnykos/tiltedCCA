# Projection of vector onto another vector

Returns the component of `vec1` that is orthogonal to `vec2` if
`orthogonal` is `TRUE`, and the component of `vec1` that is parallel to
`vec2` if `orthogonal` is `FALSE`.

## Usage

``` r
.project_vec2vec(vec1, vec2, orthogonal, tol = 1e-06)
```

## Arguments

- vec1:

  vector

- vec2:

  vector

- orthogonal:

  boolean

- tol:

  small positive number

## Value

vector
