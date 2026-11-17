# Sampling from a multivariate normal distribution with precision Q and mean equal to the solution of Q mu = b

Sampling from a multivariate normal distribution with precision Q and
mean equal to the solution of Q mu = b

## Usage

``` r
rMVNc(b, Q, is_chol = FALSE)
```

## Arguments

- b:

  a vector

- Q:

  precision matrix

- is_chol:

  logical

## Value

One sample from that distribution

## Examples

``` r
rMVNc(c(0,0), diag(2))
#> [1] 0.2121665 0.7099192
```
