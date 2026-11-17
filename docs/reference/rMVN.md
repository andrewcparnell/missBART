# Sampling once from a multivariate normal distribution with mean mu and precision Q

Sampling once from a multivariate normal distribution with mean mu and
precision Q

## Usage

``` r
rMVN(mu, Q)
```

## Arguments

- mu:

  mean vector

- Q:

  precision matrix

## Value

One sample from a multivariate normal distribution with mean mu and
precision Q

## Examples

``` r
rMVN(c(0,0), diag(2))
#> [1] -0.2457808  0.4639714
```
