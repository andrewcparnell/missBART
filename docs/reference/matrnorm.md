# Sampling from a matrix normal distribution MV(mu, U, V)

Sampling from a matrix normal distribution MV(mu, U, V)

## Usage

``` r
matrnorm(mu, U, V)
```

## Arguments

- mu:

  mean matrix (real nxp matrix)

- U:

  scale matrix (positive-definite real nxn matrix)

- V:

  scale matrix (positive-definite real pxp matrix)

## Value

one sample from a MV(mu, U, V) distribution

## Examples

``` r
matrnorm(matrix(0, nrow=5, ncol=2), diag(5), diag(2))
#>              [,1]       [,2]
#> [1,] -1.400043517  1.1484116
#> [2,]  0.255317055 -1.8218177
#> [3,] -2.437263611 -0.2473253
#> [4,] -0.005571287 -0.2441996
#> [5,]  0.621552721 -0.2827054
```
