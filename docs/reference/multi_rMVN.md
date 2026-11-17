# Sampling multiple times from a multivariate normal distribution

Sampling multiple times from a multivariate normal distribution

## Usage

``` r
multi_rMVN(mean_mat, precision)
```

## Arguments

- mean_mat:

  matrix containing n mean vectors

- precision:

  precision matrix

## Value

a matrix where each row corresponds to a draw from a MVN distribution
with mean equal to the rows of `mean_mat` and precision `precision`

## Examples

``` r
multi_rMVN(matrix(0, ncol=2, nrow=10), diag(2))
#>             [,1]       [,2]
#>  [1,]  0.2190175  1.7693545
#>  [2,] -1.1977317 -1.5669185
#>  [3,]  0.1319101 -0.8649049
#>  [4,] -0.4226498  1.0042502
#>  [5,] -1.7167643  0.7990631
#>  [6,]  1.2498249 -0.2669456
#>  [7,] -1.2147225  0.1601833
#>  [8,]  0.4441200  1.1724409
#>  [9,]  1.1455691 -0.5824341
#> [10,]  1.6181296 -0.9992989
```
