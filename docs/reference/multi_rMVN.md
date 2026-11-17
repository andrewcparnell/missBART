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
#>              [,1]       [,2]
#>  [1,] -1.36979327 -1.0047133
#>  [2,] -1.08743262 -0.5804202
#>  [3,]  0.40813887  0.5149419
#>  [4,] -0.02564099  1.0283955
#>  [5,] -0.32516783  1.0633289
#>  [6,]  0.50440683  0.2341026
#>  [7,]  0.19177092  1.5400044
#>  [8,] -2.03238685  0.7516753
#>  [9,] -1.26225023 -0.3909854
#> [10,]  0.71354954  2.1613693
```
