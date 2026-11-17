# BART prior parameters

BART prior parameters

## Usage

``` r
hypers_list(
  mu0 = 0,
  kappa = NULL,
  alpha = NULL,
  V = NULL,
  df = 10,
  q = 0.75,
  tau_b = 100,
  ...
)
```

## Arguments

- mu0:

  hyperparameter

- kappa:

  hyperparameter

- alpha:

  hyperparameter

- V:

  hyperparameter

- df:

  hyperparameter

- q:

  hyperparameter

- tau_b:

  not used; kept for compatibility with earlier versions

- ...:

  Catches unused arguments

## Value

list of hyperparameters

## Examples

``` r
hypers_list(mu0 = 5, kappa = 2, alpha = 1, V = 2)
#> $mu0
#> [1] 5
#> 
#> $kappa
#> [1] 2
#> 
#> $alpha
#> [1] 1
#> 
#> $V
#> [1] 2
#> 
#> $df
#> [1] 10
#> 
#> $q
#> [1] 0.75
#> 
#> $tau_b
#> [1] 100
#> 
```
