# Prior for the probit regression coefficients of missBART1

In [`missBART1()`](missBART.md) the probit coefficients have prior
precisions `tau_0` (intercept), `tau_X` (covariates) and `tau_Y`
(responses) with `Gamma(a, b)` priors (shape `a`, rate `b`). A prior
that favours large `tau_Y` shrinks the coefficients of the responses
towards zero, i.e. favours data that are missing at random.

## Usage

``` r
probit_prior_list(a0 = 2, b0 = 1, ax = NULL, bx = NULL, ay = NULL, by = NULL)
```

## Arguments

- a0, b0:

  shape and rate for the intercept precision.

- ax, bx:

  shape and rate for the covariate precision. The defaults are
  `ax = q + 1` and `bx = p + q + 1`, with `q` covariates and `p`
  responses.

- ay, by:

  shape and rate for the response precision. The defaults are `ay = 1`
  and `by = 1`.

## Value

list of prior parameters

## Details

The defaults are those of the package code before version 0.2.0. Section
5.1 of the paper instead states `(ax, bx) = (1 + q, 1)` and
`(ay, by) = (1 + p + q, 1)`; use
`probit_prior_list(bx = 1, ay = 1 + p + q)` for that setting.

## Examples

``` r
probit_prior_list(ay = 10)
#> $a0
#> [1] 2
#> 
#> $b0
#> [1] 1
#> 
#> $ax
#> NULL
#> 
#> $bx
#> NULL
#> 
#> $ay
#> [1] 10
#> 
#> $by
#> NULL
#> 
```
