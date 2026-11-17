# Agreement between chains

For a fit with `n_chains` greater than one, compares summaries between
the chains. Selection models can have several modes, in which case
chains started from different values give different answers; see the
real data article on the package website.

## Usage

``` r
chain_diagnostics(object)
```

## Arguments

- object:

  an object of class `"missBART"` fitted with `n_chains > 1`.

## Value

A list with `by_chain`, a matrix of the posterior mean of each quantity
in each chain, and `rhat`, the potential scale reduction factor of
Gelman and Rubin for each quantity (values well above 1 mean that the
chains disagree). The quantities are the mean imputed value of each
response, the residual variances and, for
[`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md),
the probit coefficients of the responses.
