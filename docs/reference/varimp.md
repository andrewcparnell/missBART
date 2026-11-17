# Variable importance

The proportion of splitting rules that use each predictor, averaged over
the posterior samples. For the missingness trees of
[`missBART2()`](missBART.md), splits concentrated on the responses point
towards data that are missing not at random.

## Usage

``` r
varimp(object, trees = c("missingness", "regression"))
```

## Arguments

- object:

  an object of class `"missBART"`.

- trees:

  which set of trees: `"missingness"` (only for
  [`missBART2()`](missBART.md) fits) or `"regression"`.

## Value

A named vector of proportions.
