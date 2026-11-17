# Variable interactions

Counts how often two predictors are used in consecutive splitting rules
on the same branch of a tree (a parent and one of its children), as a
measure of how often they are used together.

## Usage

``` r
interactions(object, trees = c("regression", "missingness"))
```

## Arguments

- object:

  an object of class `"missBART"` fitted with `keep_trees = TRUE`.

- trees:

  which set of trees: `"regression"` or `"missingness"` (only for
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  fits).

## Value

A symmetric matrix. The off-diagonal entries are the average number of
parent and child pairs per posterior sample that use the two predictors,
and the diagonal entries the average number of splitting rules that use
the predictor.
