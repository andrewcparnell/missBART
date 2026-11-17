# Simulating data from trees

Simulating data from trees

## Usage

``` r
sim_data_trees(
  n,
  p,
  q,
  min_x = 0,
  max_x = 1,
  trees = 1,
  ome = NULL,
  kappa = NULL,
  splits = NULL,
  ...
)
```

## Arguments

- n:

  number of observations

- p:

  number of responses

- q:

  number of covariates

- min_x:

  min_x

- max_x:

  max_x

- trees:

  number of trees

- ome:

  precision matrix

- kappa:

  prior scaling factor

- splits:

  tree splits

- ...:

  Catches unused arguments

## Value

A list with the following components

- y:

  responses

- x:

  covariates

- ome:

  precision matrix

- true_trees:

  trees

## Examples

``` r
# data <- sim_data_trees(n = 100, p = 3, q = 4, trees = 6)
# y <- y_original <- data$y
# x <- data$x
# ome <- data$ome
# true_trees_data <- data$true_trees
```
