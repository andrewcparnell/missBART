# Simulate missing data from a BART model

Simulate missing data from a BART model

## Usage

``` r
sim_missing_trees(
  x,
  y,
  trees = 1,
  include_x = FALSE,
  include_y = TRUE,
  min_missing_prop = 0.6,
  max_missing_prop = 0.9,
  ...
)
```

## Arguments

- x:

  regression covariates

- y:

  regression response

- trees:

  number of trees

- include_x:

  logical. Include x in probit model?

- include_y:

  logical. Include y in probit model?

- min_missing_prop:

  minimum proportion of missingness

- max_missing_prop:

  maximum proportion of missingness

- ...:

  Catches unused arguments

## Value

A list with the following components

- m:

  missingness indicators

- missing_y:

  inputted y with missing values deleted

- missing_prop:

  proportion of missingness per response variable

- missing_index:

  indices of missing cells in `missing_y`

- obs_index:

  indices of observed cells in `missing_y`

- true_trees:

  trees

- z:

  latent variable

## Examples

``` r
# data <- sim_data_trees(n = 100, p = 3, q = 4, trees = 6)
# missing <- sim_missing_trees(x = data$x, y = data$y, include_x = FALSE,
#                              include_y = TRUE, max_missing_prop = 0.99,
#                              min_node = 5, min_missing_prop = 0.7, trees = 1)
```
