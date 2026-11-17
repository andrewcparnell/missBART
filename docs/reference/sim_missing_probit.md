# Simulating missing values from a probit regression model

Simulating missing values from a probit regression model

## Usage

``` r
sim_missing_probit(
  x,
  y,
  include_x = FALSE,
  include_y = FALSE,
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

- B:

  probit coefficients

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

- corR:

  correlation matrix

- zmod:

  probit latent variable

## Examples

``` r
# data <- sim_data_friedman(n = 1000, p = 3)
# missing <- sim_missing_probit(x = data$x, y = data$y, include_x = TRUE,
#                               include_y = TRUE, max_missing_prop = 0.99,
#                               min_node = 10, min_missing_prop = 0.7, intercept = TRUE)
```
