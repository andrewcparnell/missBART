# Simulate

Simulate

## Usage

``` r
sim_data_friedman(n, p = 1, scale_par = 1, omega_diag = 1, Omega = NULL)
```

## Arguments

- n:

  number of observations

- p:

  number of responses

- scale_par:

  scale_par

- omega_diag:

  precision ~ rWishart(1, p + 1, diag(omega_diag, p))

- Omega:

  allows user-supplied precision matrix

## Value

A list with the following components

- y:

  responses

- x:

  covariates

- ome:

  precision matrix

- mean:

  mean vector

- p:

  number of responses

- q:

  number of covariates (always 5)

## Examples

``` r
# data <- sim_data_friedman(n = 1000, p = 3)
# true_trees_data <- data$true_trees
# y <- y_original <- data$y
# x <- data$x
# ome <- data$ome
```
