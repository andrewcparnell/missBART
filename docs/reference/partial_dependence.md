# Partial dependence and individual conditional expectation curves

Shows how the predictions of a fitted model change as one predictor
varies and the others are held at their observed values.
`type = "response"` gives the mean of each response as a covariate
varies. `type = "missingness"` gives the probability that each response
is missing as a covariate or a response varies.

## Usage

``` r
partial_dependence(
  object,
  variable,
  type = c("response", "missingness"),
  grid = NULL,
  n_grid = 20,
  n_ice = 200,
  draws = 200,
  probs = c(0.025, 0.975)
)

# S3 method for class 'missBART_pdp'
plot(x, ...)
```

## Arguments

- object:

  an object of class `"missBART"`. Except for the missingness model of a
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  fit, it must have been fitted with `keep_trees = TRUE`.

- variable:

  name or column number of the predictor to vary: a column of `x`, or
  for `type = "missingness"` also a response.

- type:

  which model to use, see above.

- grid:

  values of the predictor at which to evaluate the curves, on the scale
  of the data. The default is `n_grid` values spanning its observed
  range.

- n_grid:

  number of grid values if `grid` is not given.

- n_ice:

  number of observations (chosen at random) for which individual curves
  are computed and over which the partial dependence is averaged.

- draws:

  number of posterior samples to use (evenly spaced).

- probs:

  two probabilities giving the posterior interval of the partial
  dependence curve.

- x:

  an object returned by `partial_dependence()`.

- ...:

  passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

An object of class `"missBART_pdp"`: a list with the `grid`, the partial
dependence curve `pdp` (an array: grid values, responses, and
mean/lower/upper) and the individual curves `ice` (observations, grid
values, responses; posterior means). The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draws
one panel per response.

## Details

For the missingness model, missing responses are replaced by their
posterior means before the curves are computed. Probabilities are for
each response on its own, given the predictors.

## Examples

``` r
set.seed(1)
x <- matrix(runif(300), ncol = 3)
y <- cbind(sin(2 * pi * x[, 1]), x[, 2]) + matrix(rnorm(200, sd = 0.2), ncol = 2)
y[runif(100) < pnorm(2 * y[, 1] - 1), 1] <- NA
fit <- missBART2(x, y, n_reg_trees = 10, n_class_trees = 5, burn = 50, iters = 50,
                 thin = 1, mice_impute = FALSE, show_progress = FALSE, keep_trees = TRUE)
plot(partial_dependence(fit, "X1"))

plot(partial_dependence(fit, "Y1", type = "missingness"))
```
