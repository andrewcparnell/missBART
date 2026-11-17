# Multivariate BART for complete responses

Fits the multivariate BART model used for the responses in
[`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
and
[`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
to data without missing responses.

## Usage

``` r
mvBART(
  x,
  y,
  x_predict = NULL,
  n_trees = 100,
  burn = 1000,
  iters = 1000,
  thin = 2,
  predict = TRUE,
  tree_prior_params = tree_list(...),
  hypers = hypers_list(...),
  scale = "center",
  show_progress = TRUE,
  progress_every = 10,
  save_y_pred = TRUE,
  keep_trees = FALSE,
  n_chains = 1,
  cores = 1,
  ...
)
```

## Arguments

- x:

  numeric matrix (or data frame) of covariates. Missing covariate values
  are allowed in
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  and `mvBART()`, where they are handled inside the trees;
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  needs complete covariates.

- y:

  numeric matrix (or vector) of responses without missing values.

- x_predict:

  optional covariates at which to make out-of-sample predictions.

- n_trees:

  number of trees.

- burn:

  number of burn-in iterations.

- iters:

  number of posterior samples kept after burn-in and thinning.

- thin:

  thinning interval. The sampler runs `burn + thin * iters` iterations.

- predict:

  logical. Make predictions at `x_predict`?

- tree_prior_params:

  tree prior and proposal settings, see
  [`tree_list()`](https://andrewcparnell.github.io/missBART/reference/tree_list.md).

- hypers:

  hyperparameters of the BART model for the responses, see
  [`hypers_list()`](https://andrewcparnell.github.io/missBART/reference/hypers_list.md).

- scale:

  how the responses are standardised before fitting; see Details.
  `"center"` (the default) subtracts the mean of the observed values of
  each response and divides by their range, `"range"` (or `TRUE`) maps
  each response to `[-0.5, 0.5]` as in Chipman et al. (2010) and
  versions before 0.2.0, and `"none"` (or `FALSE`) leaves them alone.
  Unless `scale = "none"`,
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  also scales the covariates to `[-0.5, 0.5]`, so `B_post` refers to
  standardised variables. Everything else is returned on the original
  scale.

- show_progress:

  logical. Print the iteration count?

- progress_every:

  how often to update the iteration count.

- save_y_pred:

  logical. Keep posterior predictive draws for the training data
  (`y_pred`)?

- keep_trees:

  logical. Keep the trees of every posterior sample? This is needed for
  [predict()](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  at new covariates,
  [`partial_dependence()`](https://andrewcparnell.github.io/missBART/reference/partial_dependence.md)
  and
  [`interactions()`](https://andrewcparnell.github.io/missBART/reference/interactions.md),
  and can take a lot of memory.

- n_chains:

  number of chains to run from different starting values. The samples of
  all chains are returned together, with `chain` giving the chain of
  each sample; see
  [`chain_diagnostics()`](https://andrewcparnell.github.io/missBART/reference/chain_diagnostics.md).

- cores:

  number of chains to run in parallel (not on Windows).

- ...:

  further arguments passed to
  [`tree_list()`](https://andrewcparnell.github.io/missBART/reference/tree_list.md)
  and
  [`hypers_list()`](https://andrewcparnell.github.io/missBART/reference/hypers_list.md).
  Arguments of earlier versions that are no longer used are ignored.

## Value

An object of class `"missBART"`, see
[`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md).

## Examples

``` r
set.seed(1)
data <- sim_data_friedman(n = 100, p = 2)
fit <- mvBART(data$x, data$y, n_trees = 20, burn = 100, iters = 100,
              thin = 1, show_progress = FALSE)
plot(data$mean[, 1], fitted(fit)[, 1])
```
