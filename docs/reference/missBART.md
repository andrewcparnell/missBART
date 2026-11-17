# Joint BART models for responses with non-ignorable missingness

`missBART1()` and `missBART2()` fit the two selection models of Goh,
Soh, Parnell and Murphy (2024). Both use multivariate BART for the
responses. `missBART1()` models the missingness indicators with a
multivariate probit regression on `(1, x, y)`; `missBART2()` uses a
multivariate probit BART model on `(x, y)`. Missing responses are
imputed within the sampler. `missBARTprobit()` is the former name of
`missBART1()`.

## Usage

``` r
missBART1(
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
  scale = TRUE,
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = c(-0.5, 0.5),
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  probit_prior = probit_prior_list(),
  save_y_pred = TRUE,
  save_z = TRUE,
  ...
)

missBARTprobit(
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
  scale = TRUE,
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = c(-0.5, 0.5),
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  probit_prior = probit_prior_list(),
  save_y_pred = TRUE,
  save_z = TRUE,
  ...
)

missBART2(
  x,
  y,
  x_predict = NULL,
  n_reg_trees = 100,
  n_class_trees = 100,
  burn = 1000,
  iters = 1000,
  thin = 2,
  predict = TRUE,
  MH_sd = NULL,
  tree_prior_params = tree_list(...),
  hypers = hypers_list(...),
  scale = TRUE,
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = c(-0.5, 0.5),
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  save_y_pred = TRUE,
  save_z = TRUE,
  ...
)
```

## Arguments

- x:

  numeric matrix (or data frame) of covariates. Missing covariate values
  are allowed in `missBART2()` and [`mvBART()`](mvBART.md), where they
  are handled inside the trees; `missBART1()` needs complete covariates.

- y:

  numeric matrix (or vector) of responses, with `NA` for missing values.

- x_predict:

  optional covariates at which to make out-of-sample predictions.

- n_trees, n_reg_trees:

  number of trees in the BART model for the responses.

- burn:

  number of burn-in iterations.

- iters:

  number of posterior samples kept after burn-in and thinning.

- thin:

  thinning interval. The sampler runs `burn + thin * iters` iterations.

- predict:

  logical. Make predictions at `x_predict`?

- tree_prior_params:

  tree prior and proposal settings, see [`tree_list()`](tree_list.md).

- hypers:

  hyperparameters of the BART model for the responses, see
  [`hypers_list()`](hypers_list.md).

- scale:

  logical. Scale each response to `[-0.5, 0.5]` before fitting?
  `missBART1()` then also scales the covariates, so `B_post` refers to
  the scaled variables. Everything else is returned on the original
  scale.

- include_x, include_y:

  logical. Use `x` and/or `y` as predictors of the missingness?
  `include_y = FALSE` assumes the responses are missing at random.

- show_progress:

  logical. Print the iteration count?

- progress_every:

  how often to update the iteration count.

- pdp_range, make_pdp:

  for `missBART2()` with a univariate response: if `make_pdp = TRUE`,
  the latent missingness predictor is evaluated over a grid of (scaled)
  response values spanning `pdp_range` at every kept iteration and
  returned in `pdp_out`; apply
  [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html) to get the
  probability of being observed. Ignored by `missBART1()`, where the
  same curves follow directly from `B_post`.

- mice_impute:

  logical. Initialise the missing responses with
  [`mice::mice()`](https://amices.org/mice/reference/mice.html)?
  Otherwise they start at the centre of the observed range.

- ymis_update:

  how the Gaussian update of the missing responses treats rows that are
  only partly missing. `"conditional"` (the default) conditions on the
  observed responses in the same row. `"marginal"` draws from the
  corresponding margin of the full-row distribution, as written in
  Appendix C.1 of the paper and done by versions before 0.2.0. The two
  agree for univariate responses. Used by `missBART1()`, and by
  `missBART2()` when `include_y = FALSE`.

- probit_prior:

  gamma prior settings for the precisions of the probit coefficients in
  `missBART1()`, see [`probit_prior_list()`](probit_prior_list.md).

- save_y_pred:

  logical. Keep posterior predictive draws for the training data
  (`y_pred`)?

- save_z:

  logical. Keep the draws of the latent probit variables (`z_post`)?

- ...:

  further arguments passed to [`tree_list()`](tree_list.md) and
  [`hypers_list()`](hypers_list.md). Arguments of earlier versions that
  are no longer used are ignored.

- n_class_trees:

  number of trees in the probit BART model for the missingness
  indicators.

- MH_sd:

  standard deviation of the random-walk proposal for missing responses
  in `missBART2()`. Defaults to `0.5 / ncol(y)`.

## Value

An object of class `"missBART"`: a list with, among others,

- y_post:

  list of `iters` matrices with the BART mean function for each training
  row.

- y_pred:

  list of posterior predictive draws for the training rows.

- y_impute:

  list of vectors of imputed values, ordered as `which(is.na(y))`.

- new_y_hat, new_y_post:

  lists with the mean function and posterior predictive draws at
  `x_predict` (or `NA`).

- sigma_post:

  list of residual covariance matrices. `omega_post` is the same object
  under its old name.

- B_post, R_post, tau_b:

  `missBART1()` only: probit coefficients, correlation matrices, and
  prior precisions of the coefficients.

- var_imp:

  `missBART2()` only: matrix with the number of splitting rules on each
  predictor in the missingness trees, one row per kept iteration.

- var_imp_reg:

  the same for the regression trees.

- z_post:

  list of latent probit variables.

- y_miss_accept:

  `missBART2()` only: Metropolis-Hastings acceptance rate for each
  missing value.

- tree_accept:

  acceptance rates of the tree proposals.

Use [fitted()](missBART-methods.md), [predict()](missBART-methods.md),
[`imputed_values()`](imputed_values.md),
[`completed_data()`](imputed_values.md), [coef()](missBART-methods.md)
and [`varimp()`](varimp.md) to summarise a fit.

## References

Goh, Y. C., Soh, W. K., Parnell, A. C. and Murphy, K. (2024). Joint
models for handling non-ignorable missing data using Bayesian additive
regression trees: application to leaf photosynthetic traits data.
<doi:10.48550/arXiv.2412.14946>

## Examples

``` r
set.seed(1)
n <- 200
x <- matrix(runif(n * 3), ncol = 3)
y_true <- cbind(sin(2 * pi * x[, 1]), x[, 2]^2) + matrix(rnorm(2 * n, sd = 0.2), ncol = 2)
# large values of the first response are more likely to be missing
y <- y_true
y[runif(n) < pnorm(2 * y_true[, 1] - 1), 1] <- NA

fit <- missBART1(x, y, n_trees = 20, burn = 100, iters = 100, thin = 1,
                 mice_impute = FALSE, show_progress = FALSE)
fit
#> missBART1 fit: 200 observations, 2 responses
#> Missing responses: 53 of 400 (13.2%)
#> Posterior samples: 100 (burn-in 100, thinning 1)
#> Tree acceptance rate: regression 0.14
head(imputed_values(fit))
#>   row column      mean        sd      lower    upper
#> 1   1     Y1 0.7951750 0.2789181 0.28832332 1.290783
#> 2   5     Y1 0.8680268 0.2504344 0.44719060 1.335738
#> 3  10     Y1 0.6088303 0.2555736 0.05654392 1.042045
#> 4  11     Y1 0.9759243 0.2448018 0.49286643 1.434412
#> 5  12     Y1 0.7663950 0.2482717 0.28646051 1.246861
#> 6  22     Y1 0.7867143 0.2395778 0.37375796 1.200557
coef(fit)
#> , , mean
#> 
#>                   Y1         Y2
#> Intercept  1.6870607  4.6919918
#> X1         0.7832562  0.4425826
#> X2         1.3711738 -0.5957830
#> X3        -1.3374684  0.5115927
#> Y1        -4.1586286 -0.5503308
#> Y2         0.3426845 -3.3238790
#> 
#> , , lower
#> 
#>                   Y1        Y2
#> Intercept  1.0919603  2.526577
#> X1        -0.8079257 -1.760812
#> X2         0.3633939 -2.607402
#> X3        -2.2994658 -1.212002
#> Y1        -6.3317119 -2.812361
#> Y2        -1.1812488 -7.373290
#> 
#> , , upper
#> 
#>                   Y1        Y2
#> Intercept  2.1490277 6.9990208
#> X1         2.2655507 2.1798183
#> X2         2.3506381 1.0384469
#> X3        -0.4136713 2.3724394
#> Y1        -1.9123977 1.7032939
#> Y2         2.1912406 0.9553731
#> 

fit2 <- missBART2(x, y, n_reg_trees = 20, n_class_trees = 10, burn = 100,
                  iters = 100, thin = 1, mice_impute = FALSE,
                  show_progress = FALSE)
varimp(fit2)
#>        X1        X2        X3        Y1        Y2 
#> 0.2614935 0.1336365 0.1178880 0.3044160 0.1825660 
```
