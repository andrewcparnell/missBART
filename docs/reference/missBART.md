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
  scale = "center",
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = NULL,
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  probit_prior = probit_prior_list(),
  save_y_pred = TRUE,
  save_z = TRUE,
  keep_trees = FALSE,
  n_chains = 1,
  cores = 1,
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
  scale = "center",
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = NULL,
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  probit_prior = probit_prior_list(),
  save_y_pred = TRUE,
  save_z = TRUE,
  keep_trees = FALSE,
  n_chains = 1,
  cores = 1,
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
  scale = "center",
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  pdp_range = NULL,
  make_pdp = FALSE,
  mice_impute = TRUE,
  ymis_update = c("conditional", "marginal"),
  save_y_pred = TRUE,
  save_z = TRUE,
  keep_trees = FALSE,
  n_chains = 1,
  cores = 1,
  ...
)
```

## Arguments

- x:

  numeric matrix (or data frame) of covariates. Missing covariate values
  are allowed in `missBART2()` and
  [`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md),
  where they are handled inside the trees; `missBART1()` needs complete
  covariates.

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
  Unless `scale = "none"`, `missBART1()` also scales the covariates to
  `[-0.5, 0.5]`, so `B_post` refers to standardised variables.
  Everything else is returned on the original scale.

- include_x, include_y:

  logical. Use `x` and/or `y` as predictors of the missingness?
  `include_y = FALSE` assumes the responses are missing at random.

- show_progress:

  logical. Print the iteration count?

- progress_every:

  how often to update the iteration count.

- pdp_range, make_pdp:

  for `missBART2()` with a univariate response: if `make_pdp = TRUE`,
  the latent missingness predictor is evaluated over a grid of response
  values spanning `pdp_range` (on the scale of the data; the default is
  the observed range) at every kept iteration and returned in `pdp_out`,
  with the grid in `pdp_grid`; apply
  [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html) to get the
  probability of being missing. Ignored by `missBART1()`, where the same
  curves follow directly from `B_post`.

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
  `missBART1()`, see
  [`probit_prior_list()`](https://andrewcparnell.github.io/missBART/reference/probit_prior_list.md).

- save_y_pred:

  logical. Keep posterior predictive draws for the training data
  (`y_pred`)?

- save_z:

  logical. Keep the draws of the latent probit variables (`z_post`)?

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

  list of latent probit variables (positive for missing responses).

- y_miss_accept:

  `missBART2()` only: Metropolis-Hastings acceptance rate for each
  missing value.

- tree_accept:

  acceptance rates of the tree proposals.

- reg_trees, class_trees:

  if `keep_trees = TRUE`, the regression and missingness trees of every
  sample: a matrix with one row per node giving the tree, the row of the
  left child (0 for a terminal node; the right child is the next row),
  the splitting variable and value, whether missing covariate values go
  left, and the terminal node parameters.

- chain:

  the chain each sample came from.

Use
[fitted()](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md),
[predict()](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md),
[plot()](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md),
[`imputed_values()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md),
[`completed_data()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md),
[coef()](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md),
[`varimp()`](https://andrewcparnell.github.io/missBART/reference/varimp.md),
[`partial_dependence()`](https://andrewcparnell.github.io/missBART/reference/partial_dependence.md),
[`interactions()`](https://andrewcparnell.github.io/missBART/reference/interactions.md)
and
[`chain_diagnostics()`](https://andrewcparnell.github.io/missBART/reference/chain_diagnostics.md)
to summarise a fit.

## Details

**Missingness indicator.** Write \\M\_{ij} = 1\\ if response \\j\\ of
observation \\i\\ is missing and 0 if it is observed. Both models use a
latent variable \\M^\*\_{ij}\\ that is positive exactly when \\M\_{ij} =
1\\. In `missBART1()` its mean is \\B^\top Z_i\\, so a positive
coefficient in `B_post` means that larger values of the predictor make
the response more likely to be *missing*. (Versions before 0.2.0, and
the first version of the paper, used the opposite coding.)

**Prior for the responses.** Following the `BART` package (Sparapani,
Spanbauer and McCulloch, 2021), the default `scale = "center"` fits
\\Y\_{ij} = \mu\_{0j} + f_j(X_i) + \epsilon\_{ij}\\, with \\\mu\_{0j}\\
the mean of the observed values of response \\j\\ and terminal node
parameters centred at zero with standard deviation \\(\max_j - \min_j) /
(2 k \sqrt{K})\\ for \\K\\ trees, where \\k = \sqrt{\kappa / (4K)}\\ is
1.96 by default (see
[`hypers_list()`](https://andrewcparnell.github.io/missBART/reference/hypers_list.md)).
Internally this is done by standardising each response, which is
equivalent. The probit BART model of `missBART2()` likewise has an
offset \\\mu\_{0j} = \Phi^{-1}(\bar M_j)\\, the probit of the proportion
of missing values of response \\j\\.

**Tree updates.** Each tree is updated by a grow, prune or change
proposal accepted with a Metropolis-Hastings step. The acceptance ratio
includes the proposal ratio of the grow and prune moves.

## References

Chipman, H. A., George, E. I. and McCulloch, R. E. (2010). BART:
Bayesian additive regression trees. *The Annals of Applied Statistics*,
4(1), 266-298.

Sparapani, R., Spanbauer, C. and McCulloch, R. (2021). Nonparametric
machine learning and efficient computation with Bayesian additive
regression trees: the BART R package. *Journal of Statistical Software*,
97(1), 1-66. <doi:10.18637/jss.v097.i01>

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
#> Tree acceptance rate: regression 0.15
head(imputed_values(fit))
#>   row column      mean        sd     lower    upper
#> 1   1     Y1 1.0008367 0.2717892 0.4486066 1.489274
#> 2   5     Y1 1.0989380 0.2334634 0.6631485 1.570068
#> 3  10     Y1 0.6120019 0.2457835 0.1716297 1.085211
#> 4  11     Y1 0.7820146 0.3068687 0.1985202 1.373719
#> 5  12     Y1 0.7540756 0.2887843 0.1828341 1.431624
#> 6  22     Y1 0.9899950 0.2793448 0.4979206 1.649501
coef(fit)
#> , , mean
#> 
#>                   Y1         Y2
#> Intercept -1.9656382 -3.1954136
#> X1        -0.5265129 -0.3406855
#> X2        -0.8873524  0.7625772
#> X3         1.3557372 -0.5982172
#> Y1         4.3799813  0.6950918
#> Y2        -0.7781849 -1.0852781
#> 
#> , , lower
#> 
#>                   Y1        Y2
#> Intercept -2.4819605 -4.034056
#> X1        -1.7242043 -2.254201
#> X2        -2.0192534 -0.859558
#> X3         0.4566221 -2.106015
#> Y1         3.0492314 -1.210050
#> Y2        -2.9935171 -3.431173
#> 
#> , , upper
#> 
#>                   Y1         Y2
#> Intercept -1.4733225 -2.5253265
#> X1         0.7302450  1.3884587
#> X2         0.5039515  2.5225349
#> X3         2.2404024  0.8156959
#> Y1         5.8093208  1.9576431
#> Y2         1.2527833  3.3305669
#> 

fit2 <- missBART2(x, y, n_reg_trees = 20, n_class_trees = 10, burn = 100,
                  iters = 100, thin = 1, mice_impute = FALSE,
                  show_progress = FALSE)
varimp(fit2)
#>        X1        X2        X3        Y1        Y2 
#> 0.2460134 0.2079455 0.1434260 0.2065475 0.1960676 
```
