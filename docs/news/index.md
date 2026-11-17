# Changelog

## missBART 0.2.0

### Speed

- The samplers behind
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  (formerly
  [`missBARTprobit()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)),
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  and
  [`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md)
  are now written in C++ and are roughly 80 to 350 times faster than
  version 0.1.0; see the README for timings.
- Results are reproducible with
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but the draws
  differ from those of version 0.1.0 for the same seed.

### Changes to the samplers

These change the results relative to version 0.1.0.

- **Missingness indicator.** The models are now written with `M = 1` for
  a missing response, and the latent probit variable is positive for
  missing responses. The probit coefficients of
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  (`B_post`, [`coef()`](https://rdrr.io/r/stats/coef.html)) and the
  latent variables (`z_post`, `pdp_out`) therefore have the opposite
  sign to version 0.1.0: a positive coefficient now means that larger
  values make the response more likely to be *missing*.
  [`sim_missing_probit()`](https://andrewcparnell.github.io/missBART/reference/sim_missing_probit.md)
  and
  [`sim_missing_trees()`](https://andrewcparnell.github.io/missBART/reference/sim_missing_trees.md)
  use the same coding (their `m` is 1 for missing), and their
  `min_missing_prop` and `max_missing_prop` arguments now really bound
  the proportion missing (defaults 0.1 and 0.4).
- **Standardisation of the responses.** The new default
  `scale = "center"` follows the `BART` package (Sparapani, Spanbauer
  and McCulloch, 2021): each response has an offset equal to the mean of
  its observed values, and the terminal node parameters are centred at
  zero with a standard deviation proportional to the range of the
  response. `scale = "range"` (or `TRUE`) gives the `[-0.5, 0.5]`
  scaling of version 0.1.0, and `scale = "none"` (or `FALSE`) no
  scaling. The probit BART model of
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  has an offset equal to the probit of the proportion of missing values.
- **Proposal ratio.** The Metropolis-Hastings step for the trees now
  includes the proposal ratio of the grow and prune moves. Version 0.1.0
  accepted on the ratio of prior times likelihood alone.
- **Missing responses in partly observed rows.** The Gaussian update of
  the missing responses
  ([`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md),
  and
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  with `include_y = FALSE`) now conditions on the observed responses in
  the same row. Version 0.1.0, and Appendix C.1 of the paper as written,
  draw from the margin of the full-row distribution instead, which
  ignores the residual correlation between observed and missing
  responses of the same row. `ymis_update = "marginal"` restores the old
  behaviour. The two agree for univariate responses.
- **Node marginal likelihood.** For more than one response, version
  0.1.0 computed the quadratic form in the node marginal likelihood as
  `sum(Sigma_inv %*% mu^2)` rather than `t(mu) %*% Sigma_inv %*% mu`.
  The formula in the paper is now used.
- **Inverse-Wishart degrees of freedom** in the parameter-expanded
  update of the probit correlation matrix are `n + p + 1` (Talhouk et
  al., 2012, in the usual parameterisation) instead of `n - p + 3`.
- **Prior precisions of the probit coefficients.** One precision is
  drawn for all covariates and one for all responses, as in Appendix B
  of the paper. Version 0.1.0 drew a separate value for every
  coefficient from the pooled full conditional.
- **Latent probit variables** in
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  with more than one response are updated one coordinate at a time from
  their univariate truncated normal full conditionals, instead of by an
  exact draw from the multivariate truncated normal with
  [`TruncatedNormal::rtmvnorm()`](https://rdrr.io/pkg/TruncatedNormal/man/rtmvnorm.html)
  for every row.
- **Tree proposals.** Split values are always drawn within the range of
  the covariate that is still available at the node. Version 0.1.0 used
  the interval (0, 1) after a change move at the root, and restricted
  trees whose deepest node is at depth 3 to grow and prune moves.
- Trees can send missing covariate values left or right at every split,
  and missingness indicators are added only for covariates that have
  missing values.
- [`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md)
  now respects `hypers$kappa`.

### Interface

- New name
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md),
  matching the paper;
  [`missBARTprobit()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  still works.
- Fits have class `"missBART"` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`coef()`](https://rdrr.io/r/stats/coef.html) methods, and the helpers
  [`imputed_values()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md),
  [`completed_data()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md)
  and
  [`varimp()`](https://andrewcparnell.github.io/missBART/reference/varimp.md).
- New arguments `ymis_update`, `probit_prior` (see
  [`probit_prior_list()`](https://andrewcparnell.github.io/missBART/reference/probit_prior_list.md)),
  `save_y_pred` and `save_z`.
- The posterior draws keep their names and list-of-matrices layout
  (`y_post`, `y_pred`, `y_impute`, `new_y_post`, `B_post`, `R_post`,
  `z_post`, …), so existing analysis scripts keep working, with these
  exceptions:
  - `omega_post` (also available as `sigma_post`) always holds full
    residual covariance matrices on the original scale. It used to hold
    only the variances when `scale = TRUE` and there was more than one
    response. (This is issue
    [\#1](https://github.com/yongchengoh/missBART/issues/1) on GitHub.)
  - `new_y_hat` is new and holds the BART mean function at `x_predict`;
    `new_y_post` holds posterior predictive draws, as before.
  - `var_imp` is a matrix of split counts (iterations by predictors)
    rather than a list of tables; `var_imp_reg` is the same for the
    regression trees.
  - `y_miss_accept` holds one acceptance rate per missing value rather
    than an iterations-by-values matrix.
  - `tau_b` no longer has a leading row of `NA`.
  - `pdp_range` is given on the scale of the data (it used to be on the
    `[-0.5, 0.5]` scale) and defaults to the observed range; the grid is
    returned as `pdp_grid`.
  - `pdp_out`
    ([`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md),
    univariate response) has one element per kept iteration, and no
    longer accumulates across iterations.
    [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
    no longer returns `pdp_out`; compute the curves from `B_post`.
  - The trees and terminal node parameters (`reg_trees`, `class_trees`,
    `reg_mu`, `class_mu`, `accepted_trees`) are no longer returned.
- The arguments `true_trees_data`, `true_trees_missing`,
  `true_change_points` and `true_change_points_miss` are ignored.
- Dependencies on `CholWishart`, `extraDistr`, `TruncatedNormal`,
  `Matrix`, `MASS` and `ggplot2` were dropped; `Rcpp` and
  `RcppArmadillo` were added.
- The package has a pkgdown website (built into `docs/`) with a quick
  start vignette and further articles, including a re-run of the
  simulation studies in the paper.
- The unexported plotting helpers `plot_posterior()` and
  `plot_plant_posterior()` were removed.
