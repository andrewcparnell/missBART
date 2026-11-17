# missBART 0.2.0

## Speed

* The samplers behind `missBART1()` (formerly `missBARTprobit()`),
  `missBART2()` and `mvBART()` are now written in C++ and are roughly 80 to
  350 times faster than version 0.1.0; see the README for timings.
* Results are reproducible with `set.seed()`, but the draws differ from those
  of version 0.1.0 for the same seed.

## Changes to the samplers

These change the results relative to version 0.1.0.

* **Missingness indicator.** The models are now written with
  `M = 1` for a missing response, and the latent probit variable is positive
  for missing responses. The probit coefficients of `missBART1()` (`B_post`,
  `coef()`) and the latent variables (`z_post`, `pdp_out`) therefore have the
  opposite sign to version 0.1.0: a positive coefficient now means that larger
  values make the response more likely to be *missing*. `sim_missing_probit()`
  and `sim_missing_trees()` use the same coding (their `m` is 1 for missing),
  and their `min_missing_prop` and `max_missing_prop` arguments now really
  bound the proportion missing (defaults 0.1 and 0.4).
* **Standardisation of the responses.** The new default `scale = "center"`
  follows the `BART` package (Sparapani, Spanbauer and McCulloch, 2021): each
  response has an offset equal to the mean of its observed values, and the
  terminal node parameters are centred at zero with a standard deviation
  proportional to the range of the response. `scale = "range"` (or `TRUE`)
  gives the `[-0.5, 0.5]` scaling of version 0.1.0, and `scale = "none"` (or
  `FALSE`) no scaling. The probit BART model of `missBART2()` has an offset
  equal to the probit of the proportion of missing values.
* **Proposal ratio.** The Metropolis-Hastings step for the trees now includes
  the proposal ratio of the grow and prune moves. Version 0.1.0 accepted on
  the ratio of prior times likelihood alone.
* **Missing responses in partly observed rows.** The Gaussian update of the
  missing responses (`missBART1()`, and `missBART2()` with
  `include_y = FALSE`) now conditions on the observed responses in the same
  row. Version 0.1.0, and Appendix C.1 of the paper as written, draw from the
  margin of the full-row distribution instead, which ignores the residual
  correlation between observed and missing responses of the same row.
  `ymis_update = "marginal"` restores the old behaviour. The two agree for
  univariate responses.
* **Node marginal likelihood.** For more than one response, version 0.1.0
  computed the quadratic form in the node marginal likelihood as
  `sum(Sigma_inv %*% mu^2)` rather than `t(mu) %*% Sigma_inv %*% mu`. The
  formula in the paper is now used.
* **Inverse-Wishart degrees of freedom** in the parameter-expanded update of
  the probit correlation matrix are `n + p + 1` (Talhouk et al., 2012, in the
  usual parameterisation) instead of `n - p + 3`.
* **Prior precisions of the probit coefficients.** One precision is drawn for
  all covariates and one for all responses, as in Appendix B of the paper.
  Version 0.1.0 drew a separate value for every coefficient from the pooled
  full conditional.
* **Latent probit variables** in `missBART1()` with more than one response are
  updated one coordinate at a time from their univariate truncated normal full
  conditionals, instead of by an exact draw from the multivariate truncated
  normal with `TruncatedNormal::rtmvnorm()` for every row.
* **Tree proposals.** Split values are always drawn within the range of the
  covariate that is still available at the node. Version 0.1.0 used the
  interval (0, 1) after a change move at the root, and restricted trees whose
  deepest node is at depth 3 to grow and prune moves.
* Trees can send missing covariate values left or right at every split, and
  missingness indicators are added only for covariates that have missing
  values.
* `mvBART()` now respects `hypers$kappa`.

## New features

* `keep_trees = TRUE` stores the trees of every posterior sample
  (`reg_trees`, `class_trees`). With them, `predict()` works at new
  covariates (`newdata`), `partial_dependence()` gives partial dependence and
  individual conditional expectation curves for the responses and for the
  missingness probabilities, and `interactions()` counts how often pairs of
  predictors are used in consecutive splitting rules.
* `n_chains` runs several chains from different starting values (in parallel
  with `cores`), and `chain_diagnostics()` compares them. `print()` warns
  when the chains disagree. Selection models can be multimodal, so this is
  recommended before interpreting a `missBART1()` fit.
* `plot()` compares the distributions of imputed and observed values
  (`which = "imputed"`) or plots fitted against observed values with the
  squared correlation (`which = "fit"`).
* The assumption of ignorable missingness is available as a benchmark with
  the same response model: `include_y = FALSE`.

## Interface

* New name `missBART1()`, matching the paper; `missBARTprobit()` still works.
* Fits have class `"missBART"` with `print()`, `summary()`, `fitted()`,
  `predict()` and `coef()` methods, and the helpers `imputed_values()`,
  `completed_data()` and `varimp()`.
* New arguments `ymis_update`, `probit_prior` (see `probit_prior_list()`),
  `save_y_pred` and `save_z`.
* The posterior draws keep their names and list-of-matrices layout
  (`y_post`, `y_pred`, `y_impute`, `new_y_post`, `B_post`, `R_post`,
  `z_post`, ...), so existing analysis scripts keep working, with these
  exceptions:
    * `omega_post` (also available as `sigma_post`) always holds full residual
      covariance matrices on the original scale. It used to hold only the
      variances when `scale = TRUE` and there was more than one response.
      (This is issue #1 on GitHub.)
    * `new_y_hat` is new and holds the BART mean function at `x_predict`;
      `new_y_post` holds posterior predictive draws, as before.
    * `var_imp` is a matrix of split counts (iterations by predictors) rather
      than a list of tables; `var_imp_reg` is the same for the regression
      trees.
    * `y_miss_accept` holds one acceptance rate per missing value rather than
      an iterations-by-values matrix.
    * `tau_b` no longer has a leading row of `NA`.
    * `pdp_range` is given on the scale of the data (it used to be on the
      `[-0.5, 0.5]` scale) and defaults to the observed range; the grid is
      returned as `pdp_grid`.
    * `pdp_out` (`missBART2()`, univariate response) has one element per kept
      iteration, and no longer accumulates across iterations. `missBART1()`
      no longer returns `pdp_out`; compute the curves from `B_post`.
    * The trees are returned only with `keep_trees = TRUE`, in a compact
      matrix format, and `reg_mu`, `class_mu` and `accepted_trees` are no
      longer returned.
* The arguments `true_trees_data`, `true_trees_missing`, `true_change_points`
  and `true_change_points_miss` are ignored.
* Dependencies on `CholWishart`, `extraDistr`, `TruncatedNormal`, `Matrix`,
  `MASS` and `ggplot2` were dropped; `Rcpp` and `RcppArmadillo` were added.
* The package has a pkgdown website (built into `docs/`) with a quick start
  vignette and further articles, including a re-run of the simulation studies
  in the paper.
* The unexported plotting helpers `plot_posterior()` and
  `plot_plant_posterior()` were removed.
