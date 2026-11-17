# Package index

## Models

Fit the joint models, or multivariate BART on its own.

- [`missBART1()`](missBART.md) [`missBARTprobit()`](missBART.md)
  [`missBART2()`](missBART.md) : Joint BART models for responses with
  non-ignorable missingness
- [`mvBART()`](mvBART.md) : Multivariate BART for complete responses

## Summarising a fit

- [`print(`*`<missBART>`*`)`](missBART-methods.md)
  [`summary(`*`<missBART>`*`)`](missBART-methods.md)
  [`fitted(`*`<missBART>`*`)`](missBART-methods.md)
  [`predict(`*`<missBART>`*`)`](missBART-methods.md)
  [`coef(`*`<missBART>`*`)`](missBART-methods.md) : Methods for missBART
  fits
- [`imputed_values()`](imputed_values.md)
  [`completed_data()`](imputed_values.md) : Summaries of the imputed
  responses
- [`varimp()`](varimp.md) : Variable importance

## Priors and settings

- [`hypers_list()`](hypers_list.md) : BART prior parameters
- [`tree_list()`](tree_list.md) : Tree prior parameters
- [`probit_prior_list()`](probit_prior_list.md) : Prior for the probit
  regression coefficients of missBART1

## Simulating data

- [`sim_data_friedman()`](sim_data_friedman.md) : Simulate
- [`sim_data_trees()`](sim_data_trees.md) : Simulating data from trees
- [`sim_missing_probit()`](sim_missing_probit.md) : Simulating missing
  values from a probit regression model
- [`sim_missing_trees()`](sim_missing_trees.md) : Simulate missing data
  from a BART model
- [`propose_tree()`](propose_tree.md) : A function for proposing an
  updated tree structure

## Utilities

- [`scale_bart()`](scale_bart.md) [`unscale()`](scale_bart.md) : Scaling
  responses to and from the BART scale
- [`rMVN()`](rMVN.md) : Sampling once from a multivariate normal
  distribution with mean mu and precision Q
- [`rMVNc()`](rMVNc.md) : Sampling from a multivariate normal
  distribution with precision Q and mean equal to the solution of Q mu =
  b
- [`multi_rMVN()`](multi_rMVN.md) : Sampling multiple times from a
  multivariate normal distribution
- [`matrnorm()`](matrnorm.md) : Sampling from a matrix normal
  distribution MV(mu, U, V)
- [`missBART-package`](missBART-package.md) : missBART: joint models for
  non-ignorable missing responses using BART
