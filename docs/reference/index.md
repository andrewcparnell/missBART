# Package index

## Models

Fit the joint models, or multivariate BART on its own.

- [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  [`missBARTprobit()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  : Joint BART models for responses with non-ignorable missingness
- [`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md)
  : Multivariate BART for complete responses

## Summarising a fit

- [`print(`*`<missBART>`*`)`](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  [`summary(`*`<missBART>`*`)`](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  [`fitted(`*`<missBART>`*`)`](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  [`predict(`*`<missBART>`*`)`](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  [`coef(`*`<missBART>`*`)`](https://andrewcparnell.github.io/missBART/reference/missBART-methods.md)
  : Methods for missBART fits
- [`imputed_values()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md)
  [`completed_data()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md)
  : Summaries of the imputed responses
- [`varimp()`](https://andrewcparnell.github.io/missBART/reference/varimp.md)
  : Variable importance

## Priors and settings

- [`hypers_list()`](https://andrewcparnell.github.io/missBART/reference/hypers_list.md)
  : BART prior parameters
- [`tree_list()`](https://andrewcparnell.github.io/missBART/reference/tree_list.md)
  : Tree prior parameters
- [`probit_prior_list()`](https://andrewcparnell.github.io/missBART/reference/probit_prior_list.md)
  : Prior for the probit regression coefficients of missBART1

## Simulating data

- [`sim_data_friedman()`](https://andrewcparnell.github.io/missBART/reference/sim_data_friedman.md)
  : Simulate
- [`sim_data_trees()`](https://andrewcparnell.github.io/missBART/reference/sim_data_trees.md)
  : Simulating data from trees
- [`sim_missing_probit()`](https://andrewcparnell.github.io/missBART/reference/sim_missing_probit.md)
  : Simulating missing values from a probit regression model
- [`sim_missing_trees()`](https://andrewcparnell.github.io/missBART/reference/sim_missing_trees.md)
  : Simulate missing data from a BART model
- [`propose_tree()`](https://andrewcparnell.github.io/missBART/reference/propose_tree.md)
  : A function for proposing an updated tree structure

## Utilities

- [`scale_bart()`](https://andrewcparnell.github.io/missBART/reference/scale_bart.md)
  [`unscale()`](https://andrewcparnell.github.io/missBART/reference/scale_bart.md)
  : Scaling responses to and from the BART scale
- [`rMVN()`](https://andrewcparnell.github.io/missBART/reference/rMVN.md)
  : Sampling once from a multivariate normal distribution with mean mu
  and precision Q
- [`rMVNc()`](https://andrewcparnell.github.io/missBART/reference/rMVNc.md)
  : Sampling from a multivariate normal distribution with precision Q
  and mean equal to the solution of Q mu = b
- [`multi_rMVN()`](https://andrewcparnell.github.io/missBART/reference/multi_rMVN.md)
  : Sampling multiple times from a multivariate normal distribution
- [`matrnorm()`](https://andrewcparnell.github.io/missBART/reference/matrnorm.md)
  : Sampling from a matrix normal distribution MV(mu, U, V)
- [`missBART-package`](https://andrewcparnell.github.io/missBART/reference/missBART-package.md)
  : missBART: joint models for non-ignorable missing responses using
  BART
