# missBART

**Website and documentation: <https://andrewcparnell.github.io/missBART/>**

Joint models for multivariate responses with non-ignorable missing values,
built on Bayesian additive regression trees (BART).

The responses follow a multivariate BART model. The missingness indicators
follow a model in which the responses themselves are predictors, so that data
missing completely at random, missing at random and missing not at random
(MNAR) can all be handled, and the missing responses are imputed while the
model is fitted:

| Function      | Model for the missingness                                  |
|---------------|------------------------------------------------------------|
| `missBART1()` | multivariate probit regression on `(1, x, y)`              |
| `missBART2()` | multivariate probit BART on `(x, y)`                       |
| `mvBART()`    | none: multivariate BART for complete responses             |

The methods are described in

> Goh, Y. C., Soh, W. K., Parnell, A. C. and Murphy, K. (2024). Joint models
> for handling non-ignorable missing data using Bayesian additive regression
> trees: application to leaf photosynthetic traits data.
> [arXiv:2412.14946](https://arxiv.org/abs/2412.14946)

## Installation

The samplers are written in C++ (Rcpp and RcppArmadillo), so a C++ toolchain
is needed to install from source (Rtools on Windows, Xcode command line tools
on macOS).

```r
# install.packages("remotes")
remotes::install_github("yongchengoh/missBART")
```

## Example

```r
library(missBART)

set.seed(1)
n <- 500
x <- matrix(runif(n * 5), ncol = 5)
y_true <- cbind(10 * sin(pi * x[, 1] * x[, 2]) + 10 * x[, 4], 5 * x[, 1] - 8 * x[, 2]^2) +
  matrix(rnorm(2 * n), ncol = 2)

# large values are more likely to be missing
y <- y_true
y[runif(n) < pnorm(-1 + scale(y_true[, 1])), 1] <- NA
y[runif(n) < pnorm(-1 + scale(y_true[, 2])), 2] <- NA

fit1 <- missBART1(x, y)                      # probit regression missingness model
fit2 <- missBART2(x, y, n_class_trees = 20)  # probit BART missingness model

fit1
head(imputed_values(fit1))   # posterior summaries of the missing responses
y_complete <- completed_data(fit1)
coef(fit1)                   # probit coefficients: positive = more likely to be missing
varimp(fit2)                 # which predictors the missingness trees split on
```

`vignette("missBART")` has a longer walk-through, and the [package website](https://andrewcparnell.github.io/missBART/) has further articles, including a re-run of the simulation studies in the paper. Scripts for the
simulations and the application in the paper are in [`Examples/`](Examples).

## Speed

Version 0.2.0 moved the samplers from R to C++. Seconds per 1,000 MCMC
iterations on simulated MNAR data with five covariates, 100 regression trees
and (for `missBART2`) 20 missingness trees, single core:

| n    | responses | model     | 0.1.0 (R) | 0.2.0 (C++) | speed-up |
|------|-----------|-----------|----------:|------------:|---------:|
| 500  | 1         | missBART1 |        95 |        0.57 |     170x |
| 500  | 1         | missBART2 |       110 |        0.63 |     170x |
| 500  | 1         | mvBART    |        86 |        0.48 |     180x |
| 2000 | 2         | missBART1 |      1043 |        3.0  |     350x |
| 2000 | 2         | missBART2 |       248 |        2.9  |      90x |
| 2000 | 2         | mvBART    |       186 |        2.1  |      90x |
| 2000 | 5         | missBART1 |      1553 |        5.4  |     290x |
| 2000 | 5         | missBART2 |       389 |        5.0  |      80x |
| 2000 | 5         | mvBART    |       304 |        3.5  |      90x |

The C++ timings are averages over 5,000 iterations. The R timings are averages
over the first 50 iterations only, when the trees are still small, so they
flatter the R code. At these rates the 10,000 iterations used in the paper for
`n = 2000` and five responses take about 55 seconds with `missBART1`, instead
of about 4.3 hours. [`Examples/benchmark.R`](Examples/benchmark.R) reproduces
the table.

Running both versions for the same 2,000 iterations on one simulated data set
(`n = 500`, two responses) gave:

| model     | time 0.1.0 | time 0.2.0 | imputation RMSE 0.1.0 | imputation RMSE 0.2.0 |
|-----------|-----------:|-----------:|----------------------:|----------------------:|
| missBART1 |      621 s |      1.7 s |                  1.24 |                  1.17 |
| missBART2 |      246 s |      1.5 s |                  1.16 |                  1.11 |
| mvBART    |      196 s |      1.1 s |                     – |                     – |

The two versions do not produce identical draws for the same seed, and a few
steps of the sampler were changed; see [`NEWS.md`](NEWS.md).
