# Quick start

``` r
library(missBART)
```

**missBART** fits selection models to multivariate responses whose
missing values may be *missing not at random* (MNAR): whether a value is
recorded can depend on the value itself. The responses follow a
multivariate BART model. The missingness indicators follow either

- a multivariate probit regression on the covariates and the responses
  ([`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)),
  or
- a multivariate probit BART model on the covariates and the responses
  ([`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)).

The missing responses are imputed while the model is fitted.
[`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md)
fits the multivariate BART model on its own to complete responses. The
models are described in Goh, Soh, Parnell and Murphy (2024),
[arXiv:2412.14946](https://arxiv.org/abs/2412.14946).

## Simulated data

We simulate two correlated responses from five covariates, and then
delete values with a probability that increases with the value itself.

``` r
set.seed(2024)
n <- 400
x <- matrix(runif(n * 5), ncol = 5)
f <- cbind(10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5],
           5 * x[, 1] - 8 * x[, 2]^2 + 3 * x[, 3])
y_true <- f + matrix(rnorm(2 * n), ncol = 2) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2))

y <- y_true
y[runif(n) < pnorm(-1 + 1.2 * scale(y_true[, 1])), 1] <- NA
y[runif(n) < pnorm(-1 + 1.2 * scale(y_true[, 2])), 2] <- NA
colMeans(is.na(y))
#> [1] 0.2375 0.2625
```

The observed values are a biased sample: the missing values are larger
on average than the observed ones.

``` r
rbind(observed = colMeans(y, na.rm = TRUE),
      missing = c(mean(y_true[is.na(y[, 1]), 1]), mean(y_true[is.na(y[, 2]), 2])))
#>              [,1]      [,2]
#> observed 12.77439 0.2265791
#> missing  19.15552 4.3086642
```

## missBART1: probit regression for the missingness

The chains here are short so that the vignette builds quickly; use the
defaults, or longer chains, for real analyses.

``` r
fit1 <- missBART1(x, y, n_trees = 50, burn = 500, iters = 500, thin = 1,
                  mice_impute = FALSE, show_progress = FALSE)
fit1
#> missBART1 fit: 400 observations, 2 responses
#> Missing responses: 200 of 800 (25%)
#> Posterior samples: 500 (burn-in 500, thinning 1)
#> Tree acceptance rate: regression 0.087
```

[`imputed_values()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md)
summarises the posterior of every missing response, and
[`completed_data()`](https://andrewcparnell.github.io/missBART/reference/imputed_values.md)
fills in their posterior means.

``` r
imp1 <- imputed_values(fit1)
head(imp1)
#>   row column     mean       sd    lower    upper
#> 1   1     Y1 15.85965 1.105576 13.69506 18.07238
#> 2  12     Y1 21.83528 1.047334 19.81081 23.87049
#> 3  23     Y1 22.69528 1.075083 20.80349 24.86201
#> 4  26     Y1 17.82940 1.080126 15.62825 20.11384
#> 5  35     Y1 20.57692 1.275867 18.12664 23.08140
#> 6  36     Y1 20.37656 1.091618 18.29674 22.58855
truth <- y_true[fit1$missing_index]
plot(truth, imp1$mean, xlab = "True value", ylab = "Posterior mean imputation",
     col = factor(imp1$column), pch = 19, cex = 0.6)
abline(0, 1)
```

![](missBART_files/figure-html/unnamed-chunk-5-1.png)

The probit coefficients describe the missingness mechanism. Each column
refers to the missingness of one response, and a positive coefficient
makes that response more likely to be *missing*. Posterior intervals for
the rows `Y1` and `Y2` that exclude zero indicate MNAR data. (The
coefficients refer to standardised covariates and responses; see
[`?missBART1`](https://andrewcparnell.github.io/missBART/reference/missBART.md).)

``` r
round(coef(fit1), 2)
#> , , mean
#> 
#>              Y1    Y2
#> Intercept -1.71 -1.55
#> X1        -0.72  0.04
#> X2         0.54 -0.44
#> X3        -0.48  0.31
#> X4        -0.35 -0.09
#> X5         0.08 -0.28
#> Y1         5.51  0.47
#> Y2         1.62  5.28
#> 
#> , , lower
#> 
#>              Y1    Y2
#> Intercept -2.08 -1.98
#> X1        -1.61 -1.15
#> X2        -0.73 -1.88
#> X3        -1.14 -0.39
#> X4        -1.18 -0.97
#> X5        -0.52 -0.95
#> Y1         3.96 -0.78
#> Y2        -0.12  2.71
#> 
#> , , upper
#> 
#>              Y1    Y2
#> Intercept -1.37 -1.19
#> X1         0.16  0.99
#> X2         1.97  0.92
#> X3         0.14  1.00
#> X4         0.48  0.72
#> X5         0.71  0.38
#> Y1         7.12  1.87
#> Y2         3.57  7.84
```

## missBART2: probit BART for the missingness

``` r
fit2 <- missBART2(x, y, n_reg_trees = 50, n_class_trees = 20, burn = 500, iters = 500,
                  thin = 1, mice_impute = FALSE, show_progress = FALSE)
fit2
#> missBART2 fit: 400 observations, 2 responses
#> Missing responses: 200 of 800 (25%)
#> Posterior samples: 500 (burn-in 500, thinning 1)
#> Tree acceptance rate: regression 0.089, missingness 0.197
#> Acceptance rate for missing responses: 0.19
```

There are no coefficients to inspect here. Instead,
[`varimp()`](https://andrewcparnell.github.io/missBART/reference/varimp.md)
reports how often each predictor is used in the splitting rules of the
missingness trees. Splits concentrated on the responses point towards
MNAR data.

``` r
barplot(varimp(fit2), ylab = "Proportion of splitting rules")
```

![](missBART_files/figure-html/unnamed-chunk-8-1.png)

## Comparison with ignoring the missingness

Fitting multivariate BART to the complete cases and predicting the
incomplete rows ignores the reason the values are missing.

``` r
cc <- complete.cases(y)
fit0 <- mvBART(x[cc, ], y[cc, ], x_predict = x[!cc, ], n_trees = 50, burn = 500,
               iters = 500, thin = 1, show_progress = FALSE)
pred0 <- matrix(NA, n, 2)
pred0[!cc, ] <- predict(fit0)

rmse <- function(a, b) sqrt(mean((a - b)^2))
c(complete_cases = rmse(pred0[fit1$missing_index], truth),
  missBART1 = rmse(imp1$mean, truth),
  missBART2 = rmse(imputed_values(fit2)$mean, truth))
#> complete_cases      missBART1      missBART2 
#>       1.666967       1.153195       1.132137
```

## What is in a fit

A fit is a list of posterior draws, one list element per kept iteration;
see
[`?missBART1`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
for all components. For example, `y_post` holds the BART mean function,
`y_pred` posterior predictive draws, `y_impute` the imputed values
(ordered as `which(is.na(y))`), `sigma_post` the residual covariance
matrix and, if `x_predict` is supplied, `new_y_hat` and `new_y_post` the
mean function and predictive draws for new covariates.

``` r
round(Reduce("+", fit1$sigma_post) / length(fit1$sigma_post), 2)
#>      Y1   Y2
#> Y1 0.92 0.34
#> Y2 0.34 0.78
```
