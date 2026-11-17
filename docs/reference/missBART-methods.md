# Methods for missBART fits

Methods for missBART fits

## Usage

``` r
# S3 method for class 'missBART'
print(x, ...)

# S3 method for class 'missBART'
summary(object, ...)

# S3 method for class 'missBART'
fitted(object, ...)

# S3 method for class 'missBART'
predict(object, newdata, type = c("mean", "draws"), ...)

# S3 method for class 'missBART'
plot(x, which = c("imputed", "fit"), ...)

# S3 method for class 'missBART'
coef(object, probs = c(0.025, 0.975), ...)
```

## Arguments

- x, object:

  an object of class `"missBART"` returned by
  [`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md),
  [`missBART2()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
  or
  [`mvBART()`](https://andrewcparnell.github.io/missBART/reference/mvBART.md).

- ...:

  not used.

- newdata:

  covariates at which to predict. If missing, the predictions at the
  `x_predict` supplied when fitting are returned. Predictions at
  `newdata` need a fit with `keep_trees = TRUE`.

- type:

  `"mean"` for the posterior mean of the BART mean function, or
  `"draws"` for posterior draws: the list of posterior predictive draws
  at `x_predict`, or an array (rows, responses, samples) of draws of the
  mean function at `newdata`.

- which:

  for [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
  `"imputed"` compares the distribution of the imputed values (posterior
  means) with that of the observed values of each response, and `"fit"`
  plots the fitted against the observed values with the squared
  correlation in the panel title.

- probs:

  two probabilities giving the posterior interval to report.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
posterior mean of the BART mean function for the training rows and
[`predict()`](https://rdrr.io/r/stats/predict.html) the same for new
covariates. [`coef()`](https://rdrr.io/r/stats/coef.html) returns the
posterior mean of the probit coefficients of a
[`missBART1()`](https://andrewcparnell.github.io/missBART/reference/missBART.md)
fit together with their posterior intervals (an array with one slice per
summary). A positive coefficient means that larger values of the
predictor make the response more likely to be missing.
