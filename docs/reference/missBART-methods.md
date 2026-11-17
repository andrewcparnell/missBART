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
predict(object, type = c("mean", "draws"), ...)

# S3 method for class 'missBART'
coef(object, probs = c(0.025, 0.975), ...)
```

## Arguments

- x, object:

  an object of class `"missBART"` returned by
  [`missBART1()`](missBART.md), [`missBART2()`](missBART.md) or
  [`mvBART()`](mvBART.md).

- ...:

  not used.

- type:

  `"mean"` for the posterior mean of the BART mean function at
  `x_predict`, or `"draws"` for the list of posterior predictive draws.

- probs:

  two probabilities giving the posterior interval to report.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
posterior mean of the BART mean function for the training rows and
[`predict()`](https://rdrr.io/r/stats/predict.html) the same at the
`x_predict` supplied when fitting.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the posterior mean
of the probit coefficients of a [`missBART1()`](missBART.md) fit
together with their posterior intervals (an array with one slice per
summary).
