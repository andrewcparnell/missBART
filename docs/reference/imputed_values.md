# Summaries of the imputed responses

Summaries of the imputed responses

## Usage

``` r
imputed_values(object, probs = c(0.025, 0.975))

completed_data(object)
```

## Arguments

- object:

  an object of class `"missBART"`.

- probs:

  two probabilities giving the posterior interval to report.

## Value

`imputed_values()` returns a data frame with one row per missing
response: its row and column in `y`, and the posterior mean, standard
deviation and interval of its imputations. `completed_data()` returns
`y` with the missing values replaced by their posterior means.
