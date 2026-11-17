# Tree prior parameters

Tree prior parameters

## Usage

``` r
tree_list(
  prior_alpha = 0.95,
  prior_beta = 2,
  min_node = 1,
  max_attempt = 1,
  ...
)
```

## Arguments

- prior_alpha:

  hyperparameter

- prior_beta:

  hyperparameter

- min_node:

  minimum number of observations that should fall into a single terminal
  node.

- max_attempt:

  maximum number of attempts to find a suitable tree

- ...:

  Catches unused arguments

## Value

list of hyperparameters

## Examples

``` r
tree_list(prior_alpha = 0.95, prior_beta = 3)
#> $prior_alpha
#> [1] 0.95
#> 
#> $prior_beta
#> [1] 3
#> 
#> $min_node
#> [1] 1
#> 
#> $max_attempt
#> [1] 1
#> 
```
