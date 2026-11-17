# A function for proposing an updated tree structure

A function for proposing an updated tree structure

## Usage

``` r
propose_tree(
  df,
  x,
  min_node,
  max_attempt = 10,
  i,
  probit = FALSE,
  miss_row = NA,
  cat_list = NA
)
```

## Arguments

- df:

  dataframe

- x:

  covariates

- min_node:

  minimum number of observations that should fall into a single terminal
  node.

- max_attempt:

  maximum number of attempts to find a suitable tree

- i:

  iteration number

- probit:

  logical. Indicates whether tree is for continuous or binary outcomes.

- miss_row:

  index of missing rows

- cat_list:

  list of used categorical subset

## Value

A data.frame containing details of a new tree structure

## Examples

``` r
# Create a data.frame representing the root node of a tree
df <- data.frame(matrix(ncol = 8, nrow = 1))
colnames(df) <- c("parent", "lower", "upper", "split_variable",
                  "split_value", "depth", "direction", "NA_direction")
df[1,] <- c(0,0,1,0,1,0,0, NA)
x <- matrix(stats::runif(9), ncol=3)
propose_tree(df, x, min_node = 1, max_attempt = 1, i = 1)
#> $new_df
#>   parent     lower    upper split_variable split_value depth direction
#> 1      0 0.0000000 1.000000              0     1.00000     0         0
#> 2      1 0.1010004 0.925852              3     0.90471     1         0
#> 3      1 0.1010004 0.925852              3     0.90471     1         1
#>   NA_direction
#> 1           NA
#> 2           NA
#> 3           NA
#> 
#> $MOVE
#> [1] "GROW"
#> 
#> $change_points
#> [1] 2 1 1
#> 
#> $decent_tree
#> [1] TRUE
#> 
```
