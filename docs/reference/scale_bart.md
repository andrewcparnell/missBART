# Scaling responses to and from the BART scale

`scale_bart()` maps data linearly so that `min` and `max` become -0.5
and 0.5; `unscale()` inverts the map.

## Usage

``` r
scale_bart(data, min = NULL, max = NULL)

unscale(scaled_val, min, max)
```

## Arguments

- data, scaled_val:

  numeric vector or matrix (columns are scaled separately).

- min, max:

  the minimum and maximum (one per column for a matrix). `scale_bart()`
  computes them from `data` if they are not given.

## Value

The transformed data.

## Examples

``` r
y <- matrix(rnorm(10), ncol = 2)
s <- scale_bart(y)
all.equal(unscale(s, apply(y, 2, min), apply(y, 2, max)), y)
#> [1] TRUE
```
