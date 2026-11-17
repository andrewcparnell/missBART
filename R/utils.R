#' Scaling responses to and from the BART scale
#'
#' `scale_bart()` maps data linearly so that `min` and `max` become -0.5 and
#' 0.5; `unscale()` inverts the map.
#'
#' @param data,scaled_val numeric vector or matrix (columns are scaled
#'   separately).
#' @param min,max the minimum and maximum (one per column for a matrix).
#'   `scale_bart()` computes them from `data` if they are not given.
#'
#' @return The transformed data.
#' @export
#'
#' @examples
#' y <- matrix(rnorm(10), ncol = 2)
#' s <- scale_bart(y)
#' all.equal(unscale(s, apply(y, 2, min), apply(y, 2, max)), y)
scale_bart <- function(data, min = NULL, max = NULL) {
  if(is.matrix(data)) {
    if(is.null(min)) min <- apply(data, 2, base::min, na.rm = TRUE)
    if(is.null(max)) max <- apply(data, 2, base::max, na.rm = TRUE)
    sweep(sweep(data, 2, min, "-"), 2, max - min, "/") - 0.5
  } else {
    if(is.null(min)) min <- base::min(data, na.rm = TRUE)
    if(is.null(max)) max <- base::max(data, na.rm = TRUE)
    (data - min) / (max - min) - 0.5
  }
}

#' @rdname scale_bart
#' @export
unscale <- function(scaled_val, min, max) {
  if(is.matrix(scaled_val)) {
    sweep(sweep(scaled_val + 0.5, 2, max - min, "*"), 2, min, "+")
  } else {
    (scaled_val + 0.5) * (max - min) + min
  }
}
