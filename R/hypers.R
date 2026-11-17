#' BART prior parameters
#'
#' @param mu0 hyperparameter
#' @param kappa hyperparameter
#' @param alpha hyperparameter
#' @param V hyperparameter
#' @param df hyperparameter
#' @param q hyperparameter
#' @param tau_b not used; kept for compatibility with earlier versions
#' @param ... Catches unused arguments
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples hypers_list(mu0 = 5, kappa = 2, alpha = 1, V = 2)
hypers_list <- function(mu0 = 0, kappa = NULL, alpha = NULL, V = NULL, df = 10, q = 0.75, tau_b = 100, ...) {
  return(list(mu0 = mu0, kappa = kappa, alpha = alpha, V = V, df = df, q = q, tau_b = tau_b))
}

#' Tree prior parameters
#'
#' @param prior_alpha hyperparameter
#' @param prior_beta hyperparameter
#' @param min_node minimum number of observations that should fall into a single terminal node.
#' @param max_attempt maximum number of attempts to find a suitable tree
#' @param ... Catches unused arguments
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples tree_list(prior_alpha = 0.95, prior_beta = 3)
tree_list <- function(prior_alpha = 0.95, prior_beta = 2, min_node = 1, max_attempt = 1, ...) {
  if(prior_alpha >= 1 || prior_alpha <= 0) stop("prior_alpha must be between (0,1)")
  if(prior_beta < 0) stop("prior_beta must be positive")
  return(list(prior_alpha = prior_alpha, prior_beta = prior_beta, min_node = min_node, max_attempt = max_attempt))
}

#' Prior for the probit regression coefficients of missBART1
#'
#' In [missBART1()] the probit coefficients have prior precisions
#' `tau_0` (intercept), `tau_X` (covariates) and `tau_Y` (responses) with
#' `Gamma(a, b)` priors (shape `a`, rate `b`). A prior that favours large
#' `tau_Y` shrinks the coefficients of the responses towards zero, i.e. favours
#' data that are missing at random.
#'
#' @param a0,b0 shape and rate for the intercept precision.
#' @param ax,bx shape and rate for the covariate precision. The defaults are
#'   `ax = q + 1` and `bx = p + q + 1`, with `q` covariates and `p` responses.
#' @param ay,by shape and rate for the response precision. The defaults are
#'   `ay = 1` and `by = 1`.
#'
#' @details The defaults are those of the package code before version 0.2.0.
#'   Section 5.1 of the paper instead states `(ax, bx) = (1 + q, 1)` and
#'   `(ay, by) = (1 + p + q, 1)`; use
#'   `probit_prior_list(bx = 1, ay = 1 + p + q)` for that setting.
#'
#' @return list of prior parameters
#' @export
#'
#' @examples probit_prior_list(ay = 10)
probit_prior_list <- function(a0 = 2, b0 = 1, ax = NULL, bx = NULL, ay = NULL, by = NULL) {
  vals <- c(a0, b0, ax, bx, ay, by)
  if(any(!is.finite(vals)) || any(vals <= 0)) stop("All prior parameters must be positive")
  list(a0 = a0, b0 = b0, ax = ax, bx = bx, ay = ay, by = by)
}
