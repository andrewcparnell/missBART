post_mean <- function(draws) Reduce("+", draws) / length(draws)

#' Methods for missBART fits
#'
#' @param x,object an object of class `"missBART"` returned by [missBART1()],
#'   [missBART2()] or [mvBART()].
#' @param type `"mean"` for the posterior mean of the BART mean function at
#'   `x_predict`, or `"draws"` for the list of posterior predictive draws.
#' @param probs two probabilities giving the posterior interval to report.
#' @param ... not used.
#'
#' @return `fitted()` returns the posterior mean of the BART mean function for
#'   the training rows and `predict()` the same at the `x_predict` supplied
#'   when fitting. `coef()` returns the posterior mean of the probit
#'   coefficients of a `missBART1()` fit together with their posterior
#'   intervals (an array with one slice per summary).
#' @name missBART-methods
NULL

#' @rdname missBART-methods
#' @export
print.missBART <- function(x, ...) {
  n_mis <- length(x$missing_index)
  cat(x$model, " fit: ", nrow(x$y), " observations, ", ncol(x$y), " response",
      if(ncol(x$y) > 1) "s", "\n", sep = "")
  if(x$model != "mvBART") {
    cat("Missing responses: ", n_mis, " of ", length(x$y), " (",
        format(100 * n_mis / length(x$y), digits = 3), "%)\n", sep = "")
  }
  cat("Posterior samples: ", x$iters, " (burn-in ", x$burn, ", thinning ", x$thin, ")\n", sep = "")
  cat("Tree acceptance rate: ",
      paste0(names(x$tree_accept), " ", format(x$tree_accept, digits = 2), collapse = ", "), "\n", sep = "")
  if(!is.null(x$y_miss_accept)) {
    cat("Acceptance rate for missing responses: ", format(mean(x$y_miss_accept), digits = 2), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname missBART-methods
#' @export
summary.missBART <- function(object, ...) {
  print(object)
  obs <- !is.na(object$y)
  fit <- fitted(object)
  rmse <- vapply(seq_len(ncol(object$y)), function(j) sqrt(mean((object$y[obs[, j], j] - fit[obs[, j], j])^2)), numeric(1))
  names(rmse) <- colnames(object$y)
  cat("\nIn-sample RMSE on observed responses:\n")
  print(rmse, digits = 4)
  cat("\nPosterior mean of the residual covariance:\n")
  print(post_mean(object$sigma_post), digits = 4)
  if(object$model == "missBART1") {
    cat("\nPosterior mean of the probit coefficients (rows: predictors, columns: missingness of each response):\n")
    print(post_mean(object$B_post), digits = 3)
  }
  if(object$model == "missBART2") {
    cat("\nVariable importance in the missingness trees:\n")
    print(varimp(object), digits = 3)
  }
  invisible(object)
}

#' @rdname missBART-methods
#' @export
fitted.missBART <- function(object, ...) post_mean(object$y_post)

#' @rdname missBART-methods
#' @export
predict.missBART <- function(object, type = c("mean", "draws"), ...) {
  type <- match.arg(type)
  if(!is.list(object$new_y_hat)) stop("No predictions stored: supply x_predict when fitting the model")
  if(type == "mean") post_mean(object$new_y_hat) else object$new_y_post
}

#' @rdname missBART-methods
#' @export
coef.missBART <- function(object, probs = c(0.025, 0.975), ...) {
  if(object$model != "missBART1") stop("Only missBART1 fits have probit coefficients")
  draws <- simplify2array(object$B_post)
  q <- apply(draws, c(1, 2), stats::quantile, probs = probs)
  out <- array(c(apply(draws, c(1, 2), mean), q[1, , ], q[2, , ]), dim = c(dim(draws)[1:2], 3),
               dimnames = c(dimnames(draws)[1:2], list(c("mean", "lower", "upper"))))
  out
}

#' Summaries of the imputed responses
#'
#' @param object an object of class `"missBART"`.
#' @param probs two probabilities giving the posterior interval to report.
#'
#' @return `imputed_values()` returns a data frame with one row per missing
#'   response: its row and column in `y`, and the posterior mean, standard
#'   deviation and interval of its imputations. `completed_data()` returns `y`
#'   with the missing values replaced by their posterior means.
#' @export
imputed_values <- function(object, probs = c(0.025, 0.975)) {
  stopifnot(inherits(object, "missBART"))
  idx <- object$missing_index
  if(length(idx) == 0) stop("The responses have no missing values")
  draws <- do.call(rbind, object$y_impute)
  q <- apply(draws, 2, stats::quantile, probs = probs)
  n <- nrow(object$y)
  data.frame(row = (idx - 1) %% n + 1, column = colnames(object$y)[(idx - 1) %/% n + 1],
             mean = colMeans(draws), sd = apply(draws, 2, stats::sd), lower = q[1, ], upper = q[2, ])
}

#' @rdname imputed_values
#' @export
completed_data <- function(object) {
  stopifnot(inherits(object, "missBART"))
  y <- object$y
  if(length(object$missing_index) > 0) y[object$missing_index] <- colMeans(do.call(rbind, object$y_impute))
  y
}

#' Variable importance
#'
#' The proportion of splitting rules that use each predictor, averaged over
#' the posterior samples. For the missingness trees of [missBART2()], splits
#' concentrated on the responses point towards data that are missing not at
#' random.
#'
#' @param object an object of class `"missBART"`.
#' @param trees which set of trees: `"missingness"` (only for [missBART2()]
#'   fits) or `"regression"`.
#'
#' @return A named vector of proportions.
#' @export
varimp <- function(object, trees = c("missingness", "regression")) {
  stopifnot(inherits(object, "missBART"))
  trees <- match.arg(trees)
  if(trees == "missingness" && object$model != "missBART2") trees <- "regression"
  counts <- if(trees == "missingness") object$var_imp else object$var_imp_reg
  totals <- rowSums(counts)
  colMeans(counts[totals > 0, , drop = FALSE] / totals[totals > 0])
}
