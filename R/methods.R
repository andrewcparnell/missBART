post_mean <- function(draws) Reduce("+", draws) / length(draws)

#' Methods for missBART fits
#'
#' @param x,object an object of class `"missBART"` returned by [missBART1()],
#'   [missBART2()] or [mvBART()].
#' @param newdata covariates at which to predict. If missing, the predictions
#'   at the `x_predict` supplied when fitting are returned. Predictions at
#'   `newdata` need a fit with `keep_trees = TRUE`.
#' @param type `"mean"` for the posterior mean of the BART mean function, or
#'   `"draws"` for posterior draws: the list of posterior predictive draws at
#'   `x_predict`, or an array (rows, responses, samples) of draws of the mean
#'   function at `newdata`.
#' @param which for `plot()`: `"imputed"` compares the distribution of the
#'   imputed values (posterior means) with that of the observed values of each
#'   response, and `"fit"` plots the fitted against the observed values with
#'   the squared correlation in the panel title.
#' @param probs two probabilities giving the posterior interval to report.
#' @param ... not used.
#'
#' @return `fitted()` returns the posterior mean of the BART mean function for
#'   the training rows and `predict()` the same for new covariates. `coef()` returns the posterior mean of the probit
#'   coefficients of a `missBART1()` fit together with their posterior
#'   intervals (an array with one slice per summary). A positive coefficient
#'   means that larger values of the predictor make the response more likely
#'   to be missing.
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
  if(x$n_chains > 1) {
    rhat <- max(chain_diagnostics(x)$rhat, na.rm = TRUE)
    cat("Chains: ", x$n_chains, ", largest R-hat ", format(rhat, digits = 3), "\n", sep = "")
    if(rhat > 1.1) cat("The chains disagree: see chain_diagnostics() before interpreting this fit.\n")
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
predict.missBART <- function(object, newdata, type = c("mean", "draws"), ...) {
  type <- match.arg(type)
  if(missing(newdata)) {
    if(!is.list(object$new_y_hat)) stop("No predictions stored: supply x_predict when fitting the model, or newdata here")
    return(if(type == "mean") post_mean(object$new_y_hat) else object$new_y_post)
  }
  need_trees(object)
  draws <- .forest_predict_cpp(object$reg_trees, prep_newx(object, newdata), ncol(object$y))
  for(j in seq_len(ncol(object$y))) draws[, j, ] <- draws[, j, ] * object$rng_y[j] + object$center_y[j]
  dimnames(draws) <- list(NULL, colnames(object$y), NULL)
  if(type == "mean") apply(draws, c(1, 2), mean) else draws
}

#' @rdname missBART-methods
#' @export
plot.missBART <- function(x, which = c("imputed", "fit"), ...) {
  which <- match.arg(which)
  p <- ncol(x$y)
  old <- graphics::par(mfrow = grDevices::n2mfrow(p), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old))
  if(which == "fit") {
    fit <- fitted(x)
    for(j in seq_len(p)) {
      ok <- !is.na(x$y[, j])
      graphics::plot(x$y[ok, j], fit[ok, j], xlab = "Observed", ylab = "Fitted", pch = 16, cex = 0.5, col = "#1A9993",
                     main = sprintf("%s (R squared = %.2f)", colnames(x$y)[j], stats::cor(x$y[ok, j], fit[ok, j])^2), ...)
      graphics::abline(0, 1)
    }
  } else {
    if(length(x$missing_index) == 0) stop("The responses have no missing values")
    imp <- imputed_values(x)
    for(j in seq_len(p)) {
      obs <- stats::density(x$y[!is.na(x$y[, j]), j])
      vals <- imp$mean[imp$column == colnames(x$y)[j]]
      den <- if(length(vals) > 1) stats::density(vals) else NULL
      graphics::plot(obs, main = colnames(x$y)[j], xlab = "", col = "#1A9993", lwd = 2,
                     xlim = range(obs$x, den$x, vals), ylim = range(0, obs$y, den$y), ...)
      if(!is.null(den)) graphics::lines(den, col = "#FD7446", lwd = 2)
      graphics::rug(vals, col = "#FD7446")
      if(j == 1) graphics::legend("topright", c("Observed", "Imputed"), col = c("#1A9993", "#FD7446"), lwd = 2, bty = "n")
    }
  }
  invisible(x)
}

need_trees <- function(object) {
  if(is.null(object$reg_trees)) stop("This needs the trees: refit the model with keep_trees = TRUE")
}

# Brings new covariates to the form used by the sampler.
prep_newx <- function(object, newdata) {
  newdata <- as_numeric_matrix(newdata, "newdata")
  q0 <- ncol(object$x) - length(object$na_cols)
  if(ncol(newdata) != q0) stop("newdata must have ", q0, " columns")
  if(object$model == "missBART1") {
    if(anyNA(newdata)) stop("missBART1 needs complete covariates")
    if(!is.null(object$min_x)) newdata <- range01(newdata, object$min_x, object$rng_x)
  } else if(length(object$na_cols) > 0) {
    newdata <- cbind(newdata, 1 - is.na(newdata[, object$na_cols, drop = FALSE]))
  }
  newdata
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

#' Partial dependence and individual conditional expectation curves
#'
#' Shows how the predictions of a fitted model change as one predictor varies
#' and the others are held at their observed values. `type = "response"`
#' gives the mean of each response as a covariate varies. `type =
#' "missingness"` gives the probability that each response is missing as a
#' covariate or a response varies.
#'
#' @param object an object of class `"missBART"`. Except for the missingness
#'   model of a [missBART1()] fit, it must have been fitted with
#'   `keep_trees = TRUE`.
#' @param variable name or column number of the predictor to vary: a column of
#'   `x`, or for `type = "missingness"` also a response.
#' @param type which model to use, see above.
#' @param grid values of the predictor at which to evaluate the curves, on the
#'   scale of the data. The default is `n_grid` values spanning its observed
#'   range.
#' @param n_grid number of grid values if `grid` is not given.
#' @param n_ice number of observations (chosen at random) for which individual
#'   curves are computed and over which the partial dependence is averaged.
#' @param draws number of posterior samples to use (evenly spaced).
#' @param probs two probabilities giving the posterior interval of the partial
#'   dependence curve.
#' @param x an object returned by `partial_dependence()`.
#' @param ... passed to [graphics::matplot()].
#'
#' @details For the missingness model, missing responses are replaced by their
#'   posterior means before the curves are computed. Probabilities are for
#'   each response on its own, given the predictors.
#'
#' @return An object of class `"missBART_pdp"`: a list with the `grid`, the
#'   partial dependence curve `pdp` (an array: grid values, responses, and
#'   mean/lower/upper) and the individual curves `ice` (observations, grid
#'   values, responses; posterior means). The `plot()` method draws one panel
#'   per response.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(runif(300), ncol = 3)
#' y <- cbind(sin(2 * pi * x[, 1]), x[, 2]) + matrix(rnorm(200, sd = 0.2), ncol = 2)
#' y[runif(100) < pnorm(2 * y[, 1] - 1), 1] <- NA
#' fit <- missBART2(x, y, n_reg_trees = 10, n_class_trees = 5, burn = 50, iters = 50,
#'                  thin = 1, mice_impute = FALSE, show_progress = FALSE, keep_trees = TRUE)
#' plot(partial_dependence(fit, "X1"))
#' plot(partial_dependence(fit, "Y1", type = "missingness"))
partial_dependence <- function(object, variable, type = c("response", "missingness"), grid = NULL, n_grid = 20,
                               n_ice = 200, draws = 200, probs = c(0.025, 0.975)) {
  stopifnot(inherits(object, "missBART"))
  type <- match.arg(type)
  p <- ncol(object$y)
  y_names <- colnames(object$y)
  x_names <- colnames(object$x)
  if(type == "missingness" && object$model == "mvBART") stop("mvBART fits have no missingness model")
  # predictors of the chosen model, on the scale of the data
  y_complete <- completed_data(object)
  x_data <- object$x
  if(!is.null(object$min_x)) x_data <- sweep(sweep(x_data + 0.5, 2, object$rng_x, "*"), 2, object$min_x, "+")
  preds <- if(type == "response") x_names else c(if(object$include_x) x_names, if(object$include_y) y_names)
  if(is.numeric(variable)) variable <- preds[variable]
  if(length(variable) != 1 || !variable %in% preds) stop("variable must be one of: ", paste(preds, collapse = ", "))
  is_y <- type == "missingness" && variable %in% y_names && !(variable %in% x_names && object$include_x)
  values <- if(is_y) y_complete[, variable] else x_data[, variable]
  if(is.null(grid)) grid <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE), length.out = n_grid)
  # the same values on the scale used by the sampler
  to_model <- function(v) {
    if(is_y) {
      (v - object$center_y[variable]) / object$rng_y[variable]
    } else if(!is.null(object$min_x)) {
      (v - object$min_x[match(variable, x_names)]) / object$rng_x[match(variable, x_names)] - 0.5
    } else v
  }
  rows <- if(nrow(object$x) > n_ice) sort(sample.int(nrow(object$x), n_ice)) else seq_len(nrow(object$x))
  n_draws <- length(object$y_post)
  keep <- unique(round(seq(1, n_draws, length.out = min(draws, n_draws))))
  y_model <- sweep(sweep(y_complete, 2, object$center_y, "-"), 2, object$rng_y, "/")

  if(type == "response") {
    need_trees(object)
    res <- .forest_pdp_cpp(object$reg_trees[keep], object$x[rows, , drop = FALSE], match(variable, x_names),
                           to_model(grid), p, rep(0, p), FALSE)
    for(j in seq_len(p)) {
      res$ice[, , j] <- res$ice[, , j] * object$rng_y[j] + object$center_y[j]
      res$pdp[, , j] <- res$pdp[, , j] * object$rng_y[j] + object$center_y[j]
    }
  } else if(object$model == "missBART2") {
    need_trees(object)
    Z <- cbind(if(object$include_x) object$x, if(object$include_y) y_model)[rows, , drop = FALSE]
    res <- .forest_pdp_cpp(object$class_trees[keep], Z, match(variable, preds), to_model(grid), p, object$z_offset, TRUE)
  } else {
    # probit regression: probabilities follow from the coefficients
    Z <- cbind(1, if(object$include_x) object$x, if(object$include_y) y_model)[rows, , drop = FALSE]
    col <- 1 + match(variable, preds)
    g <- to_model(grid)
    res <- list(ice = array(0, c(length(rows), length(grid), p)), pdp = array(0, c(length(keep), length(grid), p)))
    for(d in seq_along(keep)) {
      B <- object$B_post[[keep[d]]]
      base <- Z[, -col, drop = FALSE] %*% B[-col, , drop = FALSE]
      for(k in seq_along(g)) {
        pr <- stats::pnorm(sweep(base, 2, g[k] * B[col, ], "+"))
        res$ice[, k, ] <- res$ice[, k, ] + pr / length(keep)
        res$pdp[d, k, ] <- colMeans(pr)
      }
    }
  }
  pdp <- array(NA_real_, c(length(grid), p, 3), dimnames = list(NULL, y_names, c("mean", "lower", "upper")))
  pdp[, , "mean"] <- apply(res$pdp, c(2, 3), mean)
  pdp[, , "lower"] <- apply(res$pdp, c(2, 3), stats::quantile, probs = probs[1])
  pdp[, , "upper"] <- apply(res$pdp, c(2, 3), stats::quantile, probs = probs[2])
  dimnames(res$ice) <- list(NULL, NULL, y_names)
  structure(list(grid = grid, pdp = pdp, ice = res$ice, variable = variable, type = type, rows = rows),
            class = "missBART_pdp")
}

#' @rdname partial_dependence
#' @export
plot.missBART_pdp <- function(x, ...) {
  p <- dim(x$ice)[3]
  old <- graphics::par(mfrow = grDevices::n2mfrow(p), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old))
  for(j in seq_len(p)) {
    graphics::matplot(x$grid, t(x$ice[, , j]), type = "l", lty = 1, col = grDevices::grey(0.6, 0.4),
                      xlab = x$variable, main = dimnames(x$ice)[[3]][j],
                      ylab = if(x$type == "response") "Predicted response" else "Probability of being missing",
                      ylim = if(x$type == "missingness") c(0, 1) else range(x$ice[, , j], x$pdp[, j, ]), ...)
    graphics::polygon(c(x$grid, rev(x$grid)), c(x$pdp[, j, "lower"], rev(x$pdp[, j, "upper"])),
                      col = grDevices::adjustcolor("#FD7446", 0.3), border = NA)
    graphics::lines(x$grid, x$pdp[, j, "mean"], col = "#FD7446", lwd = 2)
  }
  invisible(x)
}

#' Variable interactions
#'
#' Counts how often two predictors are used in consecutive splitting rules on
#' the same branch of a tree (a parent and one of its children), as a measure
#' of how often they are used together.
#'
#' @param object an object of class `"missBART"` fitted with
#'   `keep_trees = TRUE`.
#' @param trees which set of trees: `"regression"` or `"missingness"` (only
#'   for [missBART2()] fits).
#'
#' @return A symmetric matrix. The off-diagonal entries are the average number
#'   of parent and child pairs per posterior sample that use the two
#'   predictors, and the diagonal entries the average number of splitting
#'   rules that use the predictor.
#' @export
interactions <- function(object, trees = c("regression", "missingness")) {
  stopifnot(inherits(object, "missBART"))
  trees <- match.arg(trees)
  need_trees(object)
  if(trees == "missingness" && object$model != "missBART2") stop("Only missBART2 fits have missingness trees")
  stored <- if(trees == "regression") object$reg_trees else object$class_trees
  names_ <- if(trees == "regression") colnames(object$x) else colnames(object$var_imp)
  q <- length(names_)
  counts <- numeric(q * q)
  for(tr in stored) {
    left <- tr[, 2]
    var <- tr[, 3]
    parent <- which(left > 0)
    child <- c(left[parent], left[parent] + 1)
    pv <- rep(var[parent], 2)
    internal <- var[child] > 0
    a <- pmin(pv[internal], var[child][internal])
    b <- pmax(pv[internal], var[child][internal])
    counts <- counts + tabulate((b - 1) * q + a, q * q)
  }
  out <- matrix(counts, q, q, dimnames = list(names_, names_)) / length(stored)
  same <- diag(out)
  out <- out + t(out)
  diag(out) <- colMeans(if(trees == "regression") object$var_imp_reg else object$var_imp)
  attr(out, "same_variable_pairs") <- same
  out
}

#' Agreement between chains
#'
#' For a fit with `n_chains` greater than one, compares summaries between the
#' chains. Selection models can have several modes, in which case chains
#' started from different values give different answers; see the real data
#' article on the package website.
#'
#' @param object an object of class `"missBART"` fitted with `n_chains > 1`.
#'
#' @return A list with `by_chain`, a matrix of the posterior mean of each
#'   quantity in each chain, and `rhat`, the potential scale reduction factor
#'   of Gelman and Rubin for each quantity (values well above 1 mean that the
#'   chains disagree). The quantities are the mean imputed value of each
#'   response, the residual variances and, for [missBART1()], the probit
#'   coefficients of the responses.
#' @export
chain_diagnostics <- function(object) {
  stopifnot(inherits(object, "missBART"))
  if(object$n_chains < 2) stop("This needs a fit with n_chains > 1")
  y_names <- colnames(object$y)
  draws <- sapply(object$sigma_post, diag)
  draws <- matrix(draws, ncol = length(object$sigma_post))
  rownames(draws) <- paste0("variance(", y_names, ")")
  if(length(object$missing_index) > 0) {
    col <- (object$missing_index - 1) %/% nrow(object$y) + 1
    imp <- sapply(object$y_impute, function(v) tapply(v, factor(col, levels = seq_along(y_names)), mean))
    imp <- matrix(imp, ncol = length(object$y_impute))
    rownames(imp) <- paste0("mean imputed ", y_names)
    draws <- rbind(imp[!is.na(imp[, 1]), , drop = FALSE], draws)
  }
  if(object$model == "missBART1" && object$include_y) {
    b <- sapply(object$B_post, function(B) c(B[y_names, , drop = FALSE]))
    b <- matrix(b, ncol = length(object$B_post))
    rownames(b) <- paste0("B[", rep(y_names, length(y_names)), ", ", rep(y_names, each = length(y_names)), "]")
    draws <- rbind(draws, b)
  }
  chain <- object$chain
  by_chain <- t(apply(draws, 1, function(v) tapply(v, chain, mean)))
  colnames(by_chain) <- paste("chain", sort(unique(chain)))
  n <- length(chain) / object$n_chains
  within <- rowMeans(t(apply(draws, 1, function(v) tapply(v, chain, stats::var))))
  between <- n * apply(by_chain, 1, stats::var)
  rhat <- sqrt(((n - 1) / n * within + between / n) / within)
  list(by_chain = by_chain, rhat = rhat)
}
