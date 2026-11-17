#' Joint BART models for responses with non-ignorable missingness
#'
#' `missBART1()` and `missBART2()` fit the two selection models of Goh, Soh,
#' Parnell and Murphy (2024). Both use multivariate BART for the responses.
#' `missBART1()` models the missingness indicators with a multivariate probit
#' regression on `(1, x, y)`; `missBART2()` uses a multivariate probit BART
#' model on `(x, y)`. Missing responses are imputed within the sampler.
#' `missBARTprobit()` is the former name of `missBART1()`.
#'
#' @param x numeric matrix (or data frame) of covariates. Missing covariate
#'   values are allowed in `missBART2()` and `mvBART()`, where they are handled
#'   inside the trees; `missBART1()` needs complete covariates.
#' @param y numeric matrix (or vector) of responses, with `NA` for missing
#'   values.
#' @param x_predict optional covariates at which to make out-of-sample
#'   predictions.
#' @param n_trees,n_reg_trees number of trees in the BART model for the
#'   responses.
#' @param n_class_trees number of trees in the probit BART model for the
#'   missingness indicators.
#' @param burn number of burn-in iterations.
#' @param iters number of posterior samples kept after burn-in and thinning.
#' @param thin thinning interval. The sampler runs `burn + thin * iters`
#'   iterations.
#' @param predict logical. Make predictions at `x_predict`?
#' @param tree_prior_params tree prior and proposal settings, see [tree_list()].
#' @param hypers hyperparameters of the BART model for the responses, see
#'   [hypers_list()].
#' @param scale how the responses are standardised before fitting; see
#'   Details. `"center"` (the default) subtracts the mean of the observed
#'   values of each response and divides by their range, `"range"` (or `TRUE`)
#'   maps each response to `[-0.5, 0.5]` as in Chipman et al. (2010) and
#'   versions before 0.2.0, and `"none"` (or `FALSE`) leaves them alone.
#'   Unless `scale = "none"`, `missBART1()` also scales the covariates to
#'   `[-0.5, 0.5]`, so `B_post` refers to standardised variables. Everything
#'   else is returned on the original scale.
#' @param include_x,include_y logical. Use `x` and/or `y` as predictors of the
#'   missingness? `include_y = FALSE` assumes the responses are missing at
#'   random.
#' @param show_progress logical. Print the iteration count?
#' @param progress_every how often to update the iteration count.
#' @param pdp_range,make_pdp for `missBART2()` with a univariate response:
#'   if `make_pdp = TRUE`, the latent missingness predictor is evaluated over a
#'   grid of response values spanning `pdp_range` (on the scale of the data;
#'   the default is the observed range) at every kept iteration and returned
#'   in `pdp_out`, with the grid in `pdp_grid`; apply [stats::pnorm()] to get
#'   the probability of being missing. Ignored by `missBART1()`, where the
#'   same curves follow directly from `B_post`.
#' @param mice_impute logical. Initialise the missing responses with
#'   [mice::mice()]? Otherwise they start at the centre of the observed range.
#' @param MH_sd standard deviation of the random-walk proposal for missing
#'   responses in `missBART2()`. Defaults to `0.5 / ncol(y)`.
#' @param ymis_update how the Gaussian update of the missing responses treats
#'   rows that are only partly missing. `"conditional"` (the default)
#'   conditions on the observed responses in the same row. `"marginal"` draws
#'   from the corresponding margin of the full-row distribution, as written in
#'   Appendix C.1 of the paper and done by versions before 0.2.0. The two agree
#'   for univariate responses. Used by `missBART1()`, and by `missBART2()` when
#'   `include_y = FALSE`.
#' @param probit_prior gamma prior settings for the precisions of the probit
#'   coefficients in `missBART1()`, see [probit_prior_list()].
#' @param save_y_pred logical. Keep posterior predictive draws for the
#'   training data (`y_pred`)?
#' @param save_z logical. Keep the draws of the latent probit variables
#'   (`z_post`)?
#' @param ... further arguments passed to [tree_list()] and [hypers_list()].
#'   Arguments of earlier versions that are no longer used are ignored.
#'
#' @details
#' **Missingness indicator.** Write \eqn{M_{ij} = 1} if response \eqn{j} of
#' observation \eqn{i} is missing and 0 if it is observed. Both models use a
#' latent variable \eqn{M^*_{ij}} that is positive exactly when
#' \eqn{M_{ij} = 1}. In `missBART1()` its mean is \eqn{B^\top Z_i}, so a
#' positive coefficient in `B_post` means that larger values of the predictor
#' make the response more likely to be *missing*. (Versions before 0.2.0, and
#' the first version of the paper, used the opposite coding.)
#'
#' **Prior for the responses.** Following the `BART` package (Sparapani,
#' Spanbauer and McCulloch, 2021), the default `scale = "center"` fits
#' \eqn{Y_{ij} = \mu_{0j} + f_j(X_i) + \epsilon_{ij}}, with
#' \eqn{\mu_{0j}} the mean of the observed values of response \eqn{j} and
#' terminal node parameters centred at zero with standard deviation
#' \eqn{(\max_j - \min_j) / (2 k \sqrt{K})} for \eqn{K} trees, where
#' \eqn{k = \sqrt{\kappa / (4K)}} is 1.96 by default (see [hypers_list()]).
#' Internally this is done by standardising each response, which is equivalent.
#' The probit BART model of `missBART2()` likewise has an offset
#' \eqn{\mu_{0j} = \Phi^{-1}(\bar M_j)}, the probit of the proportion of
#' missing values of response \eqn{j}.
#'
#' **Tree updates.** Each tree is updated by a grow, prune or change proposal
#' accepted with a Metropolis-Hastings step. The acceptance ratio includes
#' the proposal ratio of the grow and prune moves.
#'
#' @return An object of class `"missBART"`: a list with, among others,
#' \describe{
#' \item{y_post}{list of `iters` matrices with the BART mean function for
#'   each training row.}
#' \item{y_pred}{list of posterior predictive draws for the training rows.}
#' \item{y_impute}{list of vectors of imputed values, ordered as
#'   `which(is.na(y))`.}
#' \item{new_y_hat, new_y_post}{lists with the mean function and posterior
#'   predictive draws at `x_predict` (or `NA`).}
#' \item{sigma_post}{list of residual covariance matrices. `omega_post` is the
#'   same object under its old name.}
#' \item{B_post, R_post, tau_b}{`missBART1()` only: probit coefficients,
#'   correlation matrices, and prior precisions of the coefficients.}
#' \item{var_imp}{`missBART2()` only: matrix with the number of splitting
#'   rules on each predictor in the missingness trees, one row per kept
#'   iteration.}
#' \item{var_imp_reg}{the same for the regression trees.}
#' \item{z_post}{list of latent probit variables (positive for missing
#'   responses).}
#' \item{y_miss_accept}{`missBART2()` only: Metropolis-Hastings acceptance
#'   rate for each missing value.}
#' \item{tree_accept}{acceptance rates of the tree proposals.}
#' }
#' Use [fitted()][fitted.missBART], [predict()][predict.missBART],
#' [imputed_values()], [completed_data()], [coef()][coef.missBART] and
#' [varimp()] to summarise a fit.
#'
#' @references Chipman, H. A., George, E. I. and McCulloch, R. E. (2010). BART:
#'   Bayesian additive regression trees. *The Annals of Applied Statistics*,
#'   4(1), 266-298.
#'
#'   Sparapani, R., Spanbauer, C. and McCulloch, R. (2021). Nonparametric
#'   machine learning and efficient computation with Bayesian additive
#'   regression trees: the BART R package. *Journal of Statistical Software*,
#'   97(1), 1-66. <doi:10.18637/jss.v097.i01>
#'
#'   Goh, Y. C., Soh, W. K., Parnell, A. C. and Murphy, K. (2024).
#'   Joint models for handling non-ignorable missing data using Bayesian
#'   additive regression trees: application to leaf photosynthetic traits
#'   data. <doi:10.48550/arXiv.2412.14946>
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(runif(n * 3), ncol = 3)
#' y_true <- cbind(sin(2 * pi * x[, 1]), x[, 2]^2) + matrix(rnorm(2 * n, sd = 0.2), ncol = 2)
#' # large values of the first response are more likely to be missing
#' y <- y_true
#' y[runif(n) < pnorm(2 * y_true[, 1] - 1), 1] <- NA
#'
#' fit <- missBART1(x, y, n_trees = 20, burn = 100, iters = 100, thin = 1,
#'                  mice_impute = FALSE, show_progress = FALSE)
#' fit
#' head(imputed_values(fit))
#' coef(fit)
#'
#' fit2 <- missBART2(x, y, n_reg_trees = 20, n_class_trees = 10, burn = 100,
#'                   iters = 100, thin = 1, mice_impute = FALSE,
#'                   show_progress = FALSE)
#' varimp(fit2)
#' @name missBART
NULL

#' @rdname missBART
#' @export
missBART1 <- function(x, y, x_predict = NULL, n_trees = 100, burn = 1000, iters = 1000, thin = 2,
                      predict = TRUE, tree_prior_params = tree_list(...), hypers = hypers_list(...),
                      scale = "center", include_x = TRUE, include_y = TRUE, show_progress = TRUE,
                      progress_every = 10, pdp_range = NULL, make_pdp = FALSE,
                      mice_impute = TRUE, ymis_update = c("conditional", "marginal"),
                      probit_prior = probit_prior_list(), save_y_pred = TRUE, save_z = TRUE, ...) {
  fit_missbart(model = 1L, x = x, y = y, x_predict = x_predict, n_reg_trees = n_trees,
               n_class_trees = 0L, burn = burn, iters = iters, thin = thin, predict = predict,
               tree_prior_params = tree_prior_params, hypers = hypers, scale = scale,
               include_x = include_x, include_y = include_y, show_progress = show_progress,
               progress_every = progress_every, mice_impute = mice_impute,
               ymis_update = match.arg(ymis_update), probit_prior = probit_prior,
               save_y_pred = save_y_pred, save_z = save_z, call = match.call())
}

#' @rdname missBART
#' @export
missBARTprobit <- missBART1

#' @rdname missBART
#' @export
missBART2 <- function(x, y, x_predict = NULL, n_reg_trees = 100, n_class_trees = 100, burn = 1000,
                      iters = 1000, thin = 2, predict = TRUE, MH_sd = NULL,
                      tree_prior_params = tree_list(...), hypers = hypers_list(...), scale = "center",
                      include_x = TRUE, include_y = TRUE, show_progress = TRUE, progress_every = 10,
                      pdp_range = NULL, make_pdp = FALSE, mice_impute = TRUE,
                      ymis_update = c("conditional", "marginal"), save_y_pred = TRUE,
                      save_z = TRUE, ...) {
  if(!include_x && !include_y) stop("At least one of include_x and include_y must be TRUE")
  fit_missbart(model = 2L, x = x, y = y, x_predict = x_predict, n_reg_trees = n_reg_trees,
               n_class_trees = n_class_trees, burn = burn, iters = iters, thin = thin,
               predict = predict, tree_prior_params = tree_prior_params, hypers = hypers,
               scale = scale, include_x = include_x, include_y = include_y,
               show_progress = show_progress, progress_every = progress_every,
               mice_impute = mice_impute, ymis_update = match.arg(ymis_update), MH_sd = MH_sd,
               pdp_range = pdp_range, make_pdp = make_pdp, save_y_pred = save_y_pred,
               save_z = save_z, call = match.call())
}

#' Multivariate BART for complete responses
#'
#' Fits the multivariate BART model used for the responses in [missBART1()]
#' and [missBART2()] to data without missing responses.
#'
#' @inheritParams missBART
#' @param y numeric matrix (or vector) of responses without missing values.
#' @param n_trees number of trees.
#'
#' @return An object of class `"missBART"`, see [missBART1()].
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- sim_data_friedman(n = 100, p = 2)
#' fit <- mvBART(data$x, data$y, n_trees = 20, burn = 100, iters = 100,
#'               thin = 1, show_progress = FALSE)
#' plot(data$mean[, 1], fitted(fit)[, 1])
mvBART <- function(x, y, x_predict = NULL, n_trees = 100, burn = 1000, iters = 1000, thin = 2,
                   predict = TRUE, tree_prior_params = tree_list(...), hypers = hypers_list(...),
                   scale = "center", show_progress = TRUE, progress_every = 10, save_y_pred = TRUE, ...) {
  if(anyNA(y)) stop("mvBART needs complete responses; use missBART1() or missBART2() when y has missing values")
  fit_missbart(model = 0L, x = x, y = y, x_predict = x_predict, n_reg_trees = n_trees,
               n_class_trees = 0L, burn = burn, iters = iters, thin = thin, predict = predict,
               tree_prior_params = tree_prior_params, hypers = hypers, scale = scale,
               show_progress = show_progress, progress_every = progress_every,
               mice_impute = FALSE, save_y_pred = save_y_pred, call = match.call())
}

as_numeric_matrix <- function(a, name) {
  if(is.data.frame(a)) a <- as.matrix(a)
  if(is.null(dim(a))) a <- matrix(a, ncol = 1)
  if(!is.numeric(a)) stop(name, " must be numeric")
  storage.mode(a) <- "double"
  a
}

range01 <- function(a, min, rng) sweep(sweep(a, 2, min, "-"), 2, rng, "/") - 0.5

# Data preparation, prior calibration and post-processing shared by the three models.
fit_missbart <- function(model, x, y, x_predict, n_reg_trees, n_class_trees, burn, iters, thin,
                         predict, tree_prior_params, hypers, scale, include_x = TRUE,
                         include_y = TRUE, show_progress, progress_every, mice_impute,
                         ymis_update = "conditional", probit_prior = probit_prior_list(),
                         MH_sd = NULL, pdp_range = NULL, make_pdp = FALSE,
                         save_y_pred = TRUE, save_z = FALSE, call) {
  y <- as_numeric_matrix(y, "y")
  x <- as_numeric_matrix(x, "x")
  n <- nrow(y)
  p <- ncol(y)
  if(nrow(x) != n) stop("x and y must have the same number of rows")
  if(is.null(colnames(y))) colnames(y) <- paste0("Y", seq_len(p))
  if(is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  if(any(!is.finite(y) & !is.na(y)) || any(is.infinite(x))) stop("x and y must not contain infinite values")
  for(a in c("n_reg_trees", "burn", "iters", "thin")) {
    v <- get(a)
    if(length(v) != 1 || is.na(v) || v < (a != "burn") || v != round(v)) stop(a, " must be a ", if(a == "burn") "non-negative" else "positive", " whole number")
  }
  if(model == 2L && (length(n_class_trees) != 1 || n_class_trees < 1)) stop("n_class_trees must be a positive whole number")

  predict <- isTRUE(predict) && !is.null(x_predict) && !all(is.na(x_predict))
  if(predict) {
    x_predict <- as_numeric_matrix(x_predict, "x_predict")
    if(ncol(x_predict) != ncol(x)) stop("x_predict must have the same number of columns as x")
    colnames(x_predict) <- colnames(x)
  } else {
    x_predict <- matrix(0, nrow = 0, ncol = ncol(x))
  }

  m <- matrix(1L, nrow = n, ncol = p)
  m[is.na(y)] <- 0L
  missing_index <- which(is.na(y))
  if(any(colSums(m) < 2)) stop("Every response needs at least two observed values")
  y_original <- y

  #-- scaling
  min_y <- apply(y, 2, min, na.rm = TRUE)
  max_y <- apply(y, 2, max, na.rm = TRUE)
  if(isTRUE(scale)) scale <- "range" else if(isFALSE(scale)) scale <- "none"
  scale <- match.arg(scale, c("center", "range", "none"))
  # standardised response = (y - center_y) / rng_y
  if(scale == "none") {
    center_y <- rep(0, p)
    rng_y <- rep(1, p)
  } else {
    rng_y <- max_y - min_y
    rng_y[rng_y == 0] <- 1
    center_y <- if(scale == "center") colMeans(y, na.rm = TRUE) else min_y + rng_y / 2
  }
  y <- sweep(sweep(y, 2, center_y, "-"), 2, rng_y, "/")
  names(center_y) <- colnames(y)
  min_x <- max_x <- NULL
  if(model == 1L) {
    if(anyNA(x) || anyNA(x_predict)) stop("missBART1 needs complete covariates: impute x first, or use missBART2()")
    if(scale != "none") {
      min_x <- apply(x, 2, min)
      max_x <- apply(x, 2, max)
      rng_x <- max_x - min_x
      rng_x[rng_x == 0] <- 1
      x <- range01(x, min_x, rng_x)
      if(predict) x_predict <- range01(x_predict, min_x, rng_x)
    }
  } else if(anyNA(x)) {
    # missingness indicators for incomplete covariates, so trees can also split on them
    na_cols <- which(colSums(is.na(x)) > 0)
    ind <- 1 - is.na(x[, na_cols, drop = FALSE])
    colnames(ind) <- paste0(colnames(x)[na_cols], "_observed")
    x_predict <- cbind(x_predict, 1 - is.na(x_predict[, na_cols, drop = FALSE]))
    x <- cbind(x, ind)
    colnames(x_predict) <- colnames(x)
  }
  q <- ncol(x)

  #-- prior calibration
  n_reg_trees <- as.integer(n_reg_trees)
  kappa <- if(is.null(hypers$kappa)) 4 * stats::qnorm(0.975)^2 * n_reg_trees else hypers$kappa
  nu <- hypers$df
  sigest <- vapply(seq_len(p), function(j) {
    s <- tryCatch(summary(stats::lm(y[, j] ~ x))$sigma, error = function(e) NaN)
    if(!is.finite(s) || s <= 0) s <- stats::sd(y[, j], na.rm = TRUE)
    s
  }, numeric(1))
  lambda <- sigest^2 * stats::qchisq(1 - hypers$q, nu) / nu
  alpha <- if(is.null(hypers$alpha)) nu else hypers$alpha
  if(is.null(hypers$V)) {
    Vinv <- diag(lambda * alpha, p)
  } else {
    if(!is.matrix(hypers$V) || any(dim(hypers$V) != p)) stop("hypers$V must be a p x p matrix")
    Vinv <- solve(hypers$V)
  }
  omega_init <- stats::rgamma(p, shape = nu / 2, rate = nu * lambda / 2)

  #-- starting values for the missing responses
  y_init <- y
  if(length(missing_index) > 0) {
    y_init[missing_index] <- 0
    if(mice_impute) {
      dat <- as.data.frame(cbind(y, x))
      names(dat) <- paste0("V", seq_len(ncol(dat)))
      imputed <- tryCatch(as.matrix(mice::complete(mice::mice(dat, printFlag = FALSE)))[, seq_len(p), drop = FALSE],
                          error = function(e) NULL)
      if(is.null(imputed)) {
        warning("mice imputation failed; missing responses start at the centre of the observed range")
      } else {
        ok <- missing_index[!is.na(imputed[missing_index])]
        y_init[ok] <- imputed[ok]
      }
    }
  }

  pp <- probit_prior
  if(is.null(pdp_range)) pdp_range <- c(min_y[1], max_y[1])
  pdp_grid <- seq(pdp_range[1], pdp_range[2], length.out = if(include_x) 20 else 50)
  ctrl <- list(n_reg_trees = n_reg_trees, n_class_trees = as.integer(n_class_trees),
               burn = as.integer(burn), iters = as.integer(iters), thin = as.integer(thin),
               kappa = kappa, kappa_class = 4 / 9 * max(n_class_trees, 1),
               mu0 = rep_len(hypers$mu0, p), wish_df = alpha, Vinv = Vinv, omega_init = omega_init,
               prior_alpha = tree_prior_params$prior_alpha, prior_beta = tree_prior_params$prior_beta,
               min_node = as.integer(max(tree_prior_params$min_node, p + 1)),
               max_attempt = as.integer(max(tree_prior_params$max_attempt, 1)),
               include_x = isTRUE(include_x), include_y = isTRUE(include_y),
               conditional = identical(ymis_update, "conditional"),
               MH_sd = if(is.null(MH_sd)) 0.5 / p else MH_sd,
               y_min = center_y - rng_y / 2, y_rng = rng_y,   # C++ maps v to (v + 0.5) * y_rng + y_min
               z_offset = stats::qnorm((colSums(m == 0L) + 0.5) / (n + 1)),
               a0 = pp$a0, b0 = pp$b0,
               ax = if(is.null(pp$ax)) q + 1 else pp$ax, bx = if(is.null(pp$bx)) p + q + 1 else pp$bx,
               ay = if(is.null(pp$ay)) 1 else pp$ay, by = if(is.null(pp$by)) 1 else pp$by,
               save_y_pred = isTRUE(save_y_pred), save_z = isTRUE(save_z),
               make_pdp = isTRUE(make_pdp),
               pdp_grid = (pdp_grid - center_y[1]) / rng_y[1],
               show_progress = isTRUE(show_progress), progress_every = as.integer(progress_every))

  out <- .missbart_cpp(x, y_init, m, x_predict, as.integer(model), ctrl)

  #-- tidy up
  y_names <- colnames(y)
  out$y_post <- lapply(out$y_post, function(a) { colnames(a) <- y_names; a })
  out$y_impute <- lapply(seq_len(nrow(out$y_impute)), function(i) out$y_impute[i, ])
  out$sigma_post <- lapply(out$sigma_post, function(a) { dimnames(a) <- list(y_names, y_names); a })
  out$omega_post <- out$sigma_post
  colnames(out$var_imp) <- colnames(x)
  out$var_imp_reg <- out$var_imp
  out$var_imp <- NULL
  if(!predict) out$new_y_post <- out$new_y_hat <- NA
  if(!ctrl$save_y_pred) out$y_pred <- NULL
  if(model == 1L) {
    z_names <- c("Intercept", if(include_x) colnames(x), if(include_y) y_names)
    out$B_post <- lapply(out$B_post, function(a) { dimnames(a) <- list(z_names, y_names); a })
    out$R_post <- lapply(out$R_post, function(a) { dimnames(a) <- list(y_names, y_names); a })
    colnames(out$tau_b) <- z_names
    out$var_imp_miss <- out$y_miss_accept <- NULL
  } else {
    out$B_post <- out$R_post <- out$tau_b <- NULL
    if(model == 2L) {
      colnames(out$var_imp_miss) <- c(if(include_x) colnames(x), if(include_y) y_names)
      out$var_imp <- out$var_imp_miss
      if(!include_y || length(missing_index) == 0) out$y_miss_accept <- NULL
    } else {
      out$y_miss_accept <- out$z_post <- NULL
      out$tree_accept <- out$tree_accept["regression"]
    }
    out$var_imp_miss <- NULL
  }
  if(model == 1L) out$tree_accept <- out$tree_accept["regression"]
  if(!ctrl$save_z) out$z_post <- NULL
  if(length(out$pdp_out) == 0) out$pdp_out <- NA else out$pdp_grid <- pdp_grid

  structure(c(out, list(model = c("mvBART", "missBART1", "missBART2")[model + 1L], call = call,
                        burn = burn, iters = iters, thin = thin, scale = scale,
                        min_y = min_y, max_y = max_y, center_y = center_y, min_x = min_x, max_x = max_x,
                        y = y_original, missing_index = missing_index,
                        include_x = include_x, include_y = include_y,
                        kappa = kappa, MH_sd = ctrl$MH_sd)),
            class = "missBART")
}
