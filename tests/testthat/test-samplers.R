sim <- function(n = 150, p = 2, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * 3), ncol = 3)
  f <- cbind(4 * sin(pi * x[, 1]) + 2 * x[, 2], 3 * x[, 3]^2 - 2 * x[, 1])[, seq_len(p), drop = FALSE]
  y_true <- f + matrix(rnorm(n * p, sd = 0.3), ncol = p)
  y <- y_true
  for(j in seq_len(p)) y[runif(n) < pnorm(-0.7 + scale(y_true[, j])), j] <- NA
  list(x = x, y = y, y_true = y_true, f = f)
}
run <- function(fun, ..., burn = 150, iters = 100) {
  fun(..., burn = burn, iters = iters, thin = 1, show_progress = FALSE, mice_impute = FALSE)
}

test_that("missBART1 returns a well-formed fit", {
  d <- sim()
  x_new <- matrix(runif(30), ncol = 3)
  fit <- run(missBART1, d$x, d$y, x_predict = x_new, n_trees = 20)
  n_mis <- sum(is.na(d$y))
  expect_s3_class(fit, "missBART")
  expect_length(fit$y_post, 100)
  expect_equal(dim(fit$y_post[[1]]), dim(d$y))
  expect_equal(dim(fit$y_pred[[1]]), dim(d$y))
  expect_length(fit$y_impute[[1]], n_mis)
  expect_equal(dim(fit$B_post[[1]]), c(1 + 3 + 2, 2))
  expect_equal(unname(diag(fit$R_post[[1]])), c(1, 1))
  expect_equal(dim(fit$tau_b), c(100, 6))
  expect_equal(dim(fit$new_y_post[[1]]), c(10, 2))
  expect_equal(dim(predict(fit)), c(10, 2))
  expect_equal(dim(fitted(fit)), dim(d$y))
  expect_equal(dim(coef(fit)), c(6, 2, 3))
  iv <- imputed_values(fit)
  expect_equal(nrow(iv), n_mis)
  expect_equal(which(is.na(d$y)), fit$missing_index)
  cd <- completed_data(fit)
  expect_false(anyNA(cd))
  expect_equal(cd[!is.na(d$y)], d$y[!is.na(d$y)])
  expect_true(all(is.finite(unlist(fit$y_post))))
  expect_output(print(fit), "missBART1")
  expect_output(summary(fit), "probit coefficients")
  # the responses drive their own missingness, so the imputations should beat the observed mean
  naive <- colMeans(d$y, na.rm = TRUE)[col(d$y)[fit$missing_index]]
  expect_lt(mean((iv$mean - d$y_true[fit$missing_index])^2), mean((naive - d$y_true[fit$missing_index])^2))
})

test_that("missBART2 returns a well-formed fit", {
  d <- sim()
  fit <- run(missBART2, d$x, d$y, n_reg_trees = 20, n_class_trees = 10)
  expect_equal(fit$model, "missBART2")
  expect_equal(dim(fit$var_imp), c(100, 5))
  expect_equal(sum(varimp(fit)), 1)
  expect_length(fit$y_miss_accept, sum(is.na(d$y)))
  expect_true(all(fit$y_miss_accept >= 0 & fit$y_miss_accept <= 1))
  expect_identical(fit$new_y_post, NA)
  expect_error(predict(fit), "x_predict")
  expect_error(coef(fit))
  expect_equal(dim(fit$z_post[[1]]), dim(d$y))
  # the latent variables agree in sign with the missingness indicators
  expect_true(all((fit$z_post[[100]] > 0) == !is.na(d$y)))
  expect_lt(sqrt(mean((fitted(fit) - d$f)^2)), 1)
})

test_that("fits are reproducible and options work", {
  d <- sim(n = 100)
  set.seed(10)
  a <- run(missBART1, d$x, d$y, n_trees = 10, burn = 20, iters = 20)
  set.seed(10)
  b <- run(missBART1, d$x, d$y, n_trees = 10, burn = 20, iters = 20)
  expect_identical(a$y_post, b$y_post)
  expect_identical(a$y_impute, b$y_impute)
  set.seed(10)
  expect_identical(run(missBARTprobit, d$x, d$y, n_trees = 10, burn = 20, iters = 20)$B_post, a$B_post)

  f1 <- run(missBART1, d$x, d$y, n_trees = 10, burn = 20, iters = 20, ymis_update = "marginal",
            include_x = FALSE, save_y_pred = FALSE, save_z = FALSE, scale = FALSE)
  expect_equal(rownames(f1$B_post[[1]]), c("Intercept", "Y1", "Y2"))
  expect_null(f1$y_pred)
  expect_null(f1$z_post)
  f2 <- run(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, burn = 20, iters = 20, include_y = FALSE)
  expect_equal(colnames(f2$var_imp), paste0("X", 1:3))
  expect_null(f2$y_miss_accept)
  expect_error(missBART2(d$x, d$y, include_x = FALSE, include_y = FALSE))
  thinned <- missBART1(d$x, d$y, n_trees = 5, burn = 7, iters = 9, thin = 3, show_progress = FALSE, mice_impute = FALSE)
  expect_length(thinned$y_post, 9)
})

test_that("univariate responses and partial dependence work", {
  d <- sim(p = 1)
  f1 <- run(missBART1, d$x, d$y[, 1], n_trees = 10, burn = 50, iters = 50)
  expect_equal(dim(f1$B_post[[1]]), c(5, 1))
  f2 <- run(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, burn = 50, iters = 50, make_pdp = TRUE)
  expect_equal(dim(f2$pdp_out[[1]]), c(nrow(d$x), 20))
  f3 <- run(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, burn = 50, iters = 50,
            make_pdp = TRUE, include_x = FALSE)
  expect_equal(dim(f3$pdp_out[[1]]), c(1, 50))
})

test_that("missing covariates are handled by the trees", {
  d <- sim()
  x <- d$x
  x[sample(length(x), 40)] <- NA
  x_new <- matrix(runif(30), ncol = 3)
  x_new[1, 2] <- NA
  fit <- run(missBART2, x, d$y, x_predict = x_new, n_reg_trees = 10, n_class_trees = 5, burn = 30, iters = 30)
  expect_true(all(is.finite(unlist(fit$y_post))))
  expect_true(all(is.finite(predict(fit))))
  expect_true(all(c("X1_observed", "X2_observed", "X3_observed") %in% colnames(fit$var_imp_reg)))
  expect_error(missBART1(x, d$y), "complete covariates")
  ok <- complete.cases(d$y)
  fit0 <- run(mvBART, x[ok, ], d$y[ok, ], n_trees = 10, burn = 30, iters = 30)
  expect_equal(fit0$model, "mvBART")
})

test_that("mvBART recovers a smooth signal and input is checked", {
  set.seed(2)
  data <- sim_data_friedman(n = 200, p = 2)
  fit <- run(mvBART, data$x, data$y, n_trees = 30, burn = 300, iters = 200)
  expect_gt(cor(fitted(fit)[, 1], data$mean[, 1]), 0.9)
  expect_null(fit$B_post)
  y <- data$y
  y[1, 1] <- NA
  expect_error(mvBART(data$x, y), "complete responses")
  expect_error(missBART1(data$x[-1, ], data$y), "same number of rows")
  expect_error(missBART1(data$x, data$y, x_predict = matrix(0, 2, 2)), "same number of columns")
  expect_error(missBART1(data$x, data$y, burn = -1))
})

test_that("mice initialisation runs", {
  skip_on_cran()
  d <- sim(n = 80)
  fit <- missBART1(d$x, d$y, n_trees = 5, burn = 5, iters = 5, thin = 1, show_progress = FALSE)
  expect_length(fit$y_post, 5)
})

test_that("simulation helpers run", {
  set.seed(3)
  data <- sim_data_trees(n = 200, p = 2, q = 3, trees = 3)
  expect_equal(dim(data$y), c(200, 2))
  miss <- sim_missing_probit(data$x, data$y, include_x = TRUE, include_y = TRUE,
                             min_missing_prop = 0.3, max_missing_prop = 0.95)
  expect_equal(dim(miss$m), c(200, 2))
  miss2 <- sim_missing_trees(data$x, data$y, min_missing_prop = 0.2, max_missing_prop = 0.98)
  expect_equal(dim(miss2$m), c(200, 2))
})
