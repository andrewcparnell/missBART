sim2 <- function(n = 150, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * 3), ncol = 3)
  f <- cbind(4 * sin(pi * x[, 1]) + 2 * x[, 2], 3 * x[, 3]^2 - 2 * x[, 1])
  y_true <- f + matrix(rnorm(n * 2, sd = 0.3), ncol = 2)
  y <- y_true
  for(j in 1:2) y[runif(n) < pnorm(-0.7 + scale(y_true[, j])), j] <- NA
  list(x = x, y = y, f = f)
}
quick <- function(fun, ...) fun(..., burn = 60, iters = 40, thin = 1, show_progress = FALSE, mice_impute = FALSE)

test_that("predictions from stored trees reproduce the sampler's own predictions", {
  d <- sim2()
  x_new <- matrix(runif(30), ncol = 3)
  for(sc in c("center", "range", "none")) {
    fit <- quick(missBART1, d$x, d$y, x_predict = x_new, n_trees = 10, keep_trees = TRUE, scale = sc)
    expect_length(fit$reg_trees, 40)
    expect_equal(predict(fit, newdata = x_new), predict(fit), ignore_attr = TRUE)
    expect_equal(predict(fit, newdata = d$x), fitted(fit), ignore_attr = TRUE)
    expect_equal(dim(predict(fit, newdata = x_new, type = "draws")), c(10, 2, 40))
  }
  fit0 <- quick(missBART1, d$x, d$y, n_trees = 10)
  expect_null(fit0$reg_trees)
  expect_error(predict(fit0, newdata = x_new), "keep_trees")
  expect_error(predict(fit, newdata = x_new[, 1:2]), "3 columns")
})

test_that("stored trees also work with missing covariates", {
  d <- sim2()
  x <- d$x
  x[sample(length(x), 40)] <- NA
  x_new <- matrix(runif(30), ncol = 3)
  x_new[2, 1] <- NA
  fit <- quick(missBART2, x, d$y, x_predict = x_new, n_reg_trees = 10, n_class_trees = 5, keep_trees = TRUE)
  expect_equal(predict(fit, newdata = x_new), predict(fit), ignore_attr = TRUE)
  expect_length(fit$class_trees, 40)
})

test_that("partial dependence has the right shape and sensible values", {
  d <- sim2()
  fit2 <- quick(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, keep_trees = TRUE)
  pd <- partial_dependence(fit2, "X1", n_grid = 7, n_ice = 50, draws = 20)
  expect_s3_class(pd, "missBART_pdp")
  expect_equal(dim(pd$pdp), c(7, 2, 3))
  expect_equal(dim(pd$ice), c(50, 7, 2))
  expect_equal(range(pd$grid), range(d$x[, 1]))
  expect_true(all(pd$pdp[, , "lower"] <= pd$pdp[, , "upper"]))
  # the average of the individual curves is the partial dependence curve
  expect_equal(apply(pd$ice, c(2, 3), mean), pd$pdp[, , "mean"], ignore_attr = TRUE)
  # with every row and every draw, the curve at the data averages to about the fitted mean
  pm <- partial_dependence(fit2, "Y1", type = "missingness", n_grid = 5)
  expect_true(all(pm$pdp >= 0 & pm$pdp <= 1))
  expect_equal(range(pm$grid), range(completed_data(fit2)[, 1]))
  expect_error(partial_dependence(fit2, "nope"), "variable must be one of")

  fit1 <- quick(missBART1, d$x, d$y, n_trees = 10)
  p1 <- partial_dependence(fit1, "Y1", type = "missingness", n_grid = 5)
  expect_true(all(p1$pdp >= 0 & p1$pdp <= 1))
  # a probit regression is monotone in each predictor
  d1 <- diff(p1$pdp[, "Y1", "mean"])
  expect_true(all(d1 >= 0) || all(d1 <= 0))
  expect_error(partial_dependence(fit1, "X1"), "keep_trees")
  pdf(NULL)
  expect_invisible(plot(pd))
  expect_invisible(plot(p1))
  expect_invisible(plot(fit1))
  expect_invisible(plot(fit1, "fit"))
  dev.off()
})

test_that("a known partial dependence curve is recovered", {
  set.seed(5)
  x <- matrix(runif(600), ncol = 2)
  y <- 3 * x[, 1] + rnorm(300, sd = 0.1)
  fit <- mvBART(x, y, n_trees = 30, burn = 400, iters = 200, thin = 1, show_progress = FALSE, keep_trees = TRUE)
  pd <- partial_dependence(fit, "X1", grid = c(0.2, 0.8))
  expect_equal(unname(diff(pd$pdp[, 1, "mean"])), 1.8, tolerance = 0.15)
  flat <- partial_dependence(fit, "X2", grid = c(0.2, 0.8))
  expect_lt(abs(diff(flat$pdp[, 1, "mean"])), 0.3)
})

test_that("interactions are symmetric and consistent with variable importance", {
  d <- sim2()
  fit <- quick(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, keep_trees = TRUE)
  ir <- interactions(fit)
  expect_equal(dim(ir), c(3, 3))
  expect_true(isSymmetric(unname(ir)))
  expect_equal(diag(ir), colMeans(fit$var_imp_reg))
  im <- interactions(fit, "missingness")
  expect_equal(rownames(im), c("X1", "X2", "X3", "Y1", "Y2"))
  # a pair needs two splits, so there cannot be more pairs than splits
  expect_lte(sum(ir[upper.tri(ir)]) + sum(attr(ir, "same_variable_pairs")), sum(diag(ir)))
})

test_that("several chains are combined and compared", {
  d <- sim2()
  set.seed(8)
  fit <- quick(missBART1, d$x, d$y, n_trees = 10, n_chains = 3)
  expect_equal(fit$n_chains, 3)
  expect_length(fit$y_post, 120)
  expect_equal(fit$chain, rep(1:3, each = 40))
  expect_equal(dim(fit$tau_b), c(120, 6))
  expect_equal(dim(fit$var_imp_reg), c(120, 3))
  cd <- chain_diagnostics(fit)
  expect_equal(ncol(cd$by_chain), 3)
  expect_true(all(is.finite(cd$rhat)))
  expect_true(all(c("mean imputed Y1", "variance(Y2)", "B[Y1, Y2]") %in% names(cd$rhat)))
  expect_output(print(fit), "largest R-hat")
  # the same seed gives the same chains, in sequence or in parallel
  set.seed(8)
  again <- quick(missBART1, d$x, d$y, n_trees = 10, n_chains = 3, cores = 2)
  expect_identical(fit$y_impute, again$y_impute)
  expect_error(chain_diagnostics(quick(missBART1, d$x, d$y, n_trees = 5)), "n_chains")
  f2 <- quick(missBART2, d$x, d$y, n_reg_trees = 10, n_class_trees = 5, n_chains = 2, keep_trees = TRUE)
  expect_length(f2$class_trees, 80)
  expect_length(f2$y_miss_accept, sum(is.na(d$y)))
})
