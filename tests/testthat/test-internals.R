test_that("prior helpers validate their input", {
  expect_length(tree_list(), 4)
  expect_error(tree_list(prior_alpha = -0.1))
  expect_error(tree_list(prior_beta = -1))
  expect_named(hypers_list(), c("mu0", "kappa", "alpha", "V", "df", "q", "tau_b"))
  expect_error(probit_prior_list(a0 = -1))
  expect_equal(probit_prior_list(ay = 10)$ay, 10)
})

test_that("scale_bart and unscale are inverses", {
  y <- matrix(rnorm(20), ncol = 2)
  s <- scale_bart(y)
  expect_equal(range(s[, 1]), c(-0.5, 0.5))
  expect_equal(unscale(s, apply(y, 2, min), apply(y, 2, max)), y)
  expect_equal(unscale(scale_bart(y[, 1]), min(y[, 1]), max(y[, 1])), y[, 1])
})

test_that("truncated normal draws respect the truncation and have the right mean", {
  set.seed(1)
  for(mu in c(-6, -1, 0, 0.3, 2)) {
    pos <- missBART:::.rtnorm_cpp(rep(mu, 2e4), 1, TRUE)
    neg <- missBART:::.rtnorm_cpp(rep(mu, 2e4), 1, FALSE)
    expect_true(all(pos > 0))
    expect_true(all(neg <= 0))
    expect_equal(mean(pos), mu + dnorm(mu) / pnorm(mu), tolerance = 0.03)
    expect_equal(mean(neg), mu - dnorm(mu) / pnorm(-mu), tolerance = 0.03)
  }
  z <- missBART:::.rtnorm_cpp(rep(1, 2e4), 2, TRUE)
  expect_equal(mean(z), 1 + 2 * dnorm(0.5) / pnorm(0.5), tolerance = 0.03)
})

test_that("Wishart and inverse-Wishart draws have the right mean", {
  set.seed(2)
  S <- matrix(c(2, 0.5, 0.5, 1), 2)
  W <- Reduce("+", lapply(1:4000, function(i) missBART:::.rwish_cpp(7, S, FALSE))) / 4000
  expect_equal(W, 7 * S, tolerance = 0.05)
  IW <- Reduce("+", lapply(1:4000, function(i) missBART:::.rwish_cpp(7, S, TRUE))) / 4000
  expect_equal(IW, S / (7 - 2 - 1), tolerance = 0.08)
})

test_that("terminal nodes agree with the R tree code, including missing covariates", {
  set.seed(3)
  x <- matrix(runif(300), ncol = 3)
  x[sample(300, 30)] <- NA
  # root splits on variable 2 at 0.5 (missing go right); its left child splits
  # on variable 1 at 0.3 (missing go left)
  df <- data.frame(parent = c(0, 1, 1, 2, 2), lower = 0, upper = 1,
                   split_variable = c(0, 2, 2, 1, 1), split_value = c(1, 0.5, 0.5, 0.3, 0.3),
                   depth = c(0, 1, 1, 2, 2), direction = c(0, 0, 1, 0, 1),
                   NA_direction = c(NA, 0, 1, 1, 0))
  r_nodes <- missBART:::get_change_points(df, x)
  cpp_nodes <- missBART:::.descend_cpp(x, left = c(1L, 3L, -1L, -1L, -1L), var = c(1L, 0L, 0L, 0L, 0L),
                                       val = c(0.5, 0.3, 0, 0, 0), na_left = c(FALSE, TRUE, FALSE, FALSE, FALSE))
  # same partition of the observations (the labels differ)
  tab <- table(r_nodes, cpp_nodes)
  expect_equal(dim(tab), c(3L, 3L))
  expect_equal(sum(tab > 0), 3)
})

test_that("node marginal likelihood matches the formula in the paper", {
  set.seed(4)
  ref <- function(res, kappa, omega, mu0) {
    n <- nrow(res)
    p <- ncol(res)
    Sinv <- n * omega + diag(kappa, p)
    S <- solve(Sinv)
    mu <- S %*% (omega %*% colSums(res) + kappa * mu0)
    as.numeric(0.5 * (-n * p * log(2 * pi) + p * log(kappa) + n * determinant(omega)$modulus +
                        determinant(S)$modulus - kappa * sum(mu0^2) + t(mu) %*% Sinv %*% mu -
                        sum((res %*% omega) * res)))
  }
  for(p in 1:3) {
    omega <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    mu0 <- rnorm(p, sd = 0.1)
    res <- matrix(rnorm(40 * p), ncol = p)
    left <- res[1:15, , drop = FALSE]
    right <- res[16:40, , drop = FALSE]
    cpp <- function(a) missBART:::.leaf_loglik_cpp(nrow(a), colSums(a), omega, 3, mu0)
    # terms common to every tree cancel in the difference used by the sampler
    expect_equal(cpp(left) + cpp(right) - cpp(res),
                 ref(left, 3, omega, mu0) + ref(right, 3, omega, mu0) - ref(res, 3, omega, mu0))
  }
})

test_that("Gaussian update of missing responses has the right moments", {
  set.seed(5)
  p <- 3
  Q <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  b <- c(1, -2, 0.5)
  nrep <- 2e4
  m <- matrix(c(0L, 1L, 0L), nrow = nrep, ncol = p, byrow = TRUE)
  y0 <- matrix(c(0, 0.7, 0), nrow = nrep, ncol = p, byrow = TRUE)
  B <- matrix(b, nrow = nrep, ncol = p, byrow = TRUE)
  mis <- c(1, 3)

  cond <- missBART:::.draw_ymis_cpp(y0, m, Q, B, TRUE)
  expect_equal(cond[, 2], rep(0.7, nrep))
  cov_c <- solve(Q[mis, mis])
  mean_c <- cov_c %*% (b[mis] - Q[mis, 2] * 0.7)
  expect_equal(colMeans(cond[, mis]), as.numeric(mean_c), tolerance = 0.03)
  expect_equal(cov(cond[, mis]), cov_c, tolerance = 0.05)

  marg <- missBART:::.draw_ymis_cpp(y0, m, Q, B, FALSE)
  S <- solve(Q)
  expect_equal(colMeans(marg[, mis]), as.numeric(S %*% b)[mis], tolerance = 0.03)
  expect_equal(cov(marg[, mis]), S[mis, mis], tolerance = 0.05)
})
