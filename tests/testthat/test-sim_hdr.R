
test_that("glmnet_cvmin works", {
  n <- 50
  p <- 100
  s0 <- 5
  x <- matrix(rnorm(n*p), nrow = n)
  beta <- c(rep(1, s0), rep(0, p - s0))
  y <- x %*% beta + rnorm(n)
  fitted_glmnet_cvmin <- fit_glmnet_cvmin(x, y, beta) #(x, y, beta, fitargs = NULL)
  expect_true(inherits(fitted_glmnet_cvmin$cv_fit, "cv.glmnet"))
  expect_true(inherits(fitted_glmnet_cvmin$glmnet_fit, "glmnet"))
  output_glmnet_cvmin <- post_glmnet_cvmin(fitted_glmnet_cvmin, x, y, beta)
  expect_equal(beta, output_glmnet_cvmin$true_beta[-1])
})

test_that("simulate_hdr works", {
  n <- 50
  p <- 100
  s0 <- 5
  niters <- 4
  result1 <- simulate_hdr(niters, n, p, s0, cores = 1, seed = 42)
  expect_length(result1, niters)

  set.seed(1)
  result2 <- simulate_hdr(niters, n, p, s0, cores = 1, seed = 42)
  expect_length(result2, niters)
  expect_true(identical(result1, result2))

  result3 <- simulate_hdr(niters, n, p, s0, cores = 1, seed = 1)
  expect_length(result3, niters)
  expect_false(identical(result1, result3))
})

test_that("simulate_hdr works with doRNG", {
  n <- 50
  p <- 100
  s0 <- 5
  niters <- 4
  result1 <- simulate_hdr(niters, n, p, s0, cores = 2, seed = 42)
  expect_length(result1, niters)

  set.seed(1)
  result2 <- simulate_hdr(niters, n, p, s0, cores = 2, seed = 42)
  expect_length(result2, niters)
  expect_true(identical(result1, result2))

  result3 <- simulate_hdr(niters, n, p, s0, cores = 2, seed = 1)
  expect_length(result3, niters)
  expect_false(identical(result1, result3))
})


