
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


