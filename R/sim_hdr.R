
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Simulation for high-dimensional linear regression model variable selection
#'
#' @description
#' `this` repeats `that` multiple times.
#'
#' @details
#' After generating simulated data this calls `something` for inference.
#'
#' @param n Sample size, number of rows in design matrix.
#' @param p Number of random predictor variables or columns in design matrix.
#' @param s0 Sparsity or number of nonzero coefficients in the true linear model.
#' @param xtype,btype,permuted,x.par Linear model data generation parameters, see \link[hdi]{rXb} for details.
#' @param yfun,yargs Function and (optional) arguments for generating outcome variable.
#' @param fitfun,fitargs Function and (optional) arguments for fitting models to simulated data.
#' @param postfun,postargs Function and (optional) arguments for post-processing fitted models.
#' @return Outputs from `postfun` after applying it to the model fitted by `fitfun` on one instance of simulated data.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
instance_hdr <- function(
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
                   "toeplitz"  = 1/3,
                   "equi.corr" = 1/20,
                   "exp.decay" = c(0.4, 5)
    ),
    yfun,
    yargs = NULL,
    fitfun,
    fitargs = NULL,
    postfun,
    postargs = NULL)
{
  xtype <- tolower(match.arg(xtype))

  sim_data <- rXb(n = n, p = p, s0 = s0, xtype = xtype,
                  btype = btype, permuted = permuted, x.par = x.par)
  x <- sim_data$x
  beta <- sim_data$beta
  if (is.null(yargs)) {
    y <- yfun(x, beta)
  } else {
    y <- yfun(x, beta, yargs)
  }
  # min_beta <- min(abs(beta[1:s0]))
  # max_corr <- max(sapply(1:s0,
  #                        function(j) {
  #                          max(abs(cor(x = x[, j], y = x[, (s0+1):p])))
  #                        }))

  if (is.null(fitargs)) {
    fit_obj <- fitfun(x, y, beta)
  } else {
    fit_obj <- fitfun(x, y, beta, fitargs)
  }
  if (is.null(postargs)) {
    return(postfun(fit_obj, x, y, beta))
  } else {
    return(postfun(fit_obj, x, y, beta, postargs))
  }
}

#' Simulation for high-dimensional linear regression model variable selection
#'
#' @description
#' `this` repeats `that` multiple times.
#'
#' @details
#' After generating simulated data this calls `something` for inference.
#'
#' @param niters Number of simulation iterations.
#' @inheritParams instance_hdr
#' @param cores Number of cores for parallel computation, passed to \link[parallel]{makeCluster}. Defaults to half of \link[parallel]{detectCores} when not specified.
#' @return List of outputs from `postfun` for each of `niters` instances.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
simulate_hdr <- function(
    niters,
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
                   "toeplitz"  = 1/3,
                   "equi.corr" = 1/20,
                   "exp.decay" = c(0.4, 5)
    ),
    yfun = y_standard_linear,
    yargs = NULL,
    fitfun = fit_glmnet_cvmin,
    fitargs = NULL,
    postfun = post_glmnet_cvmin,
    postargs = NULL,
    cores = NULL
)
{

  xtype <- tolower(match.arg(xtype))
  time_start <- Sys.time()
  local_files <- list.files(path = paste0(getwd(), "/R"), pattern = ".R")
  if (is.null(cores)) cores <- floor(detectCores()/2)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  #clusterCall(cl, function() library(selectiveInference))
  #clusterCall(cl, function() library(RPtests))
  #clusterCall(cl, function() library(unbiasedgoodness))
  clusterCall(cl, function() source(paste0(getwd(), "/R/sim_hdr.R")))

  # test_run <- instance_hdr(
  #   n,
  #   p,
  #   s0,
  #   sigma,
  #   xtype = c("toeplitz", "exp.decay", "equi.corr"),
  #   btype = "U[-2,2]",
  #   permuted = TRUE,
  #   x.par = switch(xtype,
  #                  "toeplitz"  = 1/3,
  #                  "equi.corr" = 1/20,
  #                  "exp.decay" = c(0.4, 5)
  #   ),
  #   fitfun,
  #   fitargs
  # )

  output <- foreach(icount(niters)) %dopar% {
    instance_hdr(
      n,
      p,
      s0,
      xtype = c("toeplitz", "exp.decay", "equi.corr"),
      btype = "U[-2,2]",
      permuted = TRUE,
      x.par = switch(xtype,
                     "toeplitz"  = 1/3,
                     "equi.corr" = 1/20,
                     "exp.decay" = c(0.4, 5)
      ),
      yfun,
      yargs,
      fitfun,
      fitargs,
      postfun,
      postargs
    )
  }
  stopImplicitCluster()
  time_end <- Sys.time()

  print(time_end - time_start)
  return(output)
}


y_standard_linear <- function(x, beta, yargs = NULL) {
  sigma <- 1
  if (!is.null(yargs)) sigma <- yargs$sigma
  x %*% beta + sigma * rnorm(nrow(x))
}

fit_cvglmnet <- function(x, y, beta, fitargs = NULL) {
  if (is.null(fitargs)) {
    return(cv.glmnet(x, y))
  } else {
    return(cv.glmnet(x, y, fitargs))
  }
}

post_cvglmnet <- function(fit_obj, x, y, beta, postargs = NULL) {
  if (!inherits(fit_obj, "cv.glmnet")) {
    warning("fit_obj should be of class cv.glmnet")
  }
}

fit_glmnet_cvmin <- function(x, y, beta, fitargs = NULL) {
  if (is.null(fitargs)) {
    cv_fit <- cv.glmnet(x, y)
    glmnet_fit <- glmnet(x, y)
  } else {
    cv_fit <- cv.glmnet(x, y, fitargs)
    glmnet_fit <- glmnet(x, y, fitargs)
  }
  list(cv_fit = cv_fit, glmnet_fit = glmnet_fit)
}

post_glmnet_cvmin <- function(fit_obj, x, y, beta, postargs = NULL) {
  lambda_type <- "lambda.1se"
  if (!is.null(postargs)) lambda_type <- postargs$lambda
  lambda_min <- fit_obj$cv_fit[[lambda_type]]
  list(true_beta = beta,
       beta_hat = coef.glmnet(fit_obj$glmnet_fit, s = lambda_min))
}

# test run
#source("R/simulation_tools.R")
# one_thing  <- instance_hdr(100, 200, 5, yfun = y_standard_linear, yargs = NULL, fitfun = fit_glmnet_cvmin,fitargs = NULL, postfun = post_glmnet_cvmin,postargs = NULL)
# things  <- simulate_hdr(10, 100, 200, 5)
