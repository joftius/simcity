
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# tidy() summarizes information about model components
# glance() reports information about the entire model
# augment() adds information about observations to a dataset

#' instance_hdr
#' Simulate a single instance for high-dimensional linear regression models.
#'
#' @description
#' Simulates data from a high-dimensional linear model and then applies user-specified functions to fit a model and process the results.
#'
#' @details
#' 1. Generates a random design matrix and coefficient vector using [hdi::rXb()].
#' 2. Generates an outcome variable using `yfun`.
#' 3. Fits a model to the data using the method specified with `fitfun`.
#' 4. Processes the resulting model with `postfun`.
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
#' \dontrun{
#' instance_hdr(n = 100, p = 200, s0 = 5)
#' }
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
    yfun = y_standard_linear,
    yargs = NULL,
    fitfun = fit_glmnet_cvmin,
    fitargs = NULL,
    postfun = post_glmnet_cvmin,
    postargs = NULL)
{
  xtype <- tolower(match.arg(xtype))

  sim_data <- hdi::rXb(n = n, p = p, s0 = s0, xtype = xtype,
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

#' simulate_hdr
#' Simulation for high-dimensional linear regression models.
#'
#' @inherit instance_hdr description
#' @inherit instance_hdr details
#'
#' @param niters Number of simulation iterations.
#' @inheritParams instance_hdr
#' @param cores Number of cores for parallel computation, passed to \link[parallel]{makeCluster}. Defaults to half of \link[parallel]{detectCores} when not specified.
#' @param verbose If `TRUE` will print elapsed time.
#' @return List of outputs from `postfun` for each of `niters` instances.
#' @export
#' @examples
#' \dontrun{
#' simulate_hdr(niters = 500, n = 100, p = 200, s0 = 5)
#' }
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
    cores = NULL,
    seed = NULL,
    verbose = FALSE
)
{

  xtype <- tolower(match.arg(xtype))
  time_start <- Sys.time()
  local_files <- list.files(path = paste0(getwd(), "/R"), pattern = ".R")
  if (is.null(cores)) cores <- floor(detectCores()/2)
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  registerDoRNG(seed)
  #clusterCall(cl, function() library(glmnet))
  pb <- txtProgressBar(0, niters, style = 3)
  output <- foreach(i = icount(niters)) %dorng% {
    setTxtProgressBar(pb, i)
    instance_hdr(
      n,
      p,
      s0,
      xtype = xtype,
      btype = btype,
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
  stopCluster(cl)
  time_end <- Sys.time()

  if (verbose) print(time_end - time_start)
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
  beta_hat <- coef.glmnet(fit_obj$glmnet_fit,
                          s = lambda_min,
                          exact = TRUE,
                          x = x,
                          y = y)
  terms <- rownames(beta_hat)
  output <- data.frame(
    term = terms,
    estimate = beta_hat[1:nrow(beta_hat)],
    true_beta = if (terms[1] == "(Intercept)") c(0, beta) else beta)
 attributes(output)$lambda <- lambda_min
 attributes(output)$lambda_type <- lambda_type
 output
}

#' simmary_coefs
#' Compute simulation summary using true and estimated linear model coefficients.
#'
#' @description
#' Compute simulation summary using true and estimated linear model coefficients.
#'
#' @details
#'
#' @param list_of_instances The output of a simulation function.
#' @returns
#' A `data.frame` with one row for each simulation instance and several performance metrics:
#'
#' * `screened` Boolean indicating whether the selected model's support is a superset of the true support, i.e. if the selected variables include all the ones with nonzero coefficients in the true simulated model.
#' * `mse` The mean-square error of estimating the true coefficients.
#' * `beta_min` The smallest nonzero coefficient in the true model (in absolute value).
#' @export
#' @examples
#' \dontrun{
#' instances <- simulate_hdr(niters = 500, n = 100, p = 200, s0 = 5)
#' simmary_coefs(instances)
#' }
simmary_coefs <- function(list_of_instances) {
  list_of_instances |> purrr::map_dfr(
    function(lfit) {
      true_support <- which(lfit$true_beta != 0)
      estimated_support <- which(lfit$estimate != 0)
      recovered <- true_support %in% estimated_support
      data.frame(
        true_sparsity = length(true_support),
        estimated_sparsity = length(estimated_support),
        false_positives = length(setdiff(estimated_support, true_support)),
        true_positives = sum(recovered),
        screened = all(recovered),
        mse = mean((lfit$estimate - lfit$true_beta)^2),
        beta_min = min(abs(lfit$true_beta[which(lfit$true_beta !=0)]))
        )
      })
}

# test run
#source("R/simulation_tools.R")
# one_thing  <- instance_hdr(100, 200, 5, yfun = y_standard_linear, yargs = NULL, fitfun = fit_glmnet_cvmin,fitargs = NULL, postfun = post_glmnet_cvmin,postargs = NULL)
# things  <- simulate_hdr(10, 100, 200, 5)
