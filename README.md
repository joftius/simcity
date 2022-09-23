
# simcity

<!-- badges: start -->
<!-- badges: end -->

The goal of simcity is to streamline simulating data to evaluate and compare methods for some common statistical tasks. 

## Installation

You can install the development version of simcity from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("joftius/simcity")
```

## Example

The first motivating example for this package was simulating high-dimensional linear regression models, fitting them with the lasso, and analyzing various performance metrics.

``` r
library(simcity)
n <- 100
p <- 200
s0 <- 5
one_lasso_fit <- instance_hdr(n, p, s0)
one_lasso_fit |> head()
which(one_lasso_fit$true_beta != 0)
which(one_lasso_fit$estimate != 0)
```

The above example generates one instance of simulated data, fits a regression model using `glmnet::cv.glmnet` and the `lambda.1se` option by default. The `simcity` package streamlines doing processes like this many times, and typically finishes in about half the time (or less) by using parallel processing.

``` r
niters <- 200
many_lasso_fits <- simulate_hdr(niters, n, p, s0, cores = 4)
sim_summary <- simmary_coefs(many_lasso_fits)
head(sim_summary)
ggplot(sim_summary, aes(screened, beta_min)) + geom_boxplot()
ggplot(sim_summary, aes(beta_min, mse)) + geom_point()
```

The above example generates a number `niters` of instances of simulated data, fits each of them using `cv.glmnet`, and computes some interesting summaries about the overall results.
