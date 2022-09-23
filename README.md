
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

This is a basic example which shows you how to solve a common problem:

``` r
library(simcity)
n <- 100
p <- 200
s0 <- 5
instance_hdr(n, p, s0)

niters <- 100
simulate_hdr(niters, n, p, s0, cores = 4)
```

