# emhawkes

`emhawkes` specifies, simulates, and fits exponential Hawkes models,
including marked, state-dependent, and multi-kernel models.

## Installation

Install the CRAN release:

``` r

install.packages("emhawkes")
```

For the development version and the features described in the guides:

``` r

# install.packages("remotes")
remotes::install_github("ksublee/emhawkes")
```

Source installation requires a C++ compiler and Rcpp; on Windows,
install Rtools.

## Quick start

This example requires version 1.0.0 or later.

``` r

library(emhawkes)

h <- new("hspec", mu = 0.3, alpha = 1.2, beta = 1.5)
set.seed(1107)
path <- hsim(h, size = 1000)
fit <- hfit(h, inter_arrival = path$inter_arrival, optimizer = "stats")
summary(fit)
```

[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) returns
an `hreal` path. The example uses R’s optimizers and returns an
`emhawkes_stats_fit`. By default,
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) uses
the optional `maxLik` package and returns a `maxLik` fit; install it
with `install.packages("maxLik")` to use that provider. The first path
row is the initial state, so `size` includes that row.

Simulation and fitting default to `lambda_component0 = "auto"` and
report the selected initialization. Supported stable constant linear
models use mean intensity components; other models use zero components
and recommend supplying an initial state. These components exclude the
baseline intensity.

## Documentation

- [Get
  started](https://ksublee.github.io/emhawkes/articles/example.html):
  the basic workflow, input data and a first diagnostic.
- [Fitting and
  inference](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html):
  parameter sharing, initial states, fitted models, residuals and
  standard errors.
- [Extended
  examples](https://ksublee.github.io/emhawkes/articles/extended_example.html):
  marked, state-dependent, multi-kernel and flexible residual models.
- [Computation](https://ksublee.github.io/emhawkes/articles/computation.html):
  R/C++ engines, gradients, Jacobians and performance.
- [Reference](https://ksublee.github.io/emhawkes/reference/index.html)
  and [News](https://ksublee.github.io/emhawkes/news/index.html):
  function details and version changes.

The basic vignette is also available in R with
[`vignette("example", package = "emhawkes")`](https://ksublee.github.io/emhawkes/articles/example.md).
The website is built from this repository; the version in its header
identifies the documented package.
