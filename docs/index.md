# emhawkes

`emhawkes` specifies, simulates, and fits exponential Hawkes models,
including marked, state-dependent, and multi-kernel models.

## Installation

Install the CRAN release:

[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"emhawkes"``)`

For the development version and the features described in the guides:

`# install.packages("remotes")`` ``remotes``::``install_github``(``"ksublee/emhawkes"``)`

Source installation requires a C++ compiler and Rcpp; on Windows,
install Rtools.

## Quick start

[`library`](https://rdrr.io/r/base/library.html)`(`[`emhawkes`](https://github.com/ksublee/emhawkes)`)`` `` ``h`` ``<-`` ``new``(``"hspec"``, mu ``=`` ``0.3``, alpha ``=`` ``1.2``, beta ``=`` ``1.5``)`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1107``)`` ``path`` ``<-`` `[`hsim`](https://ksublee.github.io/emhawkes/reference/hsim.md)`(``h``, size ``=`` ``1000``)`` ``fit`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h``, inter_arrival ``=`` ``path``$``inter_arrival``)`` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit``)`

[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) returns
an `hreal` path;
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) returns
a `maxLik` fit. The first path row is the initial state, so `size`
includes that row.

## Documentation

- [Get
  started](https://ksublee.github.io/emhawkes/articles/example.html):
  the basic workflow, input data and a first diagnostic.
- [Extended
  examples](https://ksublee.github.io/emhawkes/articles/extended_example.html):
  marked, state-dependent, multi-kernel and flexible residual models.
- [Fitting and
  inference](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html):
  parameter sharing, initial states, fitted models, residuals and
  standard errors.
- [Computation](https://ksublee.github.io/emhawkes/articles/computation.html):
  R/C++ engines, gradients, Jacobians and performance.
- [Reference](https://ksublee.github.io/emhawkes/reference/index.html)
  and [News](https://ksublee.github.io/emhawkes/news/index.html):
  function details and version changes.

The basic vignette is also available in R with
[`vignette("example", package = "emhawkes")`](https://ksublee.github.io/emhawkes/articles/example.md).
The website is built from this repository; the version in its header
identifies the documented package.
