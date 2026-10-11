# emhawkes

`emhawkes` specifies, simulates, and fits exponential Hawkes models,
including marked, state-dependent, and multi-kernel models.

``` r

library(emhawkes)

h <- new("hspec",
  mu = c(0.2, 0.2),
  alpha = matrix(c(0.5, 0.9, 0.9, 0.5), 2, byrow = TRUE),
  beta = matrix(2.25, 2, 2)
)
set.seed(1107)
path <- hsim(h, size = 1000)
fit <- hfit(h, path$inter_arrival, path$type)
summary(fit)
```

By default, `reduced = TRUE` treats equal starting values **within each
parameter matrix** as one parameter. This example estimates four
parameters: one for `mu`, two for `alpha`, and one for `beta`.
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md)
explains this rule once before optimization.

Use `reduced = FALSE` to estimate the ten matrix entries separately,
even when their starting values are equal. Function components share
parameters by name under either setting. A zero in a numeric matrix is
also a starting value; use a component function to fix selected entries
at zero.

The first data row is the initial state. `inter_arrival`, `type`, and
`mark` start with zero, followed by the observed events. `hsim(size)`
includes this initial row. `N0` and `Nc0` specify initial counts and
mark totals; `Nc0` must be named in
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

Numeric Hawkes specifications use a compiled Rcpp engine for likelihood
evaluation and inference, including multi-kernel models, unequal decay
rates, and numeric mark effects. Function-valued parameters,
user-defined impact or mark densities, residual-distribution models, and
numeric inputs with custom classes use the R implementation. Simulation
continues to use R.

For comparison or troubleshooting, set `options(emhawkes.engine = "r")`
to use the R calculation, or `"auto"` to restore automatic selection.
The development setting `"cpp"` requires a supported numeric
specification and reports unsupported features. Installing from source
requires a C++ compiler and Rcpp; Windows source installations require
Rtools.

See the [basic
vignette](https://ksublee.github.io/emhawkes/articles/example.html) and
[extended
examples](https://ksublee.github.io/emhawkes/articles/extended_example.html).

The website is generated with `Rscript tools/build-docs.R` after
installing the current package and `pkgdown`. The script includes the
extended article from `vignettes_ext`, which is excluded from the CRAN
package. The [GitHub Markdown
example](https://github.com/ksublee/emhawkes/blob/master/github_document/example.md)
uses the same extended source.
