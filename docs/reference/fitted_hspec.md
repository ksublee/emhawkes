# Restore an Estimated Hawkes Specification

Reconstruct an `hspec` object from the coefficients and parameter map
saved by
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md). The
original specification and fit are not modified.

## Usage

``` r
fitted_hspec(fit)
```

## Arguments

- fit:

  A result from the current
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

## Value

An `hspec` object with the estimated parameter values.

## Details

The saved map preserves matrix sharing under both `reduced` settings and
named parameters shared between functions. Function bodies and
environments are retained; their `param` defaults are updated. Fixed
constants remain fixed. The original kernel mapping, `rmark`, and
`rresidual` are retained. Named default parameters in these generators
are updated when they match coefficient names; other generator
parameters keep their original defaults.

The fit keeps its `maxLik` classes and methods. Its `emhawkes` attribute
stores a versioned list with `specification`, `parameter_map`,
`initial_state`, `reduced`, and `custom_loglik`. Older fits without this
metadata cannot be reconstructed by this function. Function environments
remain references to the original environments; external state is not
frozen. Use [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) to retain the fit's
metadata between sessions. Supply observations again to
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md) and
keep them separately. The stored optimization objective can retain data
in its function environment.

Restoration does not assert convergence or stationarity. Check the fit's
convergence information and whether the model is suitable for
simulation. A custom `mylogLik` is not part of the restored `hspec`.

## See also

[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md),
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md),
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md),
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md)

## Examples

``` r
h <- new("hspec", mu = 0.3, alpha = 0, beta = 1)
fit <- hfit(h, c(0, 0.2, 0.3, 0.5, 1), lambda_component0 = 0,
            fixed = c(FALSE, TRUE, TRUE))
estimated <- fitted_hspec(fit)
initial <- fit_initial_state(fit)
do.call(logLik, c(list(object = estimated,
                      inter_arrival = c(0, 0.2, 0.3, 0.5, 1)), initial))
#> loglikelihood 
#>     -1.227411 
```
