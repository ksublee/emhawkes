# emhawkes

`emhawkes` specifies, simulates, and fits exponential Hawkes models, including marked, state-dependent, and multi-kernel models.

```r
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

Use `fitted_hspec(fit)` to restore the estimated model and
`fit_initial_state(fit)` to retrieve the initial-state arguments:

```r
estimated <- fitted_hspec(fit)
initial <- fit_initial_state(fit)
inferred <- do.call(logLik, c(list(object = estimated,
  inter_arrival = path$inter_arrival, type = path$type, mark = path$mark,
  infer = TRUE), initial))
```

The result remains a `maxLik` object, so `coef()`, `summary()`, and `vcov()` keep
their existing behavior. The saved parameter map preserves shared parameters and
fixed constants in component functions. Simulation generators are retained, and
their named parameters are updated when they match estimated coefficient names.
Automatic initial intensity remains `NULL` and is recalculated for the restored
model. Initial counts are retained. A custom objective is not incorporated into
the restored model.

Save the fit with `saveRDS()` to retain this information between sessions.
Keep the observed data separately and supply them again for inference. The stored
optimization objective can retain data in its function environment.
The saved initial state describes the start of the observation window. Pass it
to `hsim()` to start a new path from the estimated model. Continuing after the
last observed event requires that event's state. Fits created without the new
metadata need to be refitted before using the reconstruction helpers.

By default, `hfit()` uses the final Hessian and standard errors from `maxLik`.
It does not run the additional stabilization checks. If numerical Hawkes
standard errors contain `NaN`, a short message suggests trying
`final_hessian = "auto"`. This option selects difference steps for stabilization
when requested. It does not automatically retry the default calculation.
Request `final_hessian = "auto"` to check numerical Hawkes curvature with scaled central
differences and Richardson extrapolation. This preserves the optimizer's
coefficients and log-likelihood. `summary()` and `vcov()` then use the checked
covariance; free-parameter standard errors are `NA` when stable inversion is
unavailable. Inspect `fit_diagnostics(fit)` for the requested checks.
Supplied Hessians, disabled final curvature, BHHH, custom objectives, and residual
models retain their calculation choices.

A stable Hessian does not guarantee reliable Wald intervals. Nearby kernel decay
rates, weak excitation, boundaries, short samples, and nonsmooth callbacks can
still make inference unreliable. The [Hawkes standard-error report](https://github.com/ksublee/emhawkes/blob/main/tools/hawkes-standard-errors.md)
compares unchanged estimates and corrected standard errors, and discusses profile
likelihood and bootstrap. The [earlier diagnostic report](https://github.com/ksublee/emhawkes/blob/main/tools/fit-diagnostics.md)
records the original Hawkes and FRPP issues.

By default, `reduced = TRUE` treats equal starting values **within each parameter matrix** as one parameter. This example estimates four parameters: one for `mu`, two for `alpha`, and one for `beta`. `hfit()` explains this rule once before optimization.

Use `reduced = FALSE` to estimate the ten matrix entries separately, even when their starting values are equal. Function components share parameters by name under either setting. A zero in a numeric matrix is also a starting value; use a component function to fix selected entries at zero.

The first data row is the initial state. `inter_arrival`, `type`, and `mark` start with zero, followed by the observed events. `hsim(size)` includes this initial row. `N0` and `Nc0` specify initial counts and mark totals; `Nc0` must be named in `hfit()`.

The likelihood ends at the last supplied event. `N0` and `Nc0` accept a scalar
or one value per event type. Supplied `N` and `Nc` matrices must agree with the
event data and initial counts. Explicit `type_col_map` requires initial kernel
components for simulation and the built-in likelihood.

`logLik(h, ..., infer = TRUE)` returns a list with the likelihood, intensity
paths, and residual vectors. Use it in place of the deprecated `infer_lambda()`
and `residual_process()` helpers. Simulation returns an `hreal` object.

Numeric Hawkes specifications use a compiled Rcpp engine for likelihood evaluation and inference, including multi-kernel models, unequal decay rates, and numeric mark effects. Function-valued parameters, user-defined impact or mark densities, residual-distribution models, and numeric inputs with custom classes use the R implementation. Simulation continues to use R.

`hfit()` uses this engine for its built-in likelihood. Gradient and Hessian
calculations remain numerical unless functions are supplied through `grad` and
`hess`. The Rcpp engine does not change the parameter-sharing rules.

For comparison or troubleshooting, set `options(emhawkes.engine = "r")` to use the R calculation, or `"auto"` to restore automatic selection. The development setting `"cpp"` requires a supported numeric specification and reports unsupported features. Installing from source requires a C++ compiler and Rcpp; Windows source installations require Rtools.

See the [basic vignette](https://ksublee.github.io/emhawkes/articles/example.html) and [extended examples](https://ksublee.github.io/emhawkes/articles/extended_example.html).

The website is generated with `Rscript tools/build-docs.R` after installing the current package and `pkgdown`. The script includes the extended article from `vignettes_ext`, which is excluded from the CRAN package. The [GitHub Markdown example](https://github.com/ksublee/emhawkes/blob/main/github_document/example.md) uses the same extended source.
