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
`final_hessian = "stable"`. Request this option to check numerical Hawkes
curvature with scaled central differences and Richardson extrapolation.
This preserves the optimizer's
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

`hfit()` uses this engine for its built-in likelihood. Numerical differentiation
remains the default. For numeric Hawkes models, request the Rcpp analytic gradient:

```r
fit <- hfit(h, path$inter_arrival, path$type, gradient = "analytic")
```

This option supports shared and fixed parameters, multi-kernel mappings, numeric
mark effects, and derivatives of automatic initial states. It computes the
likelihood and score together. State-dependent baselines and impacts that use latent intensities,
custom likelihoods, mark densities, residual distributions, equality constraints,
and BHHH are not yet supported. Singular automatic initialization requires explicit
`lambda_component0`. A supplied `grad` takes precedence. A supplied `hess` is
retained; otherwise `maxLik` can differentiate the analytic gradient for curvature.
The optimizer's path and last digits may differ from numerical differentiation.
The [gradient report](https://github.com/ksublee/emhawkes/blob/main/tools/gradient-results.md)
describes accuracy, full-fit timings, and the state-dependent extension.

Parameter-only `mu`, `alpha`, `beta`, and `eta` functions can also use the analytic
option when supplied with their exact Jacobians. For example, this function fixes
one excitation entry at zero:

```r
alpha_fun <- function(param = c(a = 0.2, b = 0.1)) {
  matrix(c(param["a"], 0, param["b"], param["a"]), 2)
}
alpha_jac <- cbind(a = c(1, 0, 0, 1), b = c(0, 0, 1, 0))
h_fun <- new("hspec", mu = h@mu, alpha = alpha_fun, beta = h@beta)
fit_fun <- hfit(h_fun, path$inter_arrival, path$type,
  gradient = "analytic", jacobian = list(alpha = alpha_jac))
```

Jacobian rows follow `as.vector()` of the returned matrix; columns are named with
its `param` names. For nonlinear transformations, supply a `function(param)` that
returns the Jacobian at that parameter vector. Components and Jacobians must be
deterministic functions of `param` with fixed matrix dimensions. The package does
not derive or verify the formulas automatically. A component with
`param = numeric(0)` needs no Jacobian. Shared parameter names receive the sum of
all component contributions. Missing or malformed Jacobians produce an error.
The analytic fit evaluates component functions in R and the likelihood and scores
in Rcpp. Default fits and standalone `logLik()` retain the R implementation.

Observed-state `impact` functions can use the analytic option with
`jacobian = list(impact = impact_jac)`. They may depend on observed `N`, `Nc`,
marks, types, and times, and on the resolved `mu`, `alpha`, and `beta` matrices.
Counts at row `n` include that event. The event's impact affects subsequent
intensities in the kernel columns selected by its type.

For `jacobian$impact`, `param` is the **full named optimization vector**. Return a
matrix with one row per impact entry in column-major order and one named column
per fitted parameter, including zero columns. Include derivatives through any
core matrices used by the impact. The impact function keeps receiving its local
parameters. See `?hfit` for an example and the callback contract. Callbacks that
read latent intensities remain unsupported. Default fits retain the R path.

For comparison or troubleshooting, set `options(emhawkes.engine = "r")` to use the R calculation, or `"auto"` to restore automatic selection. The development setting `"cpp"` reports unsupported features instead of falling back to R. Parameter-only functions use C++ through the analytic adapter described above. Installing from source requires a C++ compiler and Rcpp; Windows source installations require Rtools.

See the [basic vignette](https://ksublee.github.io/emhawkes/articles/example.html) and [extended examples](https://ksublee.github.io/emhawkes/articles/extended_example.html).

The website is generated with `Rscript tools/build-docs.R` after installing the current package and `pkgdown`. The script includes the extended article from `vignettes_ext`, which is excluded from the CRAN package. The [GitHub Markdown example](https://github.com/ksublee/emhawkes/blob/main/github_document/example.md) uses the same extended source.
