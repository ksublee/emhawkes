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

Use `fitted_hspec(fit)` to restore the estimated model and
`fit_initial_state(fit)` to retrieve the initial-state arguments:

``` r

estimated <- fitted_hspec(fit)
initial <- fit_initial_state(fit)
inferred <- do.call(logLik, c(list(object = estimated,
  inter_arrival = path$inter_arrival, type = path$type, mark = path$mark,
  infer = TRUE), initial))
```

The result remains a `maxLik` object, so
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) keep their existing
behavior. The saved parameter map preserves shared parameters and fixed
constants in component functions. Simulation generators are retained,
and their named parameters are updated when they match estimated
coefficient names. Automatic initial intensity remains `NULL` and is
recalculated for the restored model. Initial counts are retained. A
custom objective is not incorporated into the restored model.

Save the fit with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to
retain this information between sessions. Keep the observed data
separately and supply them again for inference. The stored optimization
objective can retain data in its function environment. The saved initial
state describes the start of the observation window. Pass it to
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) to
start a new path from the estimated model. Continuing after the last
observed event requires that event’s state. Fits created without the new
metadata need to be refitted before using the reconstruction helpers.

By default,
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) uses
the final Hessian and standard errors from `maxLik`. It does not run the
additional stabilization checks. If numerical Hawkes standard errors
contain `NaN`, a short message suggests trying
`final_hessian = "stable"`. Request this option to check numerical
Hawkes curvature with scaled central differences and Richardson
extrapolation. This preserves the optimizer’s coefficients and
log-likelihood. [`summary()`](https://rdrr.io/r/base/summary.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) then use the checked
covariance; free-parameter standard errors are `NA` when stable
inversion is unavailable. Inspect `fit_diagnostics(fit)` for the
requested checks. Supplied Hessians, disabled final curvature, BHHH,
custom objectives, and residual models retain their calculation choices.

The fit also records its calculation methods:

``` r

calculation <- fit_diagnostics(fit)$computation
calculation$engine   # Requested engine, actual engines and auto fallback reason.
calculation$gradient # Numeric, analytic or user-supplied score.
calculation$hessian # Requested and applied final-Hessian methods.
```

An auto fit can record both `"cpp"` and `"r"` if it switches during
estimation. The record includes final curvature evaluations and survives
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Older fits without
this record return `recorded = FALSE` and `NA` for unrecorded methods.
The existing curvature checks remain available at the top level.

A stable Hessian does not guarantee reliable Wald intervals. Nearby
kernel decay rates, weak excitation, boundaries, short samples, and
nonsmooth callbacks can still make inference unreliable. The [Hawkes
standard-error
report](https://github.com/ksublee/emhawkes/blob/main/tools/hawkes-standard-errors.md)
compares unchanged estimates and corrected standard errors, and
discusses profile likelihood and bootstrap. The [earlier diagnostic
report](https://github.com/ksublee/emhawkes/blob/main/tools/fit-diagnostics.md)
records the original Hawkes and FRPP issues.

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

The likelihood ends at the last supplied event. `N0` and `Nc0` accept a
scalar or one value per event type. Supplied `N` and `Nc` matrices must
agree with the event data and initial counts. Explicit `type_col_map`
requires initial kernel components for simulation and the built-in
likelihood.

`logLik(h, ..., infer = TRUE)` returns a list with the likelihood,
intensity paths, and residual vectors. Use it in place of the deprecated
[`infer_lambda()`](https://ksublee.github.io/emhawkes/reference/infer_lambda.md)
and
[`residual_process()`](https://ksublee.github.io/emhawkes/reference/residual_process.md)
helpers. Simulation returns an `hreal` object.

With automatic engine selection, numeric Hawkes specifications use a
compiled Rcpp engine for likelihood evaluation and inference, including
multi-kernel models, unequal decay rates, and numeric mark effects.
Standalone likelihood evaluation retains the R implementation for
function-valued parameters and impacts. Mark densities,
residual-distribution models, and numeric inputs with custom classes use
R. Simulation continues to use R. The fitting adapters below also select
C++ automatically for supported callback models.

[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) uses
this engine for its built-in likelihood. Numerical differentiation
remains the default. For numeric Hawkes models, request the Rcpp
analytic gradient:

``` r

fit <- hfit(h, path$inter_arrival, path$type, gradient = "analytic")
```

This option supports shared and fixed parameters, multi-kernel mappings,
numeric mark effects, and derivatives of automatic initial states. It
computes the likelihood and score together. State-dependent baselines
and impacts that use latent intensities, custom likelihoods, mark
densities, residual distributions, equality constraints, and BHHH are
not yet supported. Singular automatic initialization requires explicit
`lambda_component0`. A supplied `grad` takes precedence. A supplied
`hess` is retained; otherwise `maxLik` can differentiate the analytic
gradient for curvature. The optimizer’s path and last digits may differ
from numerical differentiation. The [gradient
report](https://github.com/ksublee/emhawkes/blob/main/tools/gradient-results.md)
describes accuracy, full-fit timings, and the state-dependent extension.

Parameter-only `mu`, `alpha`, `beta`, and `eta` functions can also use
the analytic option when supplied with their exact Jacobians. For
example, this function fixes one excitation entry at zero:

``` r

alpha_fun <- function(param = c(a = 0.2, b = 0.1)) {
  matrix(c(param["a"], 0, param["b"], param["a"]), 2)
}
alpha_jac <- cbind(a = c(1, 0, 0, 1), b = c(0, 0, 1, 0))
h_fun <- new("hspec", mu = h@mu, alpha = alpha_fun, beta = h@beta)
fit_fun <- hfit(h_fun, path$inter_arrival, path$type,
  gradient = "analytic", jacobian = list(alpha = alpha_jac))
```

Jacobian rows follow [`as.vector()`](https://rdrr.io/r/base/vector.html)
of the returned matrix; columns are named with its `param` names. For
nonlinear transformations, supply a `function(param)` that returns the
Jacobian at that parameter vector. Components and Jacobians must be
deterministic functions of `param` with fixed matrix dimensions. The
package does not derive or verify the formulas automatically. A
component with `param = numeric(0)` needs no Jacobian. Shared parameter
names receive the sum of all component contributions. Missing or
malformed Jacobians produce an error. The analytic fit evaluates
component functions in R and the likelihood and scores in Rcpp. With the
default `"auto"` engine, supported function-valued models also use C++
with numerical differentiation, without Jacobians:

``` r

fit_numeric <- hfit(h_fun, path$inter_arrival, path$type)
```

The functions must return plain numeric matrices of constant dimensions
and be deterministic in `param`. Observed-state `impact` functions also
use this numerical path without Jacobians, including count- and
mark-dependent impacts and the Hawkes flocking example. Numerical C++
fitting also supports observed baselines and sequential impacts that
explicitly request current or historical intensity states, including the
bid-ask example. Mark densities and residual distributions remain
unsupported. Supplied `grad` and `hess` are retained. Automatic
selection uses R for unsupported specifications; forcing `"cpp"` reports
unsupported features. Callback errors propagate. Callbacks must be
deterministic with fixed captures; detecting an unsupported state read
through `...` can require a retry in R. Standalone
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
keeps its numeric-only C++ contract; select the R engine when evaluating
a restored function-valued model.

Observed-state `impact` functions can use the analytic option with
`jacobian = list(impact = impact_jac)`. They may depend on observed `N`,
`Nc`, marks, types, and times, and on the resolved `mu`, `alpha`, and
`beta` matrices. Counts at row `n` include that event. The event’s
impact affects subsequent intensities in the kernel columns selected by
its type.

For `jacobian$impact`, `param` is the **full named optimization
vector**. Return a matrix with one row per impact entry in column-major
order and one named column per fitted parameter, including zero columns.
Include derivatives through any core matrices used by the impact. The
impact function keeps receiving its local parameters. See
[`?hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md) for an
example and the callback contract. Callbacks that read latent
intensities remain unsupported by the analytic option. Default fits
retain the R path. For numerical C++ fitting, declare
`lambda_component_n` to receive the current pre-jump component matrix.
Declare `lambda` or `lambda_component` only when full pre-jump histories
are needed. Observed baselines may use counts and marks, return a
dimensions-by-one matrix, and remain constant between events. They
cannot depend on latent intensities. See the [bid-ask
benchmark](https://github.com/ksublee/emhawkes/blob/main/tools/bid-ask-numeric-results.md)
for the state convention, compatibility checks and timings.

See the [count, mark and flocking
benchmark](https://github.com/ksublee/emhawkes/blob/main/tools/impact-numeric-results.md)
for complete-fit timings and accuracy checks.

For comparison or troubleshooting, set `options(emhawkes.engine = "r")`
to use the R calculation, or `"auto"` to restore automatic selection.
The development setting `"cpp"` reports unsupported features instead of
falling back to R. Parameter-only functions and observed-state impacts
can use the numerical or analytic fitting adapters described above.
Installing from source requires a C++ compiler and Rcpp; Windows source
installations require Rtools.

See the [basic
vignette](https://ksublee.github.io/emhawkes/articles/example.html) and
[extended
examples](https://ksublee.github.io/emhawkes/articles/extended_example.html).

The website is generated with `Rscript tools/build-docs.R` after
installing the current package and `pkgdown`. The script includes the
extended article from `vignettes_ext`, which is excluded from the CRAN
package. The [GitHub Markdown
example](https://github.com/ksublee/emhawkes/blob/main/github_document/example.md)
uses the same extended source.
