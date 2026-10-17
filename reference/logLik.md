# Compute the Log-Likelihood Function

Evaluate a specified model's likelihood, with optional intensity and
residual paths.

## Usage

``` r
# S4 method for class 'hspec'
logLik(
  object,
  inter_arrival,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  N0 = NULL,
  Nc0 = NULL,
  lambda_component0 = "auto",
  infer = FALSE,
  ...
)
```

## Arguments

- object:

  An
  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
  object containing parameter values for computing the log-likelihood.

- inter_arrival:

  A vector of inter-arrival times for events across all dimensions,
  starting with zero.

- type:

  A vector indicating the dimensions, represented by numbers (1, 2, 3,
  etc.), starting with zero.

- mark:

  A vector of mark (jump) sizes, starting with zero.

- N:

  Optional cumulative event-count matrix, with one row per input row and
  one column per event type. Each row contains counts after that event.

- Nc:

  Optional cumulative mark-total matrix with the same shape as `N`.

- N0:

  Initial event counts: a finite scalar or one value per event type. If
  omitted, use the first row of supplied `N`, or zero when `N` is
  absent.

- Nc0:

  Initial mark totals: a finite scalar or one value per event type. If
  omitted, use the first row of supplied `Nc`, or zero when `Nc` is
  absent.

- lambda_component0:

  Initial intensity components \\\lambda\_{ij}(0)\\, excluding the
  baseline `mu`. A finite scalar is expanded to every kernel component.
  A matrix must match the shape of `alpha` and `beta`. `"auto"` (also
  `NULL`) selects unconditional component means for supported stable
  constant numeric linear models, and zero otherwise. `"mean"` requires
  supported finite means; `"zero"` fixes zero components. The selected
  rule is reported once. Mean initialization excludes nonzero `eta`,
  impact callbacks and custom residual distributions. Function-valued
  core parameters use zero in `"auto"`.

- infer:

  Logical; if `TRUE`, return intensity paths and residual processes
  together with the likelihood.

- ...:

  Additional arguments passed to or from other methods.

## Value

With `infer = FALSE`, a length-one numeric vector named `loglikelihood`.
With `infer = TRUE`, a list containing:

- `loglikelihood`: the same likelihood value.

- `lambda` and `rambda`: intensity matrices immediately before and after
  events.

- `lambda_component` and `rambda_component`: the corresponding kernel
  components, with columns ordered by intensity row, then kernel column.

- `res_process1`, ..., `res_processD`: residuals between successive
  events of each type. Types with fewer than two events have empty
  residual vectors.

Intensity matrices include the initial row. For residual-distribution
models, `res_process` values use the specified distribution and
`eres_process1`, ..., `eres_processD` give their unit-exponential
transforms.

## Details

The first row represents the initial state and is not an event. Omit
`type` for a univariate model and `mark` for unit marks. The observation
vectors have the same length. The likelihood ends at the last supplied
event. Counting-process rows contain values after the corresponding
event. Supplied `N` and `Nc` matrices must agree with `type`, `mark`,
and any supplied initial counts. Event order is preserved, including
zero intervals. Component functions receive the arguments they declare;
functions with `...` receive the available event context. Intensity
histories are `NULL` before initialization. Baseline functions use
`n - 1` to access the preceding event. If initial intensity components
are omitted, they are determined from the supplied model. For a fitted
model, use
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)
with the arguments returned by
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)
to preserve the fit's initialization setting.

With the default automatic selection, standalone `logLik()` uses the
Rcpp engine for numeric Hawkes specifications and parameter-only
component functions. These functions have only `param` as an argument
and return plain numeric matrices. Their defaults are evaluated as in
the R engine; no Jacobian is required. This supports likelihood values
and, with `infer = TRUE`, intensity paths and residuals. Captured inputs
must remain fixed and callbacks must be deterministic. Observed-state
`impact` functions may also use `n`, `N`, `Nc`, marks, types, times and
resolved core matrices. They must return finite plain numeric matrices
matching the kernel dimensions. Each function retains its own default
parameters; standalone evaluation does not apply fitting's sharing
rules. Observed baseline functions may use event context, excluding
latent intensities, and return finite nonnegative dimensions-by-one
matrices. Sequential impacts may explicitly request `lambda`,
`lambda_component`, or `lambda_component_n`. They receive the state
before the jump. The last post-event intensity uses the next baseline,
as in R. Mark densities, residual distributions, unsupported function
signatures, and classed numeric inputs use the R engine. `rmark` and
`rresidual` are simulation generators and do not by themselves prevent
numeric likelihood acceleration. Set `options(emhawkes.engine = "r")`
for the R calculation, `"auto"` for automatic selection, or `"cpp"` to
require a supported specification. Setting `"cpp"` is not needed for
automatic acceleration. The engine supports unequal decay rates within a
row and multi-kernel mappings. During fitting,
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) also
supports C++ adapters for the callback models documented there. Engine
selection does not change simulation or enable analytic gradients.
Automatic capability fallback may evaluate impact callbacks again. User
errors propagate; deterministic callbacks with fixed captured inputs are
required.

## See also

[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md),
[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md),
[`hfit,hspec-method`](https://ksublee.github.io/emhawkes/reference/hfit.md)
