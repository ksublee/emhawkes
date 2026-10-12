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
  lambda_component0 = NULL,
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

  Initial values for the lambda component \\\lambda\_{ij}\\. A finite
  scalar is expanded to every kernel component. A matrix must match the
  shape of `alpha` and `beta`. Supply this argument when `type_col_map`
  is set.

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

Numeric Hawkes specifications use the Rcpp engine. Function-valued
parameters, impact functions, mark densities, residual distributions,
and classed numeric inputs use the R engine. `rmark` and `rresidual` are
simulation generators and do not by themselves prevent numeric
likelihood acceleration. Set `options(emhawkes.engine = "r")` for the R
calculation, `"auto"` for automatic selection, or `"cpp"` to require a
supported numeric specification. The engine supports unequal decay rates
within a row and multi-kernel mappings.

## See also

[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md),
[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md),
[`hfit,hspec-method`](https://ksublee.github.io/emhawkes/reference/hfit.md)
