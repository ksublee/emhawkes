# Changelog

## emhawkes 1.0.1.9000

- Add a model-independent
  [`parameter_plan()`](https://ksublee.github.io/emhawkes/reference/parameter_plan.md)
  registration/restoration contract. Keep the existing hspec parameter
  map, sharing rules and saved-fit format.

- Add
  [`exp_kernel_plan()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_plan.md)
  to validate fixed numeric kernel parameters and baseline once for
  repeated deterministic step/wait calculations. Changing components and
  arguments are still checked; plans retain no event state or RNG.

- Add deterministic
  [`exp_kernel_step()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md)
  and
  [`exp_kernel_wait()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md)
  interfaces for event-wise state updates and inversion of supplied
  residual targets. They support external FRPP simulation without
  replaying growing path prefixes; existing Hawkes and legacy FRPP event
  loops remain unchanged.

- Move detailed FRPP examples to the frpp package documentation. Retain
  a short compatibility example and describe the planned functionality
  migration; existing residual-law functions and Hawkes residual
  diagnostics are unchanged.

- Add experimental numeric
  [`exp_kernel_spec()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md),
  [`as_exp_kernel()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md)
  and
  [`exp_kernel_replay()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_replay.md)
  interfaces for shared state calculations. Replay returns kernel-only
  interval integrals and pre-event sums, with optional component
  histories. Numeric marks, mapped kernels and signed states are
  supported; row-wise equal decay rates are required. Existing event
  loops are unchanged.

## emhawkes 1.0.1

CRAN release: 2026-10-18

- Correct the stats-only subprocess test when maxLik is installed in a
  standard or site library. Always check stats fitting and saved-result
  restoration without loading maxLik; check the missing-provider error
  only when maxLik is genuinely unavailable. Verify the package path
  used by fresh-session tests.
- Give callbacks with `...` consistent count paths and intensity
  histories with and without inference output, correcting likelihoods
  that depended on `infer`.
- Support [`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) on
  `hspec` objects. Lazily load maxLik’s methods when saved maxLik fits
  are first printed, summarized or queried; keep their classes and
  preserve the stats optimizer’s independence from maxLik.
- Integrate density-based waiting-time moments in log time, including
  finite truncation limits, to avoid missing narrow densities or
  accepting large relative errors through the integrator’s default
  absolute tolerance. Handle singular residual densities at times or
  cumulative hazards that round to zero during quadrature; diagnose
  underflow when its error bound exceeds tolerance.
- Return empty legacy residual vectors for types with fewer than two
  events. Report unavailable stationary mean-intensity systems in
  [`hvol()`](https://ksublee.github.io/emhawkes/reference/hvol.md)
  explicitly.
- Stabilize the trapezoid-exponential coefficients and log-tail
  probabilities near its small-`a`, `ell = 1` limit.
- Avoid quantized or zero arrival times from subnormal decay
  intermediates in the R and C++ equal-decay simulators. Keep the RNG
  draw order unchanged.
- Add default R CMD check jobs with `_R_CHECK_SUGGESTS_ONLY_=false` and
  maxLik exposed through a site library, alongside the existing
  `--as-cran` jobs. Retain a separate hard-dependencies-only CI job that
  requires maxLik’s absence.

## emhawkes 1.0.0

CRAN release: 2026-10-17

### Initialization and estimation

- Default simulation, likelihood and fitting to
  `lambda_component0 = "auto"` (`NULL` is an alias). Report whether mean
  or zero intensity components were selected. Mean initialization
  applies to supported stable constant linear models, including mapped
  kernels; other models use zero and recommend an explicit initial
  state. Components exclude the baseline intensity.
- Correct the univariate half-mean initialization. Preserve explicit
  initial states, hold the selected fitting rule fixed and differentiate
  mean components. Store the rule and actual components in fit
  diagnostics; recover the actual initial states when replaying older
  saved fits. Automatic initialization can change simulated paths and
  fitted estimates relative to previous releases.
- Add `hfit(optimizer = "stats")`, using R’s
  [`optim()`](https://rdrr.io/r/stats/optim.html) or
  [`constrOptim()`](https://rdrr.io/r/stats/constrOptim.html) without
  maxLik. Support fixed/shared parameters, derivatives, linear
  inequalities and box bounds, with separate result methods and
  fitted-model reconstruction. Make maxLik optional while retaining it
  as the default provider; install it explicitly when using the default.
  Keep numerical score steps within linear constraints and reject
  unsupported controls and equality constraints.
- Accept both `constraint` and `constraints`, preserving positional
  arguments. Reject simultaneous non-NULL values and normalize either
  spelling for maxLik.
- Add
  [`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
  [`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)
  and
  [`fit_diagnostics()`](https://ksublee.github.io/emhawkes/reference/fit_diagnostics.md)
  to restore fitted specifications and inspect optimizer, engine and
  curvature choices.
- Add `final_hessian = "stable"`, with scaled differences,
  analytic-score differentiation where available and step-stability
  checks. Withhold Wald covariance when curvature or constraint checks
  fail. Stable curvature is optional for maxLik and the default for
  stats fits; it does not establish finite-sample confidence-interval
  coverage.

### Computation

- Add registered Rcpp likelihood and intensity/residual inference
  engines for supported numeric, parameter-only and observed-state
  specifications. Preserve callback defaults and event timing;
  unsupported models retain R evaluation. `emhawkes.engine` selects
  automatic, reference R or required C++ likelihoods. Source
  installation now requires Rcpp and a C++ compiler.
- Add opt-in analytic scores for numeric models and supported component
  functions through `gradient = "analytic"` and `jacobian`. Preserve
  shared/fixed parameters and initial-state derivatives; reduce
  observed-baseline Jacobian storage. Numerical optimizer
  differentiation remains the default.
- Add a complete native simulation event loop for supported numeric and
  callback specifications, including multivariate and mapped
  multi-kernel models. Reduce history copying and reuse equal decay
  factors within each row. Preserve R fallback, explicit initial states,
  callback behavior and random-draw order when comparing the engines.
  The likelihood engine option does not select the simulation path.
- Preserve native RNG state around callback and warning handlers,
  including errors, nested simulation and handlers that draw random
  numbers. Add checks for older R/Rcpp versions, garbage collection,
  interruption and memory access.

### Corrections and diagnostics

- Register S4 initialization and display methods through the package
  namespace. Preserve construction and saved-fit restoration in sessions
  where methods is installed but not attached, including isolated
  stats-optimizer sessions.
- Validate observations, initial components, decay rates, mappings,
  counting processes and callback parameters consistently. Correct
  multi-kernel event types, univariate multi-kernel likelihoods and
  constant-model stability checks. Preserve admissible signed components
  and stop when no further event occurs.
- Correct event-side baseline recording and return empty residuals for
  event types with fewer than two events. Use stable exponential
  integrals and log-scale residual distributions. Correct
  trapezoid-exponential generation and add log/tail distribution
  options. Paths can change where event selection or random generation
  was corrected.
- Correct
  [`expected_tau()`](https://ksublee.github.io/emhawkes/reference/expected_tau.md)
  for widely separated time scales and a positive probability of no
  further event. Correct mark handling in
  [`hvol()`](https://ksublee.github.io/emhawkes/reference/hvol.md) and
  reject models outside the bivariate moment formula’s scope.
- Preserve display options when printing paths, honor zero displayed
  rows and retain single-row matrix formatting. Make initialization,
  callback, dimension, early-stopping and unavailable-curvature messages
  actionable.

## emhawkes 0.9.8

CRAN release: 2025-08-26

Fixed some bugs

### Breaking changes

- Flexible residual point process with self/mutually exciting and
  exponential decay model is available.

## emhawkes 0.9.7

CRAN release: 2023-02-02

Fixed some bugs, improved efficiency in this version and several
features have been added.

### Breaking changes

- The `lambda0` argument name used in previous versions has been changed
  to `lambda_component0` in this version. This is to clearly indicate
  the meaning of the argument and to avoid confusion.

- The name of method `vol` is changed to `hvol` and this feature is
  currently experimental.

- The Vignette file contains more examples and explanations.

## emhawkes 0.9.6

Fixed some bugs, improved efficiency in this version and several
features have been added.

### Breaking changes

- Slot `eta` is introduced which represents the constant part of
  `impact`.

- The concepts of `rambda` and `rambda_component` are introduced. They
  are closely related to the right-continuous version of the intensity
  process.

- For inference of intensity and goodness of fit, `infer_lambda` and
  `residual_process` functions are implemented.

- The method `vol` to measure the volatility is introduced and this
  feature is currently experimental.
