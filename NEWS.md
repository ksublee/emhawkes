# emhawkes 1.0.2

* Match native simulation and likelihood accumulation to the precision used by
  the running R build, including R configured with --disable-long-double.
  Preserve event/RNG order and callback state consistency without changing
  existing numerical test thresholds. Retain extended precision for derivative
  calculations that do not reproduce an R accumulation.
* Add a no-long-double CI check alongside the standard platform checks.

# emhawkes 1.0.1

* Correct the stats-only subprocess test when maxLik is installed in a standard
  or site library. Always check stats fitting and saved-result restoration
  without loading maxLik; check the missing-provider error only when maxLik
  is genuinely unavailable. Verify the package path used by fresh-session tests.
* Give callbacks with `...` consistent count paths and intensity histories with
  and without inference output, correcting likelihoods that depended on `infer`.
* Support `stats::logLik()` on `hspec` objects. Lazily load maxLik's methods when
  saved maxLik fits are first printed, summarized or queried; keep their classes
  and preserve the stats optimizer's independence from maxLik.
* Integrate density-based waiting-time moments in log time, including finite
  truncation limits, to avoid missing narrow densities or accepting large
  relative errors through the integrator's default absolute tolerance. Handle
  singular residual densities at times or cumulative hazards that round to zero
  during quadrature; diagnose underflow when its error bound exceeds tolerance.
* Return empty legacy residual vectors for types with fewer than two events.
  Report unavailable stationary mean-intensity systems in `hvol()` explicitly.
* Stabilize the trapezoid-exponential coefficients and log-tail probabilities
  near its small-`a`, `ell = 1` limit.
* Avoid quantized or zero arrival times from subnormal decay intermediates in
  the R and C++ equal-decay simulators. Keep the RNG draw order unchanged.
* Add default R CMD check jobs with `_R_CHECK_SUGGESTS_ONLY_=false` and maxLik
  exposed through a site library, alongside the existing `--as-cran` jobs.
  Retain a separate hard-dependencies-only CI job that requires maxLik's absence.

# emhawkes 1.0.0

## Initialization and estimation

* Default simulation, likelihood and fitting to `lambda_component0 = "auto"`
  (`NULL` is an alias). Report whether mean or zero intensity components were
  selected. Mean initialization applies to supported stable constant linear
  models, including mapped kernels; other models use zero and recommend an
  explicit initial state. Components exclude the baseline intensity.
* Correct the univariate half-mean initialization. Preserve explicit initial
  states, hold the selected fitting rule fixed and differentiate mean components.
  Store the rule and actual components in fit diagnostics; recover the actual
  initial states when replaying older saved fits. Automatic initialization can
  change simulated paths and fitted estimates relative to previous releases.
* Add `hfit(optimizer = "stats")`, using R's `optim()` or `constrOptim()` without
  maxLik. Support fixed/shared parameters, derivatives, linear inequalities and
  box bounds, with separate result methods and fitted-model reconstruction.
  Make maxLik optional while retaining it as the default provider; install it
  explicitly when using the default. Keep numerical score steps within linear
  constraints and reject unsupported controls and equality constraints.
* Accept both `constraint` and `constraints`, preserving positional arguments.
  Reject simultaneous non-NULL values and normalize either spelling for maxLik.
* Add `fitted_hspec()`, `fit_initial_state()` and `fit_diagnostics()` to restore
  fitted specifications and inspect optimizer, engine and curvature choices.
* Add `final_hessian = "stable"`, with scaled differences, analytic-score
  differentiation where available and step-stability checks. Withhold Wald
  covariance when curvature or constraint checks fail. Stable curvature is
  optional for maxLik and the default for stats fits; it does not establish
  finite-sample confidence-interval coverage.

## Computation

* Add registered Rcpp likelihood and intensity/residual inference engines for
  supported numeric, parameter-only and observed-state specifications. Preserve
  callback defaults and event timing; unsupported models retain R evaluation.
  `emhawkes.engine` selects automatic, reference R or required C++ likelihoods.
  Source installation now requires Rcpp and a C++ compiler.
* Add opt-in analytic scores for numeric models and supported component functions
  through `gradient = "analytic"` and `jacobian`. Preserve shared/fixed parameters
  and initial-state derivatives; reduce observed-baseline Jacobian storage.
  Numerical optimizer differentiation remains the default.
* Add a complete native simulation event loop for supported numeric and callback
  specifications, including multivariate and mapped multi-kernel models. Reduce
  history copying and reuse equal decay factors within each row. Preserve R
  fallback, explicit initial states, callback behavior and random-draw order
  when comparing the engines. The likelihood engine option does not select the
  simulation path.
* Preserve native RNG state around callback and warning handlers, including
  errors, nested simulation and handlers that draw random numbers. Add checks
  for older R/Rcpp versions, garbage collection, interruption and memory access.

## Corrections and diagnostics

* Register S4 initialization and display methods through the package namespace.
  Preserve construction and saved-fit restoration in sessions where methods
  is installed but not attached, including isolated stats-optimizer sessions.
* Validate observations, initial components, decay rates, mappings, counting
  processes and callback parameters consistently. Correct multi-kernel event
  types, univariate multi-kernel likelihoods and constant-model stability checks.
  Preserve admissible signed components and stop when no further event occurs.
* Correct event-side baseline recording and return empty residuals for event
  types with fewer than two events. Use stable exponential integrals and
  log-scale residual distributions. Correct trapezoid-exponential generation
  and add log/tail distribution options. Paths can change where event selection
  or random generation was corrected.
* Correct `expected_tau()` for widely separated time scales and a positive
  probability of no further event. Correct mark handling in `hvol()` and reject
  models outside the bivariate moment formula's scope.
* Preserve display options when printing paths, honor zero displayed rows and
  retain single-row matrix formatting. Make initialization, callback, dimension,
  early-stopping and unavailable-curvature messages actionable.

# emhawkes 0.9.8

Fixed some bugs

## Breaking changes

* Flexible residual point process with self/mutually exciting and exponential decay model is available.


# emhawkes 0.9.7

Fixed some bugs, improved efficiency in this version and several features have been added.

## Breaking changes

* The `lambda0` argument name used in previous versions has been changed to `lambda_component0` in this version. This is to clearly indicate the meaning of the argument and to avoid confusion.

* The name of method `vol` is changed to `hvol` and this feature is currently experimental.

* The Vignette file contains more examples and explanations.


# emhawkes 0.9.6

Fixed some bugs, improved efficiency in this version and several features have been added.

## Breaking changes

* Slot `eta` is introduced which represents the constant part of `impact`.

* The concepts of `rambda` and `rambda_component` are introduced. They are closely related to the right-continuous version of the intensity process.

* For inference of intensity and goodness of fit, `infer_lambda` and `residual_process` functions are implemented.

* The method `vol` to measure the volatility is introduced and this feature is currently experimental.



