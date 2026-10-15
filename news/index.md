# Changelog

## emhawkes 0.9.9

- Reduce observed-baseline dispatch overhead in numerical fitting and
  standalone likelihood/inference. Evaluate the unchanged R callback
  from a native loop, preserving default parameters, event indices,
  terminal calls and capability fallback.

- Extend standalone C++ likelihood and intensity/residual inference to
  observed baselines and sequential state-dependent impacts, including
  Bid-ask models. Preserve terminal baseline evaluation, component
  attributes, callback snapshots, and independent function defaults.

- Extend standalone C++ likelihood and intensity/residual inference to
  observed-state impacts. Preserve independent callback defaults and
  permanent R fallback during fitting. Read flocking counts directly by
  row and column to avoid repeated full-column extraction.

- Accelerate standalone
  [`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
  and intensity/residual inference for parameter-only component
  functions, including restored fitted specifications. Preserve default
  evaluation, initialization, fallback and optimizer behavior.

- Add event/dimension/kernel scaling and allocation benchmarks. Use
  current components in the GitHub Bid-ask example and explain
  full-history copying costs; preserve likelihood and optimizer
  behavior.

- Align engine documentation and examples: `auto` accelerates supported
  fits by default, `r` forces the reference likelihood, and `cpp`
  requires native support. Distinguish standalone likelihood support and
  explicit analytic gradients.

- Record actual likelihood engines, auto fallback reasons, optimizer
  gradient sources and final-Hessian methods in
  `fit_diagnostics(fit)$computation`. Preserve the curvature diagnostics
  and saved-fit compatibility; report unrecorded methods without
  inferring them from current options.

- Extend automatic engine selection in
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) to
  supported parameter-only functions, observed baselines, and observed
  or sequential state-dependent impacts. Keep numerical differentiation
  and Hessian defaults. Unsupported specifications use R; forced C++ and
  user callback errors remain explicit. Standalone likelihood evaluation
  and simulation are unchanged.

- Support numerical C++ fitting for observed baselines and sequential
  state-dependent impacts, including the bid-ask price model. Preserve
  pre-jump state timing and history layout, protect retained callback
  snapshots, and keep analytic gradients unchanged. These adapters are
  now included in automatic engine selection.

- Allow observed-state `impact` functions to use the C++ likelihood with
  numerical differentiation, including count- and mark-dependent impacts
  and Hawkes flocking. Require no Jacobians and reject latent-intensity
  dependencies, including reads through `...`. These adapters are now
  included in automatic engine selection.

- Allow parameter-only component functions to use the C++ likelihood in
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) with
  `emhawkes.engine = "cpp"` and the default numerical differentiation.
  Require no Jacobians, retain supplied derivatives and the default
  Hessian. Supported fitting adapters now also use C++ with automatic
  engine selection.

- Support analytic gradients for observed-state `impact` functions with
  `jacobian$impact`. Preserve event timing, count paths, kernel
  mappings, and parameter sharing. Use reverse event-impact scores in
  Rcpp and reject latent-intensity dependencies explicitly. Keep default
  fitting unchanged.

- Extend the analytic gradient to parameter-only `mu`, `alpha`, `beta`,
  and `eta` functions through `hfit(jacobian = ...)`. Accept constant or
  parameter-dependent Jacobians, preserve fixed entries and shared
  names, and validate dimensions and parameter names. Default fitting
  remains unchanged.

- Add optional `hfit(gradient = "analytic")` for numeric Hawkes models.
  Compute likelihood and score together in Rcpp, preserve parameter
  sharing and fixed parameters, and differentiate automatic initial
  states. Keep numerical differentiation as the default and retain
  supplied derivatives.

- Name the optional final-Hessian stabilization
  `final_hessian = "stable"`. Keep `"maxLik"` as the default.

- Suggest optional `final_hessian = "stable"` when default numerical
  Hawkes standard errors contain `NaN`. Preserve the fit and respect
  disabled, supplied, and BHHH curvature. Make the disabled-Hessian test
  compatible with both NULL and NA-matrix covariance returns in maxLik.

- Keep the original `maxLik` Hessian and standard-error calculation as
  the default. Run additional Hawkes stabilization and diagnostics only
  when `final_hessian = "stable"` is requested.

- Stabilize final numerical Hawkes Hessians with scaled central
  differences, Richardson extrapolation, and checks across steps.
  Preserve optimizer estimates and likelihoods. Add
  [`fit_diagnostics()`](https://ksublee.github.io/emhawkes/reference/fit_diagnostics.md)
  and optional `final_hessian = "stable"`; withhold Wald standard errors
  when curvature or constraints prevent stable inference.

- Load `maxLik` methods when loading the package so saved fits can be
  restored in a fresh R session. Add an independent-session regression
  test and development tools for fit diagnostics.

- Preserve the original model, parameter map, and initial-state
  arguments in
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md)
  results. Add
  [`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)
  and
  [`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)
  to reuse fitted models with the existing `maxLik` interface.

- Add a registered Rcpp engine for numeric Hawkes likelihoods and
  inference, including unequal decay rates, multi-kernel mappings,
  signed components, and numeric mark effects. Retain the R engine for
  callbacks and residual-distribution models. The `emhawkes.engine`
  option selects automatic, R, or forced C++ evaluation.

- Explain shared matrix parameters once at the start of a reduced fit.
  Update examples and the documentation build to use the current
  interfaces.

- Preserve state-dependent resets and synchronized-model simulation by
  allowing signed components whose aggregate excitation satisfies a
  nonnegative, decreasing decay-rate condition.

- Resolve parameter sharing once per fit while preserving coefficient
  names and order. Validate function parameter defaults and warn about
  conflicting shared starting values.

- Share observation and initial-state validation across estimation,
  likelihood evaluation, and simulation. Check supplied counting
  processes against event data.

- Add the named `Nc0` argument to `hfit`. Select declared callback
  arguments and preserve optional defaults and `...`.

- Expand scalar initial excitation consistently. Warn that simulation’s
  legacy recycling and truncation of other initial values is deprecated.

- Correct event types for multi-kernel simulation and likelihoods for
  univariate multi-kernel models.

- Record state-dependent baselines at the correct side of each event.
  Return empty residuals for types with fewer than two events.

- Stop simulation when no further event occurs, preserving only
  generated rows.

- Use stable exponential integrals and log-scale residual densities and
  survival probabilities.

- Validate decay rates and kernel mappings. Correct stability checks for
  constant multi-kernel Hawkes models.

- Correct trapezoid-exponential random generation and validate its
  parameter domain. Add log and tail options to its distribution
  functions.

- Simulation draws can change for the same seed where event selection or
  random generation was corrected.

## emhawkes 0.9.7

CRAN release: 2023-02-02

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
