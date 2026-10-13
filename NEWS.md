# emhawkes 0.9.9

* Suggest optional `final_hessian = "auto"` when default numerical Hawkes standard errors contain `NaN`. Preserve the fit and respect disabled, supplied, and BHHH curvature. Make the disabled-Hessian test compatible with both NULL and NA-matrix covariance returns in maxLik.
* Keep the original `maxLik` Hessian and standard-error calculation as the default. Run additional Hawkes stabilization and diagnostics only when `final_hessian = "auto"` is requested.
* Stabilize final numerical Hawkes Hessians with scaled central differences, Richardson extrapolation, and checks across steps. Preserve optimizer estimates and likelihoods. Add `fit_diagnostics()` and optional `final_hessian = "auto"`; withhold Wald standard errors when curvature or constraints prevent stable inference.
* Load `maxLik` methods when loading the package so saved fits can be restored in a fresh R session. Add an independent-session regression test and development tools for fit diagnostics.
* Preserve the original model, parameter map, and initial-state arguments in `hfit()` results. Add `fitted_hspec()` and `fit_initial_state()` to reuse fitted models with the existing `maxLik` interface.
* Add a registered Rcpp engine for numeric Hawkes likelihoods and inference, including unequal decay rates, multi-kernel mappings, signed components, and numeric mark effects. Retain the R engine for callbacks and residual-distribution models. The `emhawkes.engine` option selects automatic, R, or forced C++ evaluation.
* Explain shared matrix parameters once at the start of a reduced fit. Update examples and the documentation build to use the current interfaces.
* Preserve state-dependent resets and synchronized-model simulation by allowing signed components whose aggregate excitation satisfies a nonnegative, decreasing decay-rate condition.
* Resolve parameter sharing once per fit while preserving coefficient names and order. Validate function parameter defaults and warn about conflicting shared starting values.
* Share observation and initial-state validation across estimation, likelihood evaluation, and simulation. Check supplied counting processes against event data.
* Add the named `Nc0` argument to `hfit`. Select declared callback arguments and preserve optional defaults and `...`.
* Expand scalar initial excitation consistently. Warn that simulation's legacy recycling and truncation of other initial values is deprecated.
* Correct event types for multi-kernel simulation and likelihoods for univariate multi-kernel models.
* Record state-dependent baselines at the correct side of each event. Return empty residuals for types with fewer than two events.
* Stop simulation when no further event occurs, preserving only generated rows.
* Use stable exponential integrals and log-scale residual densities and survival probabilities.
* Validate decay rates and kernel mappings. Correct stability checks for constant multi-kernel Hawkes models.
* Correct trapezoid-exponential random generation and validate its parameter domain. Add log and tail options to its distribution functions.
* Simulation draws can change for the same seed where event selection or random generation was corrected.

# emhawkes 0.9.7

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



