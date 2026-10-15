# Analytic gradients for observed baselines and current-state impacts

`hfit(..., gradient = "analytic")` now supports observed baselines and impacts
that explicitly declare `lambda_component_n`. Numerical differentiation remains
the default. The extension requires the `auto` or `cpp` engine and exact supplied
Jacobians. Full intensity histories, latent-intensity baselines, mark densities,
residual-distribution models, custom likelihoods, equality constraints, and BHHH
remain unsupported by the analytic option.

## Derivative contract

- Parameter-only core components keep their local `param` Jacobians.
- An observed baseline uses `jacobian$mu`: dimensions by all fitted parameters,
  with a column for every fitted name, including zeros. Its callback receives
  the full named parameter vector and the observed context. Include dependencies
  through resolved core components. The baseline value function keeps its local
  parameters.
- Current-component impact uses
  `jacobian$impact = list(param = ..., state = ...)`. The first matrix is impact
  entries by all fitted parameters, **holding the current state fixed**. Include
  direct dependencies through core components. The second is impact entries by
  current-component entries. Both state axes use `as.vector()` column-major order.
  A matrix or callback can supply either entry. The Jacobian callbacks receive
  the full parameter vector; the value function keeps its local parameters.

All Jacobians must be finite plain numeric matrices with fixed dimensions.
Parameter columns may be reordered, including between callback calls. The package
validates shapes and names; it does not derive or certify the supplied formulas.
Callbacks must be deterministic, with fixed captures and valid derivatives along
the evaluated parameter path. Ignored `...` is supported. Forcing unavailable
histories through `...` fails explicitly, including during value-only line searches.

Counts at row `n` include the event. A predictable baseline for that event uses
counts at `n - 1`; the example's impact uses counts at `n`. Its current component
matrix is the decayed state before the jump. Observations have zero parameter
derivative. The existing `type_col_map` selects updated columns.

## Recursion

Let `C` be the previous post-jump component matrix and `S` its full parameter
Jacobian. For an interval `d` with component decay `beta`, the pre-jump state is
`Cminus = C * exp(-beta*d)`. Its sensitivity is
`Sminus = exp(-beta*d) * S - d * Cminus * dbeta` entrywise.
The event score differentiates its log intensity; the compensator score
subtracts the derivative of `mu*d + C*(1-exp(-beta*d))/beta`.
Small decay products use the existing stable integral derivative.

For each selected kernel entry, add the jump derivative
`dalpha + (mark-1)*deta + Jparam + Jstate %*% Sminus`.
All jump derivatives use the same pre-jump state, including cross-component
feedback. Nonselected entries keep `Sminus`. An explicit initial state is
constant. Automatic initialization differentiates the existing stationarity
rule, including `mu(n=1)`; a singular initializer needs `lambda_component0`.
Additional impact does not change that initialization rule.

R resolves components, sharing, and initial derivatives. C++ reuses named R calls
for observed-baseline Jacobians and propagates state sensitivities through the
event loop. Value-only line searches defer Jacobians. Values and gradients are
cached at one parameter vector per fit. Current sensitivity storage is proportional
to components times parameters; no full derivative history is stored. An observed
baseline Jacobian path requires observations times dimensions times parameters.
Dense state-Jacobian multiplication can become costly for large models; this
small Bid-ask benchmark does not establish scaling for every callback.

## Accuracy

`tests/state-gradient.R` compares against the R likelihood and central differences
at two step sizes and a directional derivative. It covers shared and unshared
parameters, nonlinear transformations, observed-only and current-state impacts,
unequal decays, mark effects, overlapping/reordered multi-kernel mappings,
explicit and automatic initialization, zero intervals, empty paths, malformed
Jacobians, hidden history access, callback snapshots, and saved fit diagnostics.
An observed Poisson baseline also has an independent closed-form score check.

The Bid-ask fixture has 999 events plus its initial row, four dimensions, four
kernel columns, and nine fitted parameters. It uses seed 1107 and the current-state
version of the extended example. The same observations, starts, initialization,
optimizer, and default final Hessian are used in both modes.
The numerical path matches pre-change commit `80c09cf` exactly in coefficients,
log-likelihood, gradient, Hessian, return code, and iteration count in every round.

The default-start analytic fit returns code 0 and agrees with the numerical fit:
maximum log-likelihood difference **1.18e-08**, largest absolute parameter
difference **1.09e-05**, largest difference divided by `max(1, abs(reference))`
**4.94e-07**. BFGS uses 38 iterations with analytic derivatives and 41 with
numerical derivatives. Different paths and last digits are expected.
See [accuracy.csv](state-gradient-results/accuracy.csv).

Two additional starts scale every starting parameter by 0.9 and 1.1 and use
BFGS `reltol=1e-12`. All four fits return code 0. Numerical/analytic likelihood
differences are below **9.35e-11**; maximum absolute parameter differences are
**3.61e-05**, or **1.18e-06** after scaling by `max(1, abs(reference))`.
At the default start, a nearby vector, and both fitted vectors, likelihoods match
the R engine within **1e-9**. Central differences at three step sizes agree with
the analytic score within **2.31e-06** after scaling by `max(1, abs(reference))`.
See [multiple starts](state-gradient-results/multiple-starts.csv) and
[gradient checks](state-gradient-results/gradient-checks.csv).

All 25 test scripts pass with maxLik 1.5-2.2 and in `R CMD check --as-cran` with
maxLik 1.6-10. Examples and vignette rebuilding are also checked.


## Full-fit timing

Windows 11, R 4.6.0, maxLik 1.5-2.2, Rcpp 1.1.2. Three isolated rounds alternate
method order. Each worker warms up a short fit and uses the same processor
affinity mask 4. Elapsed times include fit setup, optimization, and final Hessian.
Simulation and compilation are excluded. The table reports medians.

| Method | Seconds | Relative to current numerical fit |
| --- | ---: | ---: |
| Before this change, numerical | 23.77 | 1.00x |
| Current default, numerical | 23.78 | 1.00x |
| Current analytic | 2.30 | 10.34x |

These timings describe this fixture and machine. They are not a guarantee for
other callback costs, dimensions, data, or optimizers. Faster gradients do not
establish identification or reliable Wald standard errors.
Raw [timings](state-gradient-results/timings.csv),
[summary](state-gradient-results/summary.csv), and session records are included.

## Reproduce

Install the pre-change and current packages into separate libraries with the same
dependency versions. From the repository root:

```r
# Rscript tools/benchmark-bid-ask-numeric.R prepare OLD_LIB FIXTURES
# Set EMHAWKES_BENCH_AFFINITY=4 on Windows for the recorded affinity.
# Rscript tools/benchmark-state-gradient.R worker OLD_LIB FIXTURES RESULTS old 1
# Rscript tools/benchmark-state-gradient.R worker NEW_LIB FIXTURES RESULTS numeric 1
# Rscript tools/benchmark-state-gradient.R worker NEW_LIB FIXTURES RESULTS analytic 1
# Repeat rounds 2 and 3, alternating order.
# Rscript tools/benchmark-state-gradient.R compare RESULTS
# Rscript tools/check-state-gradient.R NEW_LIB FIXTURES RESULTS
```

The complete Bid-ask Jacobians are in `tools/state-gradient-cases.R` and the
extended vignette. Requesting the analytic option without the required Jacobians
produces an error. Public defaults and numerical fitting remain unchanged.
