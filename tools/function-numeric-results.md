# Numerical fitting for parameter-only functions

> Historical benchmark: the explicit `"cpp"` setting below records the measured configuration. Current `hfit()` also selects C++ automatically for these supported models with the default `"auto"`; `"cpp"` requires support.

`hfit()` can now use the existing native Hawkes objective with numerical
differentiation when `options(emhawkes.engine = "cpp")` is selected. No
component Jacobians are needed. `mu`, `alpha`, `beta` and `eta` functions must
have exactly the `param` argument and return plain numeric matrices of fixed
dimensions. They must be deterministic, with captured constants held fixed.

The default `"auto"` engine is unchanged. State-dependent baselines, impacts,
mark densities and residual distributions remain outside this numerical path.
The analytic-gradient adapter is retained. Supplied `grad` and `hess` keep their
precedence, as does a custom likelihood. Standalone `logLik()` retains its
numeric-only C++ contract; evaluate restored function models with the R engine.

```r
previous_engine <- options(emhawkes.engine = "cpp")
fit <- hfit(h, inter_arrival, type)  # numerical differentiation; no Jacobian
options(previous_engine)
```

## Timing

All four examples contain parameter-only functions and no `impact`. They are
ordinary Hawkes models with nonlinear transformations, fixed zeros, shared
parameter names, or multiple kernels. They are distinct from the earlier
state-dependent count and mark benchmarks.

Each example has 2,000 events plus the initial row. Each version uses the same
data, starting values, parameter bounds, BFGS settings and default final Hessian.
Times cover complete `hfit()` calls, excluding simulation, loading and explicit
garbage collection. The table reports medians of three rounds after warm-up.
R and C++ runs alternate in one candidate-package session. The short univariate
C++ fit uses batches of 11 full fits; the recorded time is per fit.

| Model | CRAN (s) | Current R (s) | C++ numeric (s) | CRAN / C++ |
| --- | ---: | ---: | ---: | ---: |
| Univariate log-transformed components | 2.86 | 3.44 | 0.14 | 19.9x |
| Bivariate fixed zero | 8.16 | 8.19 | 0.76 | 10.7x |
| Bivariate shared parameter and automatic initial state | 5.83 | 5.75 | 0.56 | 10.4x |
| Bivariate multi-kernel | 14.32 | 14.18 | 2.45 | 5.8x |

The comparison uses Windows 11, R 4.6.0 and maxLik 1.5-2.2. The CRAN package is
0.9.8; the previous development baseline is `a4a638b`. Package versions are
loaded in isolated libraries and fresh processes. Raw times vary, especially
for short runs; the table records this machine's measurements rather than a
guaranteed speedup. Raw repetitions and batch sizes are included.

## Accuracy and compatibility

For every example, the current R and new C++ routes give exactly identical
estimates, log likelihood, gradient, Hessian, convergence code and iteration
count. The current default R fit also matches all six fields from the previous
development package exactly. All fits converge with code 0.

CRAN and C++ likelihoods are also compared at each starting point and at both
packages' estimated parameters. The maximum absolute difference over these
12 comparisons is 4.55e-12.

For three examples, CRAN estimate differences are at most 3.81e-7 and optimized
likelihood differences are below 1e-10. In the multi-kernel example, the small
cross effect approaches its lower parameter bound. CRAN and current R/C++
stop at slightly different values: the largest log-parameter difference is
2.61e-4 and the achieved likelihood difference is 7.40e-6. Evaluating either
estimate in both implementations still gives the same likelihood to 4.55e-12.
The comparison accepts 1e-5 for optimized likelihoods and 1e-3 for coefficients;
fixed-point likelihoods must agree within 1e-8.

Current R and C++ standard errors are exactly identical in these examples.
CRAN's default finite-difference Hessian produces materially different standard
errors and two nonfinite entries in the multi-kernel case. This is an existing
CRAN/development difference, not a change introduced by the numerical C++ route.
The timing results do not establish reliable Wald inference.

The new tests cover a Poisson model with a known optimum and variance, actual
native dispatch, unchanged automatic selection, supplied derivatives, fixed
parameters, saved-model reconstruction, parameter sharing, unequal rates,
marks, kernel mappings, invalid proposals and unsupported functions. Singular
multivariate automatic initialization preserves the existing zero-state
fallback. A nonfinite univariate automatic state returns an invalid proposal.

`R CMD check --no-manual --as-cran` passes with no errors, warnings or notes.
All 16 package test files pass; the final added singular-state assertions also
pass against maxLik 1.6-10. Documentation is updated in `hfit`, README and NEWS.

## Reproduction

Install the candidate, CRAN 0.9.8 and the previous development package in separate
libraries with the same dependency versions. From the repository root:

```sh
Rscript tools/benchmark-function-numeric.R prepare OUTPUT CANDIDATE_LIBRARY
Rscript tools/benchmark-function-numeric.R worker OUTPUT PREVIOUS_LIBRARY previous
Rscript tools/benchmark-function-numeric.R worker OUTPUT CRAN_LIBRARY cran
Rscript tools/benchmark-function-numeric.R worker OUTPUT CANDIDATE_LIBRARY candidate
Rscript tools/benchmark-function-numeric.R compare OUTPUT CANDIDATE_LIBRARY CRAN_LIBRARY
```

An optional case name after the worker kind remeasures that case while preserving
other results. The candidate worker batches fits whose warm-up time is below
0.2 seconds. Full fit snapshots and generated point models remain local RDS
files; CSV comparisons and R session information are recorded in
[function-numeric-results](function-numeric-results/).
