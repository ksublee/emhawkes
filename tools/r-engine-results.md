# R likelihood performance

This change reduces the cost of the R likelihood used by callback models in
`hfit()`. Gradient selection, optimizer settings, parameter mapping and final
Hessian selection are unchanged.

Repeated calls select their named arguments once, then invoke the callback
without rebuilding a context list at every event. Supplied values are forced
before the callback runs, preserving values captured by deferred closures.
Optional callback defaults remain unevaluated unless the callback uses them.

Decay factors and stable integration weights are computed in blocks of at most
65,536 event-component entries. Small inputs and classed numeric inputs retain
the per-event calculation. The cache preserves matrix names and handles zero
intervals, underflowed products and zero components. Ordinary nonnegative
intensities use `log()` directly; invalid proposals and custom Math methods
retain warning suppression.

## Timing

The baseline is the R likelihood at commit `63dd422`, immediately before this
change. It is not the user's older installed release. All other fitting code is
identical. Each case has 2,000 simulated events. Both functions use the package
namespace and are warmed up before three alternating measurements in one R
session. Timing covers the complete default `hfit()` call, including the final
Hessian. Data generation, function replacement and explicit garbage collection
are outside the timed interval.

Median elapsed times:

| Model | Before (s) | After (s) | Time reduction |
| --- | ---: | ---: | ---: |
| User example: univariate counts | 4.89 | 3.53 | 27.8% |
| Bivariate counts and marks | 25.54 | 19.30 | 24.4% |
| Multi-kernel counts | 10.46 | 7.95 | 24.0% |

These measurements use Windows 11, R 4.6.0 and maxLik 1.6-10. Absolute times vary
with machine load. Raw repetitions and package versions are recorded in
[r-engine-results](r-engine-results/).

## Correctness

All repetitions converge with code 0. Estimates, log likelihood, gradient,
Hessian, convergence code and iteration count are exactly identical to the
baseline for all three cases. Therefore the default Hessian-based standard
errors also remain unchanged.

`tests/r-callback-dispatch.R` checks callback argument matching, defaults,
`missing()`, dots, deferred closures and callbacks with no arguments. An
independent multivariate recursion checks the likelihood, intensity paths and
callback order. Further cases cover cache boundaries, unequal decay rates,
kernel names, tiny rates, zero and extreme intervals, custom numeric classes
and log warning suppression.

`R CMD check --no-manual --as-cran` passes with no errors, warnings or notes.
All 15 test files, examples and vignettes pass, including the existing analytic
gradient tests.

## Reproduction

Install the candidate package in a separate library, make that library visible
to R, then run from the repository root:

```sh
Rscript tools/benchmark-r-engine.R 63dd422 2000 3 results-directory
```

The script replaces only the R likelihood binding inside its R process and
restores it on exit. It writes timings, exact comparisons, estimates, fitted
result snapshots and session information. It does not install a package or
modify the user's installed library.
