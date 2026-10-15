# Automatic engine selection for fitting

`hfit()` now selects the numerical C++ likelihood automatically for supported
parameter-only functions, observed baselines, observed impacts and sequential
state-dependent impacts. Numerical differentiation and final-Hessian defaults
are unchanged. Analytic gradients remain explicitly requested.

Unsupported specifications retain R. Declared latent baselines, distributions
and classed numeric inputs are excluded before callback evaluation. Unsupported
initial component return shapes and classed impact outputs also retain R. A latent state read through `...`
is detected by an internal capability condition; after detection, the fit stays
on R. Callback errors are not treated as capability failures. Callbacks must be
deterministic with fixed captured values; capability detection can repeat a call.

Internal capability conditions distinguish unsupported callback outputs from
errors thrown by the callback itself. Native sequential validation returns a
capability result to the R adapter, which performs this distinction without
catching arbitrary native errors.

The public forced `cpp` setting still reports unsupported features; `r` always
uses R. Standalone likelihood evaluation, inference and simulation keep their
existing engine-selection behavior.

## Complete-fit verification

Five seeded fixtures were refitted with `auto`, numerical differentiation and
the default Hessian. Estimates, maximum likelihood, gradient, Hessian, convergence
code and iterations were exactly identical to the previous explicit-C++ results.
See [comparison.csv](auto-engine-results/comparison.csv).

| Model | Events | Auto elapsed seconds |
| --- | ---: | ---: |
| Bid-ask, full history | 999 | 16.55 |
| Bid-ask, current components | 999 | 15.87 |
| User count/mark example | 2,000 | 1.95 |
| Bivariate count/mark example | 2,000 | 5.37 |
| Hawkes flocking | 1,999 | 52.15 |

These are single verification runs on Windows with R 4.6.0. They are not new
speedup estimates. The change chooses the existing C++ algorithm automatically;
previous R/CRAN comparisons remain in the numerical adapter benchmark reports.

`tests/auto-function-engine.R` checks native dispatch, exact auto/forced-C++
results, R agreement, unsupported-model fallback, one-time fallback for implicit
latent reads, and propagation of user callback errors. Existing numerical adapter
tests also verify that auto actually invokes the native likelihood.

## Reproduction

Use the fixture directories from the Bid-ask and observed-impact numerical
benchmarks, including their saved previous explicit-C++ fit results:

```
Rscript --vanilla tools/benchmark-auto-engine.R LIB BID_ASK_FIXTURES IMPACT_FIXTURES OUT
```

The script stops if any of the six fit fields differs from the previous result.

All 19 standalone regression scripts passed. The final source build and
`R CMD check --no-manual --as-cran` completed with `Status: OK`, including
examples, all regression scripts and vignette rebuilds.
