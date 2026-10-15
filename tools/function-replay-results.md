# Parameter-only standalone likelihood and inference

> Stage 1 measurement snapshot. Later observed-impact support is described in
> [the stage 2 report](observed-replay-results.md) and current `logLik()` help.

Stage 1 extends automatic standalone `logLik()` and `infer = TRUE` to `mu`,
`alpha`, `beta` and `eta` functions whose sole argument is `param` and whose
values are plain numeric matrices. R evaluates the defaults using the existing
`setting()` rule; C++ computes the likelihood, intensities, integrals and residuals.
The fitting objective and optimizer are unchanged. Other callback signatures,
impacts, mark densities, residual likelihoods and classed inputs retain R.
`cpp` requires support; `r` selects the reference calculation.

## Accuracy

All 32 direct engine checks were exactly identical to the R reference, with
maximum absolute difference zero. Each native check verified that the C++ core
was called once, including inference. Every timed result also matched its R
reference. See [accuracy.csv](function-replay-results/accuracy.csv).

The new regression script checks an independent Poisson likelihood and residual
formula, unequal decay rates, marked models, non-square and overlapping kernel
mappings, unobserved event types, empty paths, explicit and automatic initial
states, singular initialization, restored fitted defaults, callback errors,
classes, input preservation and RNG preservation. It also verifies that automatic
fallback reuses resolved values instead of evaluating callbacks again. Existing
automatic-initialization evaluation counts are preserved, including the existing
second evaluation by `get_lambda0()`.

## Public-call times

Windows, Intel Core Ultra 7 356H, R 4.6.0, maxLik 1.5-2.2 and Rcpp 1.1.2.
One process pinned to logical CPU 2 (affinity mask 4), three timing rounds with
alternating R/auto order. Each interface is warmed before timing. Loading,
fixture generation, the Poisson fit and explicit GC are excluded. Timing batches
contain 50 calls at 1,000 events or 10 calls at 10,000 events; table entries are
medians of per-call elapsed time in seconds. Engine-option selection and public
observation preparation are included. Native-dispatch instrumentation is removed
before timing. These are likelihood/inference calls, not complete fitting times.

At 10,000 events:

| Model | Likelihood R | Auto | Speedup | Intensity/residual R | Auto | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Univariate | 0.090 | 0.005 | 18.0x | 0.392 | 0.005 | 78.4x |
| Marked, 2 dimensions | 0.125 | 0.010 | 12.5x | 0.442 | 0.011 | 40.2x |
| 2 dimensions, 4 kernel columns | 0.125 | 0.019 | 6.6x | 0.444 | 0.020 | 22.2x |
| Restored Poisson fit | 0.094 | 0.005 | 18.8x | 0.398 | 0.005 | 79.6x |

The shorter native calls are sensitive to timer resolution. These ratios apply
to these fixtures and this machine. The fixtures use seeded exponential event
intervals and random event types/marks; the restored Poisson model is estimated
from its fixture before replay. Parameter-only callbacks remain R functions;
their bodies are not compiled by this feature.

All sizes and individual samples are in
[summary.csv](function-replay-results/summary.csv) and
[timings.csv](function-replay-results/timings.csv). The runtime source is in this
commit, based on `46a52d3`. [Source hashes](function-replay-results/source-hashes.csv)
normalize line endings and omit roxygen comments, so documentation edits do not
change the runtime fingerprints.

## Reproduction

Install this source into an isolated library. From the repository root:

```
Rscript --vanilla tools/benchmark-function-replay.R LIB OUT
```

Optionally set `EMHAWKES_BENCH_AFFINITY=4` on Windows, using a permitted CPU mask.
The output includes accuracy, timing, source hashes and session information.
