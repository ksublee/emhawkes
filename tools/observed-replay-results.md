# Observed-impact standalone likelihood and inference

Stage 2 extends standalone `logLik()` and `infer = TRUE` to observed-state
impacts. Core functions remain parameter-only; impacts may use event index,
counts, cumulative marks, observed times/types/marks and resolved core matrices.
R evaluates each impact; C++ performs the intensity and compensator recursion.
No Jacobian is required and fitting's differentiation defaults are unchanged.

Standalone replay retains each callback's own default parameters. It does not
construct a fitting parameter map, which could merge independent same-name
defaults. Unsupported latent reads through dots and unsupported output shapes,
classes or nonfinite entries select R automatically. Ordinary callback errors
propagate. Capability fallback can evaluate callbacks again; callbacks must be
deterministic with fixed captured inputs. A resolved automatic initial state is
reused on fallback. Fitting's separate prepared adapter keeps its permanent R
fallback behavior; it does not repeatedly try the new standalone adapter.

## Accuracy

All 24 direct native checks were exactly identical to R, including all intensity,
component and residual paths: maximum absolute difference zero. Each auto/cpp
check verified one native call. All timed outputs also matched the R reference.
Original and scalar-index flocking callbacks returned exactly identical R results
at both sizes with and without inference. See
[accuracy.csv](observed-replay-results/accuracy.csv).

The regression script additionally checks independent shot-noise sums for the
likelihood, intensities, compensator and residuals; independent same-name defaults;
unequal rates and overlapping/reordered kernel mappings; N0/Nc0; integer impacts;
empty paths; hidden/declared latent states; nonfinite unused impact entries;
user errors; RNG and model preservation. Full fitting regressions cover engine
recording and permanent fallback.

## Public-call time

Same source and installed library for R and auto, Windows, Intel Core Ultra 7
356H, R 4.6.0, maxLik 1.5-2.2 and Rcpp 1.1.2. One process pinned to logical CPU 2
(mask 4), three rounds with alternating engine order. Loading, simulation,
explicit GC and warm-up are excluded. Public observation preparation, R impact
evaluation and C++ recursion are included. Batches have ten calls at 1,000 events
or three calls at 10,000 events; entries below are median per-call seconds.
These are standalone calls, not full fitting times.

At 10,000 events:

| Model/impact | Likelihood R | Auto | Speedup | Intensity/residual R | Auto | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Counts and marks, 2 dimensions | 0.240 | 0.143 | 1.67x | 0.567 | 0.147 | 3.86x |
| Flocking, full-column extraction | 2.397 | 2.427 | 0.99x | 2.203 | 1.897 | 1.16x |
| Flocking, direct scalar indexing | 0.417 | 0.323 | 1.29x | 0.770 | 0.330 | 2.33x |

The original flocking callback extracts an entire count column for every scalar
read (`N[, "N1"][n]`). Its callback cost dominates the recursion: likelihood-only
times show no clear native improvement. Using `N[n, "N1"]` returns the same scalar
without copying a full column at each event. The vignette and GitHub Markdown
example now use that indexing. The benchmark reconstructs both variants even
when the source example already uses the scalar form.

The engine and indexing improvements are separate: compare R/auto within each
row for the engine effect, and the two flocking rows within one engine for the
indexing effect. Short calls and background/CPU-frequency changes limit timing
precision. These results do not guarantee acceleration for arbitrary callbacks.
Their R bodies remain part of each replay. Impacts are temporarily stored as a
components-by-observations matrix before native recursion.

Individual samples and every tested size are in
[timings.csv](observed-replay-results/timings.csv) and
[summary.csv](observed-replay-results/summary.csv). Runtime source fingerprints
omit roxygen comments and normalize line endings:
[source-hashes.csv](observed-replay-results/source-hashes.csv).

## Reproduction

Install this source into a private library, then run from the repository root:

```
Rscript --vanilla tools/benchmark-function-replay.R LIB OUT observed
```

Optionally set `EMHAWKES_BENCH_AFFINITY=4` on Windows. The helper extracts the
vignette's flocking specification, generates a single seeded path and takes the
same prefixes for both indexing variants. The count/mark fixture uses seeded
exponential intervals and random observed types/marks. Raw fixtures remain local.
