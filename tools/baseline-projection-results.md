# Observed-baseline derivative storage

## Change

Observed-baseline analytic scores previously retained a `(D * P) by M` total
Jacobian path, where `M` includes the initial row. Only the derivative at the
event's target dimension and the dimension sum enter that event's baseline score:

`J[type[n], p, n] / intensity[n] - dt[n] * sum(J[, p, n])`.

The native adapter now stores that event derivative and a scaled dimension sum,
with a separate `D by P` Jacobian for initial-state differentiation. The sum uses
a numeric mantissa and an integer binary exponent; the event loop also scales
the interval before multiplication. This avoids overflowing a finite-Jacobian
sum or discarding subnormal inputs merely through storage conversion.

Raw array payload changes from `8 * D * P * M` bytes to
`20 * P * M + 8 * D * P`. The adapter uses the original dense path when its payload
is smaller, including dimensions one and two. Small objects also have list and
matrix overhead. A constant supplied Jacobian still uses one matrix without an
event path.

The full Jacobian callback sequence remains before the impact replay. Callback
counts, indices, parameter-column mapping, retained arguments, automatic initial
states, and cold-score caching keep their contracts. The public `jacobian$mu`
interface and default fitting options do not change. This is score storage
compression, not interleaved callback evaluation or full-history differentiation.
State-Jacobian callback matrices and baseline value paths are unchanged.

## Measurements

The comparison uses baseline `b928eca` and the updated source in separate R
libraries. Each worker is a fresh R process on the same Windows host, pinned to
processor affinity mask 4. Runs are sequential and version order alternates
across three timing rounds. Memory profiling is a separate round. Simulation,
setup, warm-up, RSS queries, and allocation-log parsing are outside timings.

Fixtures come from `benchmark-gradient-scaling.R prepare`: current-component
Markov impact, count-dependent baseline, unequal decays, and explicit initial
components. The large storage-only path repeats observed data to 100,000 events;
it is not a parameter-recovery sample. Full fits use model-generated prefixes:
750 events for the 4-dimension, 8-kernel-column, 12-parameter numerical-gradient
fit, and 1,000 events for the 8-dimension, 16-column, 24-parameter analytic fit.
Both fits use BFGS and opt-in stable curvature with the supplied exact Jacobians.

Rprofmem records cumulative R allocations and the largest block. Process peak
working set includes package loading, observations, preparation and warm-up;
it is not a native-buffer-only peak. Cumulative allocation is not peak memory.

## Results

For the 100,000-event, 8-dimension, 24-parameter storage-only operation:

| Measurement | Before | After |
| --- | ---: | ---: |
| Retained Jacobian objects | 146.49 MiB | 45.78 MiB |
| Process peak working set | 293.33 MiB | 193.15 MiB |
| Cumulative R allocations | 349.56 MiB | 249.23 MiB |
| Largest R allocation | 146.49 MiB | 18.31 MiB |

Retained storage falls **68.7%**; process peak working set falls **34.2%** in this
storage-only operation. At 10,000 events the raw array payload is approximately
**14.65 versus 4.58 MiB**. This does not predict the same decrease for whole fits.
The 10,000-event, 8-dimension value-and-score operation has about 1,399 versus
1,389 MiB of cumulative R allocations, with process peaks about 196 versus
197 MiB. Dense state-Jacobian callbacks still dominate that operation.

Three-run median elapsed seconds:

| Operation | Events | Before | After |
| --- | ---: | ---: | ---: |
| Whole fit: D4 / P12, numeric + stable | 750 | 11.61 | 9.94 |
| Whole fit: D4 / P12, mixed baseline, numeric + stable | 750 | 13.98 | 13.75 |
| Whole fit: D8 / P24, analytic + stable | 1,000 | 53.11 | 53.40 |
| Cold score: D8 / P24 | 10,000 | 1.14 | 1.16 |
| Value + score: D2 / P3, dense-path control | 10,000 | 0.22 | 0.20 |
| Value + score: D4 / P12 | 10,000 | 0.37 | 0.38 |
| Value + score: D4 / P12, mixed baseline | 10,000 | 0.49 | 0.47 |
| Value + score: D8 / P24 | 10,000 | 1.25 | 1.25 |
| Baseline Jacobian construction: D8 / P24 | 100,000 | 1.70 | 1.78 |

Whole-fit acceleration is not established. The basic D4 numerical fit is
particularly variable: before times are 11.61, 11.80 and 6.85 seconds; after
times are 6.65, 12.08 and 9.94 seconds. Its median alone is not a dependable
speedup claim. The larger analytic fit differs by about 0.5%. Construction alone
takes about 5% longer by the medians, with overlapping sample ranges; it avoids
the large retained tensor. The unchanged D2 control also shows timing variation.

The mixed-baseline fixture uses `0.8 * I + 0.2 / D` to mix baseline amplitudes.
This matrix is full rank and preserves the common generating amplitudes, while
each baseline parameter now affects every target dimension. This checks a
nontrivial dimension sum in a full fit, in addition to the sparse baseline
derivatives of the original scaling fixtures.

See [all timing samples](baseline-projection-results/timings.csv),
[medians](baseline-projection-results/summary.csv),
[memory measurements](baseline-projection-results/memory.csv),
[fit outcomes](baseline-projection-results/fits.csv),
[accuracy checks](baseline-projection-results/checks.csv), and
[estimates and standard errors](baseline-projection-results/estimates.csv).

## Accuracy and validation

All three full-fit configurations converge with code 0 and finite stable
standard errors in every round; no warnings are emitted. Estimates, likelihoods,
convergence codes and iteration counts are exactly equal before and after.
There are 48 distinct model-parameter standard errors per version, repeated
across three rounds. Their maximum relative difference is **9.12e-14**.

Maximum coordinate-scaled differences are **3.22e-15** for scores, **8.33e-14**
for Hessians, and **4.86e-15** for covariances. Coordinate scaling divides by
`max(1, abs(reference))`; the SE relative difference instead divides by the
reference SE itself. Fixed-point likelihoods also agree exactly.

Source build and `R CMD check --no-manual` pass with **Status: OK**, including
all **28 test scripts**, examples, and vignette rebuilding. There are no errors,
warnings or notes. The checked runtime files and new regression match the final
working files; the measured runtime and immutable driver match the repository.
See [source verification](baseline-projection-results/source-verification.txt),
[manifest](baseline-projection-results/manifest.txt), and
[session information](baseline-projection-results/session-info.txt).

`tests/baseline-projection.R` independently compares compact and full-path native
scores. It covers changing parameter-column order, signed dense Jacobians,
unequal rates, overlapping nonsquare kernel maps, empty and short paths, zero
intervals, retained lazy indices, malformed Jacobians, shared parameters,
explicit/automatic initial states, and current-state feedback. Public scores are
also checked against two central-difference steps of the R likelihood. Separate
Poisson cases check the integral at subnormal and near-maximum double scales.

Compression combines the baseline's dimension sum before updating the score;
last-bit differences can occur from floating-point summation. It does not alter
likelihood or optimizer settings and does not establish statistical identification
or Wald coverage.

## Reproduce

Prepare the fixtures with `tools/benchmark-gradient-scaling.R prepare LIB FIXTURES`.
For each library, run isolated workers from the repository root:

```sh
Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT n10000-d8-k16-p24-diagonal-observed score 1
Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT n10000-d8-k16-p24-diagonal-observed cold 1
Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT n10000-d8-k16-p24-diagonal-observed path 1
Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT n10000-d4-k8-p12-diagonal-observed fit_numeric 1
Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT n10000-d8-k16-p24-diagonal-observed fit_analytic 1
Rscript --vanilla tools/compare-baseline-projection.R BEFORE_OUT AFTER_OUT RESULTS
```

Use rounds 1 to 3 for timings, alternating version order. Round 0 records memory
for `score`, `cold`, and `path`; do not run timing workers concurrently. The score
study also includes the 2-dimension, 4-column, 3-parameter control and the
4-dimension, 8-column, 12-parameter case. Set `EMHAWKES_BENCH_AFFINITY=4` on Windows
to reproduce the recorded affinity. Private fixtures and RDS results are not
included; published CSV files retain the measurements and comparisons.
