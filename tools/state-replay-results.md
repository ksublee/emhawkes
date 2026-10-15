# Observed-baseline and sequential-state replay

Stage 3 extends standalone `logLik()` and `infer = TRUE` to observed baselines
and sequential impacts. R evaluates the callbacks; C++ performs the event loop,
component integration, intensity recording and residual accumulation. The default
`auto` selects this path. `r` selects the reference; `cpp` requires support.

Each callback keeps its own parameter defaults. Baselines return finite,
nonnegative dimensions-by-one matrices and cannot read latent intensities.
Impacts explicitly request `lambda`, `lambda_component`, or
`lambda_component_n` to read the pre-jump state. The baseline is evaluated for
the interval after the final event as in R, even when its value is unused by
the likelihood. In inference that value determines the last post-event intensity.

Unsupported shapes, classes and state reads fall back to R under `auto`.
Ordinary callback errors propagate. Callbacks must be deterministic with fixed
captured inputs; a capability fallback may repeat calls. Fitting's prepared
objectives, numerical differentiation and Hessian defaults are unchanged.

## Accuracy

All 24 native checks returned exactly the same complete output as R: maximum
absolute difference zero. The checks cover 1,000 and 10,000 events, observed
baselines, and both Bid-ask interfaces, with and without inference. Each check
confirmed one native call. All timed outputs also matched exactly. The full-history
and current-component Bid-ask callbacks produced exactly equal R results.
See [accuracy.csv](state-replay-results/accuracy.csv).

The regression suite also checks an independent varying-baseline likelihood and
residual formula, unequal decay rates, overlapping/reordered kernel mappings,
marks, initial counts, independent same-name defaults, automatic initialization,
empty paths, terminal baseline calls, named/attributed kernel matrices, callback
errors and nonfinite unselected impact entries. Retained lazy snapshots, garbage
collection and local callback modifications preserve the R behavior. Inference
reuses its output histories rather than allocating a second pair for callbacks.

A complete Bid-ask fit was compared with the preceding commit `cb7d116` in
separate installed libraries. Its estimates, maximum, gradient, Hessian, return
code and iteration count were exactly identical. Both reached log-likelihood
-1810.892 with code 0 and 41 iterations. Replaying `fitted_hspec()` with
`fit_initial_state()` matched that maximum and all R inference fields.
See [fit-accuracy.csv](state-replay-results/fit-accuracy.csv).

## Time

Windows, R 4.6.0, one CPU affinity mask of 4. Each entry is the median of three
rounds, using five calls per round at 1,000 events and two at 10,000. The same
installed source supplied both engines. Setup, simulation, loading, warmup and
memory profiling are excluded. These are standalone evaluations, not fitting times.

At 10,000 events:

| Model | Likelihood R / auto (s) | Speedup | Inference R / auto (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| Observed baseline, 2 dimensions | 0.250 / 0.150 | 1.67x | 0.580 / 0.145 | 4.00x |
| Bid-ask, full history | 0.815 / 0.600 | 1.36x | 1.110 / 0.595 | 1.87x |
| Bid-ask, current components | 0.690 / 0.560 | 1.23x | 1.065 / 0.555 | 1.92x |

Bid-ask remains limited by R callback execution. This change accelerates replay;
it does not imply a comparable additional improvement in full fitting, which
already used a prepared native likelihood. Small differences between inference
and likelihood-only timings reflect measurement and garbage-collection variation.
See [timings.csv](state-replay-results/timings.csv) and
[summary.csv](state-replay-results/summary.csv).

## Allocation

One separately profiled call per configuration used `Rprofmem()` after warmup.
At 10,000 events:

| Model | Mode | R allocations (MiB) | Auto allocations (MiB) |
| --- | --- | ---: | ---: |
| Bid-ask, full history | Likelihood | 9.54 | 4.81 |
| Bid-ask, current components | Likelihood | 8.32 | 3.59 |
| Bid-ask, either interface | Inference | 13.19 | 9.77 |

Recorded allocation bytes grew about tenfold between 1,000 and 10,000 events
for these callbacks, which do not retain histories. Inference returns complete
histories for either interface, explaining their equal allocation totals.
Callbacks that retain history snapshots can still cause repeated full copies;
the current-state interface avoids this when only current state is needed.

`Rprofmem()` does not include native `std::vector` storage or every small R heap
allocation. Allocation totals are cumulative bytes, not peak resident memory.
The CSV also records the R heap high-water mark after resetting GC maxima; the
same process retains reference outputs and fixtures, so that column is not an
isolated per-call memory estimate. See [memory.csv](state-replay-results/memory.csv).

## Reproduce

Install the working source into an isolated library, then run from the repository
root. On this Windows machine `EMHAWKES_BENCH_AFFINITY=4` selected one CPU;
unset it for the scheduler default.

```text
Rscript --vanilla tools/benchmark-state-replay.R LIB OUT
```

The script extracts the Bid-ask specification from the extended vignette, uses
seed 63001 for one 10,000-event simulation, and takes common prefixes for both
interfaces. The observed-baseline fixture uses seed 63002. It writes accuracy,
raw/median timings, allocation, source hashes and session information. Textual
source fingerprints were refreshed after comment edits; installed executable
bodies were checked against the final R source.

For the previous/current full-fit comparison, prepare the fixture with
`tools/benchmark-bid-ask-numeric.R prepare LIB OUT/prepare`, then run its `worker`
mode with the preceding library (`old`) and candidate library (`cpp`), case
`bid_ask_current`, one repetition and that same fixture directory. Compare the
six saved result fields with `identical()`.
