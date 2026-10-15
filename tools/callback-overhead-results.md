# Callback dispatch costs and Bid-ask fits

The baseline path was a measurable part of both fitting and standalone replay.
R previously iterated over events, selected arguments, called the baseline and
checked each result. The new native dispatcher prepares one named call and
collects the baseline path in C++. The user baseline still runs once per event
in R. Its body and the existing parameter-default binding rules are unchanged.

The likelihood kernel, parameter mapping, numerical differentiation and Hessian
defaults are unchanged. Initial and terminal baseline calls, argument names and
order, index types, callback environments, lazy captures and user condition
classes are preserved. Unsupported latent reads keep their typed capability
failure, so automatic fitting can switch permanently to R.

The new dispatcher is shared by numerical fitting and standalone likelihood and
inference. No model-specific function recognition or callback-output cache is
introduced. Functions that return integer matrices or request no runtime
arguments also work. Mark densities and residual-distribution models retain R.

## Validation

The suite includes an independent likelihood/integral reference and comparisons
of complete intensity and residual paths. Additional dispatch tests cover lazy
snapshots, local input modifications, garbage collection, unused defaults,
empty argument lists, nested calls, custom errors, terminal values and capability
failures after initialization. They compare against R and check that input
objects are unchanged. All 24 test scripts passed, and the source package check
with `--no-manual --as-cran` returned `Status: OK`. Existing parameter sharing and differentiation tests run
with the same defaults.

## Results

The previous version is `6bd5eaa`. Both isolated libraries used Windows,
R 4.6.0, Rcpp 1.1.2 and maxLik 1.5-2.2. Elapsed times below are medians of three
rounds, with one CPU affinity mask of 4. Every benchmark mode used the same
999-event observations, initial components and starting parameters.

| Calculation | Full history: old / new (s) | Current components: old / new (s) |
| --- | ---: | ---: |
| Complete `hfit()` | 30.53 / 27.09 | 28.95 / 25.67 |
| Prepared numerical objective | 0.0270 / 0.0235 | 0.0255 / 0.0225 |
| Public `logLik()` | 0.0600 / 0.0455 | 0.0580 / 0.0390 |
| Public inference | 0.0590 / 0.0515 | 0.0570 / 0.0485 |

Complete fitting improved by about 11% in both interfaces. Standalone likelihood
improved by 24% and 33%, respectively; inference improved by 13% and 15%.
These are incremental gains over the preceding native implementation, not a
comparison with CRAN. See [raw timings](callback-overhead-results/timings.csv)
and [medians](callback-overhead-results/summary.csv).

All 108 recorded checks were exact, with maximum absolute difference zero.
They compare every prepared likelihood value, complete public likelihood and
inference outputs, and six fit fields: estimates, maximum, gradient, Hessian,
return code and iteration count. For each interface, every round and version matched its first old
result. Both interfaces reached log-likelihood -1810.892, code 0 and 41
iterations. See [accuracy.csv](callback-overhead-results/accuracy.csv).

Separate profiled complete fits also matched exactly. In the current-component
fit, the baseline path accounted for 48.41% of old stack samples and 42.01% of
new samples. The R `mu_invoke` wrapper disappeared. Matrix creation still
accounted for 28.28% of new samples, inside baseline and impact callbacks.
The core likelihood call, including impact callbacks, accounted for 55.97%.
These are overlapping stack percentages, not separate additive costs. Profiled
absolute timings are excluded from the performance table. See
[profile attribution](callback-overhead-results/profile-attribution.csv) and
the [old](callback-overhead-results/profiles/old/profile-timings.csv) and
[new](callback-overhead-results/profiles/new/profile-timings.csv) profiled runs.

At 10,000 events, recorded cumulative R allocation was unchanged:

| Calculation | Full history (MiB) | Current components (MiB) |
| --- | ---: | ---: |
| Likelihood, either version | 4.81 | 3.59 |
| Inference, either version | 9.77 | 9.77 |

All eight large/prefix output comparisons also matched exactly. This optimization
reduces dispatch work without increasing recorded allocation. It does not remove
the R callback bodies or their small matrix allocations. See
[memory.csv](callback-overhead-results/memory.csv).

Further improvement would require reducing work inside the callbacks, introducing
an explicit specialized model adapter, or reducing repeated likelihood passes
with additional analytic gradients. Rebinding callback defaults remains unchanged;
its cost did not justify a broader interface change in this step.

## Reproduce

Prepare the same 999-event, four-dimensional Bid-ask fixture from the extended
vignette using an isolated previous library:

```text
Rscript --vanilla tools/benchmark-bid-ask-numeric.R prepare OLD_LIB FIXTURES
```

Install the candidate source into a second library. Run each interface in a fresh
process, alternating old/new order across three rounds. On Windows,
`EMHAWKES_BENCH_AFFINITY=4` selects one CPU; unset it for the scheduler default.
Both libraries must use the same R and dependency versions.

```text
Rscript --vanilla tools/benchmark-callback-overhead.R worker LIB FIXTURES OUT TAG ROUND CASE
Rscript --vanilla tools/benchmark-callback-overhead.R compare OUT
```

`TAG` is `old` or `new`, `ROUND` is 1, 2 or 3, and `CASE` is `bid_ask` or
`bid_ask_current`. Public likelihood and inference each use batches of 20 calls.
The prepared objective alternates nearby parameter vectors to avoid its
one-point cache. Each mode is warmed before timing; the fit warms optimizer
setup on a 100-event prefix. Fitting uses the existing default numerical
differentiation and final Hessian. Loading, setup, warmup and GC before the batch
are excluded. GC during a timed call is included. Raw result RDS files stay in
OUT and are used for exact comparisons.

Run profiling and allocation separately from timing workers:

```text
Rscript --vanilla tools/profile-callback-overhead.R LIB FIXTURES PROFILE_OUT bid_ask_current
Rscript --vanilla tools/benchmark-callback-memory.R OLD_LIB FIXTURES MEMORY_OUT old
Rscript --vanilla tools/benchmark-callback-memory.R NEW_LIB FIXTURES MEMORY_OUT new
```

The memory script simulates one 10,000-event path with seed 63001 and uses common
prefixes and both interfaces. It compares the complete returned values exactly.
`Rprofmem()` measures recorded cumulative R allocations, not peak RSS, every
small R allocation, or native vector storage. GC high-water marks include other
objects retained in the worker and are diagnostic only.

Rprof stack percentages overlap. The baseline and impact stacks contain matrix
creation and callback bodies; they cannot be added to those costs. Native `.Call`
self samples include dispatch and validation as well as arithmetic. Windows
sampling resolution also separates sampled seconds from wall time. Profiled wall
times are recorded separately and are not used as the performance benchmark.
