# Bid-ask likelihood profiling

The main cost is evaluating R callbacks, rather than the exponential recursion.
The package calculation paths and default engine selection are unchanged.

## Measurement

The fixture is the same seeded Bid-ask example used in
[the numerical likelihood benchmark](bid-ask-numeric-results.md): 999 events,
four dimensions, identical starting parameters and explicit initial components.
The installed C++ implementation is commit `63b3af7`, with R 4.6.0 on Windows.

Each profile evaluates the prepared numerical C++ objective 500 times, alternating
two nearby parameter vectors to avoid its one-point cache. Setup, optimization
and Hessian estimation are excluded. Profiling is sequential across interfaces.

| Sampled stack | Full history | Current components |
| --- | ---: | ---: |
| Baseline path, including R mu | 45.27% | 47.77% |
| Native likelihood, including R impact | 52.22% | 50.07% |
| R impact callback within native likelihood | 41.29% | 38.65% |
| Native call self samples, excluding visible R callbacks | 10.93% | 11.42% |
| matrix(), across both callbacks | 24.82% | 23.07% |

Rows overlap. The impact row is inside the native likelihood row; matrix calls
are inside the callback stacks. They must not be added together. Native self
samples also include callback dispatch and validation, so they do not isolate
exponential arithmetic. Rprof is a sampling profiler, not a native C++ profiler.

The requested sampling interval was 1 ms. Windows timer resolution causes
summaryRprof absolute sampled seconds to differ from wall time; use sample
percentages for attribution. Profiled wall time was 21.30 and 21.15 seconds,
respectively (42.60 and 42.30 ms per objective evaluation). These are starting-point
profile measurements, not replacements for complete-fit benchmark times.

## Improvement assessment

1. Precompute observed price and spread paths once per dataset in a specialized
   Bid-ask adapter. The example recomputes these from Nc in both mu and impact for
   every parameter proposal, although they are independent of fitted parameters.
   Do not cache arbitrary user callback outputs: they can depend on parameters.
2. Keep a compiled callback body and pass changing parameters explicitly, instead
   of rewriting its default arguments for every proposal. This requires preserving
   callback semantics, especially callbacks that call mu or inspect defaults.
3. A specialized native Bid-ask baseline and reset rule would remove the remaining
   per-event R calls and small-matrix allocations. It should be an explicit model
   interface, rather than recognizing user functions from their source text.
4. A state-dependent analytic gradient would reduce repeated likelihood passes.
   It must differentiate reset impacts through the current component state;
   the existing observed-only impact gradient does not cover this case.

Simple byte compilation is not a demonstrated fix. Rebinding parameter defaults
with formals(FUN) <- ... removes bytecode. Compiling each rebound callback adds
compilation work to each proposal. Exploratory timings were not used to claim a
speedup; no production change is justified by those experiments.

## Reproduction

Generate fixtures with tools/benchmark-bid-ask-numeric.R prepare first. Use the
same isolated library and fixture directory as that benchmark:

```r
# Run from the repository root in separate Rscript processes.
# Rscript --vanilla tools/profile-bid-ask.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
# Rscript --vanilla tools/profile-bid-ask-compilation.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
# Rscript --vanilla tools/profile-bid-ask-bound.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
```

Run timing experiments sequentially without competing benchmark processes.
The compilation scripts compare all returned likelihood values exactly, not
merely their elapsed times. Raw profiles are ignored; sampled-stack CSV files
and the measured wall times are included in bid-ask-profile-results/.
