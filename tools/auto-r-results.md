# Representative auto and R fitting comparison

Tested package code: `cd132769a44a5bc0edbf7f28addff2b5a2fe8a47` (2026-10-15).

## GitHub verification

The [R CMD check workflow](https://github.com/ksublee/emhawkes/actions/runs/38681705953)
succeeded in all five environments: Windows release, macOS release, and Ubuntu
release, devel and oldrel-1. The [documentation site build](https://github.com/ksublee/emhawkes/actions/runs/38681705959)
also succeeded. Exact commit, job states and URLs are recorded in
[ci.json](auto-r-results/ci.json) and [site.json](auto-r-results/site.json).

## Complete fitting times

Both engines use the same installed development version. R selects
`emhawkes.engine = "r"`; auto selects `"auto"`.
Three fresh R processes per engine and model; elapsed-time medians in seconds.
Each timed call is the complete public hfit(), including preparation, maxLik
optimization with numerical differentiation, and the default maxLik Hessian.
Each pair uses identical data, starts, parameter sharing, constraints, initial
state arguments and optimizer settings. No analytic gradient is requested.

| Model | Events | R | Auto | R / auto |
| --- | ---: | ---: | ---: | ---: |
| Univariate Hawkes | 2000 | 0.98 | 0.07 | 14.00x |
| Univariate, two kernels | 2000 | 6.45 | 0.31 | 20.81x |
| Hawkes flocking | 1999 | 129.02 | 86.32 | 1.49x |
| Bid-ask, current components | 999 | 42.10 | 30.46 | 1.38x |

Machine: Intel Core Ultra 7 356H, Windows 11, R 4.6.0, maxLik 1.5-2.2 and
Rcpp 1.1.2. Worker processes were restricted to logical CPU 2 (affinity mask 4).
Each worker verifies the mask before its warm-up. This affects only benchmark
processes. Package loading, simulation, warm-up and the explicit gc() before
fitting are excluded. A 101-observation fit with iteration limit 2 warms each
worker. Engine order alternates by repetition; runs are sequential.

Elapsed and process CPU times, affinity and warnings are recorded in
[timings.csv](auto-r-results/timings.csv). Short native fits are sensitive to timer
resolution. These ratios describe this dataset and machine, not a universal
speed guarantee. The numerical-only models benefit most from native recursion;
flocking and Bid-ask still evaluate their callbacks in R.

An initial unrestricted series had substantial timing variability. Its completed
measurements were retained in [unrestricted-summary.csv](auto-r-results/unrestricted-summary.csv)
and [unrestricted-timings.csv](auto-r-results/unrestricted-timings.csv). That series
was interrupted to verify native dispatch and resumed for missing fits. The
primary table uses the subsequent uninterrupted series on the fixed CPU.

## Numerical agreement

All six stored fields were exactly identical between auto and R in every pair:
estimates, maximum likelihood, gradient, Hessian, convergence code and iterations.
Each engine also reproduced these fields exactly across all three repetitions.
All fits converged with code 0. Field comparisons are in
[metrics.csv](auto-r-results/metrics.csv), with coefficients in
[estimates.csv](auto-r-results/estimates.csv).

At the starting vector and both fitted vectors, direct C++ and R likelihoods
also agree; see [common-points.csv](auto-r-results/common-points.csv). These native
objective checks exercise the fitting adapter for callback models, where
standalone likelihood engine selection still uses R.

## Reproduction

Install the tested package into a private library. From the repository root,
use a fresh output directory and the same library for every command:

```
Rscript --vanilla tools/benchmark-auto-r.R prepare LIB OUT
Rscript --vanilla tools/benchmark-auto-r.R run LIB OUT 3
Rscript --vanilla tools/benchmark-auto-r.R compare LIB OUT
```

On Windows, set EMHAWKES_BENCH_AFFINITY=4 before run to pin each worker to logical
CPU 2. Without that variable the operating system selects CPUs. The script checks
that the requested mask was applied. Use a mask allowed on the target machine.
The run mode resumes completed fit files; do not reuse an output directory after
changing the installed package or fixture definitions.

The preparation uses fixed seeds for the numeric models and extracts the seeded
flocking and Bid-ask examples using the existing benchmark helpers. Callback
bodies are preserved. Raw fixtures and fit snapshots remain local. Full session
information is in [session-info.txt](auto-r-results/session-info.txt).
