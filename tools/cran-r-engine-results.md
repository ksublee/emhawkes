# CRAN and current default fitting times

CRAN emhawkes 0.9.8 is compared with the current source at `b8e5c4b`.
The [CRAN package page](https://cran.r-project.org/package=emhawkes) identifies
0.9.8 as the public release. This comparison is separate from the earlier
`63dd422` versus `b8e5c4b` R likelihood optimization benchmark.

Both packages use Windows 11, R 4.6.0 and maxLik 1.5-2.2. They are installed in
separate temporary libraries and loaded in fresh R processes. Data and starting
parameters are frozen. After one warm-up fit, three complete `hfit()` calls are
measured. The table gives median elapsed seconds. Simulation, package loading
and explicit garbage collection are outside the timed interval. The default
numerical gradient and maxLik final Hessian are used in both versions.

| Model | Events | CRAN 0.9.8 | Current default |
| --- | ---: | ---: | ---: |
| Univariate count example | 2,000 | 4.26 s | 4.77 s |
| Bivariate counts and marks | 2,000 | 9.15 s | 11.26 s |
| Univariate multi-kernel counts | 2,000 | Error | 4.72 s |
| Hawkes flocking vignette | 1,999 | 73.80 s | 77.08 s |
| Bid-ask price vignette | 999 | 25.86 s | 27.72 s |

The flocking and bid-ask examples use the vignette's `size = 2000` and
`size = 1000`. Their first row records the initial state. Bid-ask uses the
maxLik final Hessian in this comparison; the vignette's optional `stable`
recalculation is excluded because CRAN has no such option. For the bivariate
count example, the marked count path is supplied to both packages because
CRAN `hfit()` has no `Nc0` argument. Initial marked counts remain 3.

These measurements do not show a default fitting speed advantage over CRAN.
The earlier 24–28% reduction refers to the immediately preceding development
version. Absolute elapsed times can vary with machine load. Another R benchmark
was active during part of the flocking and bid-ask measurements, so their small
time differences should not be treated as precise performance regressions.

## Results

All comparable successful fits converge with code 0.

| Model | Maximum coefficient difference | Log-likelihood difference | CRAN / current iterations |
| --- | ---: | ---: | ---: |
| Univariate count example | 7.63e-8 | 5.28e-11 | 48 / 49 |
| Bivariate counts and marks | 1.71e-8 | 2.86e-11 | 46 / 46 |
| Hawkes flocking | 2.78e-8 | 7.46e-11 | 84 / 84 |
| Bid-ask price | 1.10e-5 | 9.25e-9 | 40 / 41 |

The univariate multi-kernel example stops in CRAN before optimization with
`non-conformable arrays`. Its intensity branch adds a 1-by-1 baseline matrix
to a 1-by-2 kernel state. Therefore no comparable CRAN fit time is available.
The current version evaluates a scalar total intensity and fits this case
with convergence code 0.

Raw repetitions, comparison fields, coefficients, source hashes and R session
information are in [cran-r-engine-results](cran-r-engine-results/).
The package implementation is unchanged by this benchmark.

## Reproduction

Run from the repository root with R, package dependencies and an R compiler
available. The script installs both packages into temporary libraries and
checks that each worker loads the intended library.

```sh
Rscript --vanilla tools/benchmark-cran-r-engine.R install
Rscript --vanilla tools/benchmark-cran-r-engine.R prepare
Rscript --vanilla tools/benchmark-cran-r-engine.R run
Rscript --vanilla tools/benchmark-cran-r-engine.R compare
```

To measure one case separately:

```sh
Rscript --vanilla tools/benchmark-cran-r-engine.R worker cran user_example
Rscript --vanilla tools/benchmark-cran-r-engine.R worker current user_example
```

The available case names are `user_example`, `bivariate_counts_marks`,
`multi_kernel_counts`, `flocking` and `bid_ask`.
