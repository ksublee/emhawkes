# Numerical C++ fitting for the bid-ask model

> Historical benchmark: the explicit `"cpp"` setting below records the measured configuration. Current `hfit()` also selects C++ automatically for these supported models with the default `"auto"`; `"cpp"` requires support.

`hfit()` supports the bid-ask price model with `options(emhawkes.engine = "cpp")`
and the default numerical differentiation. No Jacobian is needed. Default
`"auto"` engine selection and the default maxLik Hessian are unchanged.

The model combines an observed-count baseline with intensity-dependent
post-event resets. R evaluates the baseline path once per proposal. C++
computes decay, compensator and pre-jump states, calls the bound R impact,
and applies its jump only to the selected kernel columns. A prepared named
call avoids an additional R wrapper at each event. R callbacks still account
for much of the computation, so the speedup is smaller than for numeric kernels.

## Callback contract

Observed `mu` functions may use `n`, `N`, `Nc`, `mark`, `type` and
`inter_arrival`. They return a plain numeric dimensions-by-one matrix,
constant on each interval. Index `n` supplies the baseline for the interval
ending at event `n`; use counts at `n - 1` for its predictable value.
Initialization uses `mu(n = 1)`. Baselines that use latent intensities or
change within an interval are unsupported by this adapter.

Impacts explicitly declare `lambda_component_n`, `lambda` or
`lambda_component` to select sequential evaluation. The current matrix and
history row `n` contain the state before that event's jumps. Component
histories retain the R engine's row-major order and column names; future
rows are zero. Counts at `n` already include the event. With an observed
baseline, the impact's `mu` argument is the bound function, as in R.

`lambda_component_n` avoids allocating full histories when only the current
state is needed. Full histories are allocated only when declared. The native
loop copies shared histories before writing subsequent rows; retained and
lazy callback arguments remain snapshots. Invalid matrix shapes/classes and
callback errors propagate. Nonfinite impact proposals return `-Inf`.

Callbacks must be deterministic, with captured inputs held fixed. Undeclared
latent access through `...` does not select sequential evaluation. Analytic
gradients do not cover this state feedback. Mark densities and residual
distributions remain unsupported. Standalone `logLik()` keeps its numeric-only
C++ contract; restore the R or automatic engine to evaluate callback models.

## Accuracy

The baseline is development commit `8d4f343`. For both callback interfaces,
coefficients, log-likelihood, gradient, Hessian, convergence code and iteration
count are exactly identical between the previous development package, the
current R engine and the C++ numerical fit. Every repeated fit also agrees.
The original-history and current-state variants have the same results.

CRAN 0.9.8 stops after 40 iterations and the development fits after 41.
The maximum coefficient difference is `1.10e-5`, and the achieved likelihood
difference is `9.26e-9`. Common-parameter likelihood differences are below
`1e-10`. These are numerical compatibility checks, not evidence of reliable
Wald standard errors or identification.

The vignette now uses `lambda_component_n`. With seed 1107 its 999 event
times, types and marks exactly reproduce the original vignette data.
Both callback interfaces are benchmarked on that same frozen data.

`tests/state-impact-numeric.R` covers an independent full-history likelihood
with feedback, signed components, unequal decay rates, zero intervals,
non-square overlapping mappings, count initial states, shared parameters,
public fitting and restoration, automatic initialization, snapshot retention,
lazy arguments, garbage collection, integer impacts and error propagation.
Existing numerical and analytic tests are retained. All 18 standalone tests
and all 18 tests under the source-package check pass; compiled-code, examples
and vignette checks pass too. Tests use maxLik 1.5-2.2 and the full source check
uses 1.6-10.

## Complete-fit timings

Windows 11, R 4.6.0, maxLik 1.5-2.2, 999 events plus the initial row.
Median elapsed seconds for three complete public `hfit()` calls, including
setup and the default final Hessian. Installation, simulation, loading and
explicit garbage collection are excluded. A short 100-event fit with a
two-iteration limit precedes timing. Libraries are isolated; runs are
sequential in a fresh process per package/interface. The user library is unchanged.

| Impact interface | CRAN 0.9.8 | C++ numeric | CRAN / C++ |
| --- | ---: | ---: | ---: |
| Original full history | 22.46 s | 17.36 s | 1.29x |
| Current matrix (`lambda_component_n`) | 19.16 s | 16.62 s | 1.15x |

The current matrix improves both R and C++ timings without changing the model.
Timing varies by machine. Raw measurements, estimates, fixed-point likelihoods
and session information are in `bid-ask-numeric-results/`. The current R
history timing uses three rounds; current-matrix R and previous-package
compatibility fits use one round. All CRAN and C++ medians use three rounds.
The optional `final_hessian = "stable"` shown in the vignette is excluded from
these timings so that every version uses its default final Hessian.

## Reproduction

Install CRAN 0.9.8, `8d4f343` and the working version in separate libraries,
with the same maxLik installation. From the repository root:

```text
Rscript --vanilla tools/benchmark-bid-ask-numeric.R prepare OLD_LIB OUT
Rscript --vanilla tools/benchmark-bid-ask-numeric.R worker OLD_LIB old bid_ask 1 OUT
Rscript --vanilla tools/benchmark-bid-ask-numeric.R worker CRAN_LIB cran bid_ask 3 OUT
Rscript --vanilla tools/benchmark-bid-ask-numeric.R worker NEW_LIB r bid_ask 3 OUT
Rscript --vanilla tools/benchmark-bid-ask-numeric.R worker NEW_LIB cpp bid_ask 3 OUT
```

Repeat the workers for `bid_ask_current`, then run `compare NEW_LIB OUT`.
For common-parameter checks run `points LIB ENGINE OUT` with `cpp`, `r`, and
`cran`, in that order. Raw RDS fixtures and fits remain local.

Source: [CRAN package page](https://cran.r-project.org/package=emhawkes).
