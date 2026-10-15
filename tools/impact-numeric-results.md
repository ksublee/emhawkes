# Numerical C++ fitting for observed-state impacts

`hfit()` now accepts observed-state `impact` functions with
`options(emhawkes.engine = "cpp")` and the default numerical differentiation.
No Jacobian is needed. R evaluates the event impacts; C++ evaluates the
intensity recursion, exponential decay and compensator. Numerical gradient
and default final-Hessian calculations remain with maxLik.

The default `"auto"` selection is unchanged. Core component functions must
still depend only on `param`. Impacts can use `n`, `N`, `Nc`, `mark`, `type`,
`inter_arrival`, and resolved `mu`, `alpha`, and `beta`. They must be
deterministic and return plain numeric matrices of the kernel dimensions.
Declare needed count paths explicitly or supply them. Counts at row `n`
include that event. Its impact updates the selected kernel columns after
its intensity contribution. Automatic initialization retains the existing
core-parameter rule. A custom initial state remains available.

Reading latent `lambda`, `lambda_component`, or `lambda_component_n`, even
through `...`, produces an error. Forcing all of `...` also forces these
unavailable states. Ignored `...` and optional callback defaults are retained.
State-dependent baselines, mark densities and residual distributions remain
outside this adapter. Standalone `logLik()` retains its numeric-only C++
contract; restore the R or automatic engine before evaluating a fitted
callback model.

```r
previous_engine <- options(emhawkes.engine = "cpp")
fit <- hfit(h, inter_arrival, type, mark)  # numerical differentiation
options(previous_engine)
```

## Accuracy

The previous development baseline is `fd5f6ab`. The same observations,
starting parameters, parameter sharing, bounds, count paths, initial states
and optimizer controls are used for every comparison. No analytic gradients
or Jacobians are supplied.

For each model, compare coefficients, attained log-likelihood, stored
gradient, Hessian, convergence code and iteration count. The benchmark also
checks repeated fits and evaluates common parameter points across engines.
CRAN fits may stop at slightly different coefficients; fixed-parameter
likelihood agreement is checked separately from optimization termination.
All six result fields are exactly identical between the previous development
fit, the current R fit and the C++ fit in every model. CRAN coefficient
differences are at most `7.64e-8`; attained log-likelihood differences are at
most `7.46e-11`. All fits return convergence code zero. Across nine common
parameter points, R/C++ and CRAN/C++ log-likelihood differences are below
`1e-10`. These checks do not establish parameter identifiability or reliable
standard errors.

See `impact-numeric-results/accuracy.csv`, `estimates.csv` and
`fixed-points-*.csv` for the recorded differences.

The regression test `tests/observed-impact-numeric.R` checks public dispatch,
unchanged automatic selection, unused Jacobians, local/shared parameters,
count and mark timing, optional arguments, cache reuse, state restoration,
RNG preservation and explicit latent-state rejection. An independent
full-history likelihood covers signed components, unequal decay rates,
zero intervals and overlapping, reordered nonsquare kernel mappings.
The existing analytic-impact and parameter-only numerical tests are retained.

## Complete-fit timings

Windows 11, R 4.6.0, maxLik 1.5-2.2. CRAN 0.9.8 and development 0.9.9
are installed in private libraries and loaded in fresh processes for each
model/engine. The user's installed package is unchanged. The current CRAN
version was checked at [CRAN](https://cran.r-project.org/package=emhawkes).

Times cover complete public `hfit()` calls, including setup and the default
final Hessian. They exclude installation, simulation, package loading and
explicit garbage collection. A short fit on the first 100 events with a
two-iteration limit precedes measurement. Record three full fits per engine
and report their median; the previous-development compatibility check uses
one fit. Runs are sequential. Timing varies on this machine; raw per-round
measurements are retained in `impact-numeric-results/timings.csv`.

| Model | Events | CRAN (s) | Current R (s) | C++ numeric (s) | CRAN / C++ |
| --- | ---: | ---: | ---: | ---: | ---: |
| Univariate counts | 2000 | 9.34 | 9.91 | 4.92 | 1.90x |
| Bivariate counts and marks | 2000 | 9.19 | 7.84 | 4.69 | 1.96x |
| Hawkes flocking | 1999 | 119.31 | 100.48 | 51.58 | 2.31x |

The first model has numeric core matrices and a count-dependent impact.
The second has parameter-only core functions and an impact depending on
counts, marked counts and marks. The flocking fixture comes directly from
`vignettes_ext/extended_example.Rmd`; its observed count comparison is fixed
across candidate parameters. The univariate multi-kernel count example is
excluded from the CRAN timing table because CRAN 0.9.8 cannot evaluate it.
Kernel mappings remain covered by the independent regression test.

## Package checks

All 17 standalone regression scripts pass with the timing dependencies. A
fresh source build and `R CMD check --no-manual --as-cran` also finish with
`Status: OK`, including all 17 scripts, examples and vignette rebuilds.
The full check uses maxLik 1.6-10; timings use 1.5-2.2 consistently across
all package versions.

## Reproduction

Use separate library directories for CRAN, `fd5f6ab` and the working version,
with the same maxLik installation. Run from the repository root:

```text
Rscript --vanilla tools/benchmark-impact-numeric.R prepare OLD_LIB OUT
Rscript --vanilla tools/benchmark-impact-numeric.R worker CRAN_LIB cran user_example 3 OUT
Rscript --vanilla tools/benchmark-impact-numeric.R worker OLD_LIB old user_example 1 OUT
Rscript --vanilla tools/benchmark-impact-numeric.R worker NEW_LIB r user_example 3 OUT
Rscript --vanilla tools/benchmark-impact-numeric.R worker NEW_LIB cpp user_example 3 OUT
```

Repeat workers for `bivariate_counts_marks` and `flocking`, then run:

```text
Rscript --vanilla tools/benchmark-impact-numeric.R compare NEW_LIB OUT
Rscript --vanilla tools/benchmark-impact-numeric.R points NEW_LIB cpp OUT
Rscript --vanilla tools/benchmark-impact-numeric.R points NEW_LIB r OUT
Rscript --vanilla tools/benchmark-impact-numeric.R points CRAN_LIB cran OUT
```

The `prepare` command uses seeds 9351, 9352 and the vignette's 1107. It does
not estimate parameters. Raw RDS fixtures and fits stay local; portable CSVs
and session information accompany this report.
