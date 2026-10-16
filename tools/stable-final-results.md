# Avoiding duplicate final Hessians in stable fits

## Change

For `final_hessian = "stable"`, unconstrained BFGS, CG, Nelder-Mead, and SANN
Hawkes fits omit maxLik's separate final Hessian unless `hess` or `finalHessian`
is supplied. These methods use maxLik's `maxOptim` path, where `finalHessian`
controls only the post-optimization calculation. Newton methods, constraints,
BHHH, custom likelihoods, and residual likelihoods retain their handling.
The default `final_hessian = "maxLik"` is unchanged.

Stable curvature still uses the same analytic-score or likelihood differences,
Richardson extrapolation, and positive-information and step-agreement checks.
A failed stable inversion still returns unavailable standard errors. It does not
use a different matrix merely because the optimizer's final pass was omitted.

`fit_diagnostics(fit)$computation$hessian$maxLik_skipped` records the omission.
`original_hessian` is `NULL` when that matrix was not calculated.
`computation$hessian$finalHessian` retains the caller's override, independently
of the internal skip. To retain the old two-pass stable calculation, supply
`finalHessian = TRUE`.

## Measurement

Measurements use Windows 11, R 4.6.0, and the five 2,500-event cases from the
[score-Hessian study](score-curvature-results.md). Before and after packages use
separate isolated libraries. The baseline installation's parsed runtime bodies
were verified against commit `2ac339b` before timing.

Each full public `hfit()` call runs in a fresh R process with affinity mask 4 and
a 100-event numerical warmup. The original starting values, BFGS controls,
numerical optimizer gradient, observations, initial components, mappings and
component Jacobians are unchanged. `auto` actually uses C++ in all five cases.
Versions run sequentially, with their order alternating across models. No other
benchmark, R test, or package build runs during timing.

Elapsed times include optimization, its returned score, final curvature, checks,
and fit metadata. They exclude package compilation, warmup, and serialization.
The main table uses one contemporaneous matched pair per model, rather than
comparing with historical timings. The first 4-by-8 case also has three additional
fresh pairs, with version order alternated across repeats.

After the main timings, the eligibility guard was tightened to honor abbreviated
`finalHessian` argument names accepted by maxLik. This leaves eligibility and
calculation paths unchanged for the timed arguments. The three additional pairs
and package checks use the final guard. Frozen and final source hashes are
recorded in the manifest.

## Results

All five fits retain the coefficients, likelihood, returned optimizer gradient,
convergence code, iteration count, final stable Hessian, covariance and summary
table **exactly**. All 96 standard errors are finite, positive, and bit-identical
to the baseline. Difference-step checks are also identical. Every fit returns
code 0, with no warnings. Each result also matches its saved result from the
previous score-Hessian study exactly.

| Model (dimensions x kernel columns) | Before (s) | After (s) | Speedup | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| 4 x 8, P=12: componentwise, observed baseline | 63.07 | 40.19 | 1.57x | 36.3% |
| 4 x 8, P=12: dense + alpha, observed baseline | 105.73 | 77.74 | 1.36x | 26.5% |
| 8 x 16, P=24: componentwise, constant baseline | 335.51 | 219.70 | 1.53x | 34.5% |
| 8 x 16, P=24: componentwise, observed baseline | 418.12 | 249.84 | 1.67x | 40.2% |
| 8 x 16, P=24: dense, observed baseline | 676.97 | 487.00 | 1.39x | 28.1% |

These are full-fit times. The observed reduction ranges from **26.5% to 40.2%**.
Performance depends on the fraction of fitting time spent on the omitted pass;
this change does not accelerate the optimizer iterations themselves.

The three additional final-source pairs also retain identical optimizer results,
Hessians and standard errors. Every pair is faster after the change. Paired
speedups range from **1.45x to 2.09x**, with a median of **1.68x**. Timings vary
across fresh processes; the paired values are retained below rather than treating
the first sample as a repeat-run median.

| Additional trial | Before (s) | After (s) | Speedup | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| 1 | 63.09 | 30.19 | 2.09x | 52.1% |
| 2 | 45.60 | 27.11 | 1.68x | 40.5% |
| 3 | 59.75 | 41.20 | 1.45x | 31.0% |

See [matched timings](stable-final-results/comparisons.csv),
[per-parameter comparisons](stable-final-results/coefficients.csv),
[repeat timings](stable-final-results/repeat-timings.csv),
[source manifest](stable-final-results/manifest.txt),
[hardware](stable-final-results/hardware.txt), and
[session information](stable-final-results/session-info.txt).

## Validation

Source build and `R CMD check --no-manual` pass with **Status: OK**, including
all **26 test scripts**, examples, and vignette rebuilding. There are no errors,
warnings, or notes. Updated Rd files also pass `tools::checkRd`.

Regression tests compare one-pass and explicit two-pass stable fits for BFGS, CG,
Nelder-Mead, SANN and their supported aliases, with both automatic and R engines.
They cover numerical, analytic and supplied optimizer gradients, fixed parameters,
likelihood-difference fallback, actual callback counts, unavailable curvature,
and saved diagnostics. Explicit full or abbreviated overrides, supplied Hessians,
Newton methods, constraints, custom objectives and residual likelihoods retain
their handling. Default and explicit maxLik paths retain the same evaluations.

All native sources and headers are byte-identical to the baseline. The parsed
final R bodies differ from the measured frozen bodies only by the conservative
argument-prefix guard; the timed inputs have no such override. This verification
is recorded in [source verification](stable-final-results/source-verification.txt).

The comparisons establish numerical agreement for these fits, not interval
coverage or reliable inference for unidentified or boundary parameters.

## Use

```r
fit <- hfit(h, inter_arrival, type, lambda_component0 = initial,
            jacobian = component_jacobians, final_hessian = "stable")
fit_diagnostics(fit)$computation$hessian$maxLik_skipped

# Retain the original maxLik final-Hessian pass as well.
reference <- hfit(h, inter_arrival, type, lambda_component0 = initial,
                  jacobian = component_jacobians,
                  final_hessian = "stable", finalHessian = TRUE)
```

For numeric components without callbacks, omit `jacobian`. The default Hessian
choice remains `"maxLik"`; this optimization applies only when stable curvature
is requested.
