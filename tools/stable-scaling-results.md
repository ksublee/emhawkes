# Stable Hessians on the five scaling fixtures

This study applies the existing `final_hessian = "stable"` option to the five
Markov fixtures whose default numerical-gradient fits had indefinite observed
information. It adds benchmark tools and results; it does not change package
algorithms, defaults, or the public interface.

Runtime source: `4d6ea845cbce711b3b778bde2f0decfe07a75138`. The starting repository
commit is `59814894377abc12fe4931f583b5974ea615e320`, which adds the previous scaling
study without changing package runtime code.

Measurements: 2026-10-16, Windows 11, Intel Core Ultra 7 356H, about 31.4 GiB
installed RAM, R 4.6.0, maxLik 1.5-2.2 and Rcpp 1.1.2. Raw
[session information](stable-scaling-results/session-info.txt) and
[hardware](stable-scaling-results/hardware.txt) are included.

A [follow-up study](score-curvature-results.md) compares the new analytic-score
curvature path with these likelihood-difference results.

## Design

Every case uses the same model-generated observations, starting parameters,
initial components and BFGS controls as the [scaling study](gradient-scaling-results.md).
There are 2,500 events in each fit. The five cases are:

- Four dimensions, eight kernel columns, 12 parameters, componentwise state feedback.
- Four dimensions, eight kernel columns, 12 parameters, dense feedback coupled to alpha.
- Eight dimensions, 16 kernel columns, 24 parameters, componentwise feedback.
- Eight dimensions, 16 kernel columns, 24 parameters, dense feedback.
- Eight dimensions, 16 kernel columns, 24 parameters, componentwise feedback and a parameter-only baseline.

Other cases have an observed-count-dependent baseline. Componentwise feedback
means each impact entry reads its own current component; event types can still
excite other dimensions. All impacts use current Markov state. The coefficient
and standard-error tables refer to the fitted log-scale parameters.

Both modes use the default numerical gradient and `options(emhawkes.engine = "auto")`.
Each mode/case runs once in a fresh R process, with processor affinity mask 4.
No timing worker runs alongside another benchmark or derivative check. Both modes warm up the numerical fit on the same
short prefix before timing. Timings include public `hfit()` setup, optimization,
the maxLik final Hessian, and, for `stable`, its curvature checks and replacement.
These are single matched measurements, not repeated-run medians.

No package functions are traced in the final timed fits. Earlier instrumented
stable runs are excluded from the published timing samples. Resuming the study
reuses completed untraced default runs and reruns stable in fresh processes.
Independent score-based checks run after all timings and are excluded from fit
times. The stable option retains the maxLik Hessian calculation before performing
its additional checks; its total time includes both stages.

## Checks

The tool compares both full fits, verifies the default result against the saved
scaling result, and checks that stable preserves coefficients, likelihood,
gradient, convergence code and optimizer iteration count exactly. Its stored
original Hessian is compared with the default Hessian. Standard errors come from
`vcov(fit)` and are verified against `summary(fit)`.

For an independent curvature reference at the same fitted parameters, the tool
centrally differences the supplied analytic score with step
`2e-5 * max(1, abs(parameter))`. This is separate from stable's likelihood-based
scaled differences and Richardson extrapolation. The previous study checked
this reference at two steps, with a maximum scaled matrix change of 2.94e-8.

All five stable fits report status `ok` with positive definite observed
information. All **96 fitted-parameter standard errors are finite**.
The default fits have **42 NaN standard errors in total**, and all five have
indefinite observed information. A finite standard error alone does not establish
usable curvature.

Coefficients, likelihoods, stored gradients, convergence codes and optimizer
iteration counts are **exactly identical** between each default/stable pair.
The default results also match the prior saved scaling fits exactly. Stable's
stored original Hessian matches the default Hessian. All fits return code 0.

## Standard-error accuracy

Against central differences of the analytic score at the same fitted vectors,
the largest relative standard-error discrepancy is **0.0179%** across all
96 coefficients. This measures numerical agreement with an independent
curvature calculation; it does not establish confidence-interval coverage.

| Model (dimensions x columns) | Default NaN SE | Stable NaN SE | Stable maximum SE discrepancy |
| --- | ---: | ---: | ---: |
| 4 x 8, P=12: componentwise | 0 | 0 | 0.00475% |
| 4 x 8, P=12: dense + alpha | 0 | 0 | 0.00227% |
| 8 x 16, P=24: componentwise | 20 | 0 | 0.00933% |
| 8 x 16, P=24: dense | 18 | 0 | 0.0167% |
| 8 x 16, P=24: constant baseline | 4 | 0 | 0.0179% |

For example, the componentwise 4-by-8 model's `log_alpha2` standard error changes
from **1.4899** to **0.2743**; the score-based reference is **0.2743**.
The default value was finite despite its unstable Hessian.

See [per-parameter standard errors](stable-scaling-results/coefficients.csv),
[reference comparisons](stable-scaling-results/reference-checks.csv), and
[reference standard errors](stable-scaling-results/reference-standard-errors.csv).

## Complete estimation time

Every row uses 2,500 events. Times include setup, optimization and final curvature
calculations. Each row is one measured default/stable pair. Differences
include normal run-to-run timing variation and are not guaranteed overheads.

| Model (dimensions x columns) | Default (s) | Stable (s) | Additional (s) | Time ratio |
| --- | ---: | ---: | ---: | ---: |
| 4 x 8, P=12: componentwise | 24.32 | 59.13 | +34.81 | 2.43x |
| 4 x 8, P=12: dense + alpha | 42.73 | 128.53 | +85.80 | 3.01x |
| 8 x 16, P=24: componentwise | 350.10 | 712.47 | +362.37 | 2.04x |
| 8 x 16, P=24: dense | 236.81 | 959.66 | +722.85 | 4.05x |
| 8 x 16, P=24: constant baseline | 262.17 | 332.86 | +70.69 | 1.27x |

The stable calculation uses 3 finite-difference levels in
these fits. Each level evaluates the likelihood `2 * P^2` times: 288 evaluations
for 12 parameters and 1,152 for 24 parameters. An initial center evaluation is
also made. The extra likelihood evaluations explain the substantial cost,
especially for event-by-event R callbacks. The maxLik Hessian is computed before
these additional evaluations and is included in the reported stable times.

## Practical use

The existing stable option improves numerical standard errors in all five tested
models while retaining their fitted maxima. Its cost supports keeping it an
explicit option. These timings and results apply to this model family, data and
machine; they do not establish that every Hessian failure can be repaired.

Stable fitting requires no component Jacobians. The supplied analytic Jacobians
are used only for this study's independent reference check. To request the
existing option on an eligible Hawkes fit, set `final_hessian = "stable"` and
inspect `summary(fit)` and `fit_diagnostics(fit)`.

See [fit comparisons](stable-scaling-results/comparisons.csv),
[step checks](stable-scaling-results/steps.csv), and
[captured notices](stable-scaling-results/notices.csv).

## Reproduce

Use the same isolated package library and private fixture directory prepared by
`tools/benchmark-gradient-scaling.R`. Keep sample RDS files outside the repository.
`FIXTURES` is the completed prior scaling output, including `samples/` and
`fitted-reference-checks.csv`, rather than its `fixtures/` subdirectory.
The previous study's `run`, `compare` and `check-gradient-scaling-fits.R` steps
produce these inputs. From the repository root:

```r
# Set EMHAWKES_BENCH_AFFINITY=4 on Windows for the recorded affinity.
# Rscript tools/benchmark-stable-scaling.R run LIB FIXTURES OUT
# Rscript tools/benchmark-stable-scaling.R compare LIB FIXTURES OUT
# Rscript tools/benchmark-stable-scaling.R reference LIB FIXTURES OUT
```

`run` resumes completed samples. Use an empty output directory when changing
models or timing protocols. See the CSV files for per-parameter standard errors,
step checks, messages and independent-reference comparisons.
