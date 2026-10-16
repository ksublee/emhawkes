# Stable Hessians from analytic scores

This change adds a score-difference path to `final_hessian = "stable"`.
It retains the default `gradient = "numeric"` and `final_hessian = "maxLik"`.
The five models are the same 2,500-event cases from the
[previous stable-Hessian study](stable-scaling-results.md).

## Implementation

With the `auto` or `cpp` engine, stable curvature attempts to use the built-in
analytic score. Numeric components need no user Jacobians. Function components,
observed baselines, and current-component impacts require the same Jacobians as
analytic optimization. These may be supplied with a numerical optimizer gradient.

The score is centrally differenced in scaled free-parameter coordinates. The
result is symmetrized and uses the existing Richardson extrapolation, positive
information checks, and Hessian/covariance/standard-error agreement across steps.
Fixed parameters do not move. The optimizer's coefficients, stored gradient,
likelihood, convergence code and iteration count are retained.

Current-state score evaluations obtain the likelihood and score in one native
replay when no value has already been requested. Line searches and optimization
retain their existing value-first behavior and deterministic replay checks.
Both value and score are cached for repeated requests at the same vector.

Missing support or Jacobians and failed score evaluations retain the likelihood-
difference path. A score-based matrix that fails the stability checks remains
unavailable; a second method is not searched just to obtain finite standard errors.
The explicit R engine retains likelihood differences. Supplied Hessians, BHHH,
disabled curvature, custom likelihoods, residual likelihoods, and boundary checks
retain their existing handling.

Inspect `fit_diagnostics(fit)$curvature_source`, `score_fallback_reason`, and
`computation$hessian$used`. The score path records `scaled_score_Richardson`;
the likelihood path records `scaled_central_Richardson`. The optimizer's gradient
record is independent of this final-curvature choice.

## Measurement

Measurements: 2026-10-16, Windows 11, R 4.6.0, the same machine and fixtures as the
previous study. The native core is unchanged from `4d6ea84`; the runtime changes
are in `R/hfit.R`, `R/hgradient.R`, and `R/huncertainty.R`.
[Manifest](score-curvature-results/manifest.txt),
[session information](score-curvature-results/session-info.txt), and
[hardware](score-curvature-results/hardware.txt) accompany the data.

Each model runs in a fresh R process, with affinity mask 4 and a 100-event
numerical warmup. The public fit uses the original starting vector, numerical
gradient, BFGS controls, observations, mapping, and initial components, with the
saved exact component Jacobians supplied for final curvature. No competing
benchmark, derivative check, or package build runs during timing.

After the public fit, both final-curvature paths are timed at exactly the same
fitted vector and original maxLik Hessian. Adapter construction and curvature
checks are included. The order alternates across models. These are one matched
pair per model, not repeated-run medians. Compilation, simulation, and result
serialization are excluded. No functions are traced.

An initial two-pass-score implementation was measured, then replaced with the
single-replay curvature score. Those initial samples are excluded; final timing
uses an empty output directory and a freshly installed runtime.

## Results

All five public fits retain the saved optimizer coefficients, likelihoods,
stored gradients, convergence codes and iteration counts exactly. Their original
maxLik Hessians are also exactly unchanged. All fits return code 0. The matched
likelihood-difference curvature reproduces the previous stable covariance exactly.

All **96 standard errors are finite and positive**. Against the existing stable
path, the maximum relative standard-error difference is **0.01790%**.
Against the previous small-step (`2e-5`) analytic-score reference, the maximum
relative difference is **1.08e-06%**. The latter checks numerical agreement
across difference schemes, not the supplied Jacobian formulas independently.
Existing R-reference gradient checks cover the score formulas; added tests check
one-replay scores against the existing scores and cached likelihoods.
The maximum relative covariance-matrix difference from existing stable is
0.01423% in Frobenius norm.

The following times measure **final-curvature processing only**, at the same
fitted vector. They include analytic-adapter construction and the same stability
checks. They exclude optimization and the original maxLik final Hessian.

| Model (dimensions x kernel columns) | Likelihood differences (s) | Score differences (s) | Speedup | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| 4 x 8, P=12: componentwise | 45.55 | 13.48 | 3.38x | 70.4% |
| 4 x 8, P=12: dense + alpha | 85.01 | 26.40 | 3.22x | 68.9% |
| 8 x 16, P=24: componentwise | 357.54 | 101.26 | 3.53x | 71.7% |
| 8 x 16, P=24: dense | 511.97 | 213.95 | 2.39x | 58.2% |
| 8 x 16, P=24: constant baseline | 271.34 | 89.84 | 3.02x | 66.9% |

The public full-fit times are recorded separately below. Previous full-fit times
in the raw CSV are historical measurements from the prior study, not a new
matched whole-fit benchmark. They should not be used to infer whole-fit speedup;
optimization timing varies between runs. The matched curvature table isolates
the effect of this change.

| Model | New public stable fit, including optimization (s) |
| --- | ---: |
| 4 x 8, P=12: componentwise | 62.56 |
| 4 x 8, P=12: dense + alpha | 57.30 |
| 8 x 16, P=24: componentwise | 460.09 |
| 8 x 16, P=24: dense | 687.11 |
| 8 x 16, P=24: constant baseline | 356.75 |

At each step level the score path needs `2 * P` score evaluations instead of
`2 * P^2` likelihood evaluations, plus a center score for diagnostics. A score
costs more than a likelihood evaluation, especially with dense state Jacobians;
the observed speedup is therefore smaller than the ratio of evaluation counts.

See [timing and agreement checks](score-curvature-results/comparisons.csv),
[per-parameter values](score-curvature-results/coefficients.csv), and
[difference-step checks](score-curvature-results/steps.csv).


The standard-error comparisons concern the fitted log-scale parameters. Numerical
agreement does not establish interval coverage or repair weak identification.

## Validation

Source build and `R CMD check --no-manual` pass with status `OK`, including
all 25 test scripts, examples, and vignette rebuilding. There are no errors,
warnings, or notes. Added tests cover independent Poisson and correlated-quadratic
curvature, free/fixed coordinates, reordered score names, one-replay score/value
cache agreement, missing Jacobians, runtime score failures, and step-sensitive
curvature. Default maxLik, explicit R, supplied Hessians, disabled curvature,
BHHH, custom likelihoods, and residual likelihoods retain their handling.

## Use

```r
# Numerical optimization; stable curvature can use the analytic score.
fit <- hfit(h, inter_arrival, type, lambda_component0 = initial,
            jacobian = component_jacobians, final_hessian = "stable")
fit_diagnostics(fit)$computation
fit_diagnostics(fit)$curvature_source
fit_diagnostics(fit)$score_fallback_reason
```

For numeric core matrices without callbacks, omit `jacobian`. To use the existing
R likelihood and likelihood-difference curvature explicitly, set
`options(emhawkes.engine = "r")` before fitting.

## Reproduce

Prepare the completed gradient-scaling fixtures and previous stable-scaling
samples using their tools. Keep private RDS fixtures and samples outside the
repository. Install the updated package in a separate `LIB`.

```text
Rscript --vanilla tools/benchmark-score-curvature.R run LIB FIXTURES PRIOR OUT
Rscript --vanilla tools/benchmark-score-curvature.R compare LIB FIXTURES PRIOR OUT
```

`FIXTURES` is the completed gradient-scaling output directory. `PRIOR` is the
previous stable-scaling output with full sample RDS files. On Windows, set
`EMHAWKES_BENCH_AFFINITY=4` for the recorded affinity. `run` resumes saved samples;
use a new `OUT` after changing the runtime or protocol. The public CSV files
contain coefficients, standard errors, step checks, and timings.
