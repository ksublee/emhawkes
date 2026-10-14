# Analytic gradient for numeric Hawkes models

`hfit(..., gradient = "analytic")` computes the likelihood and its score together in Rcpp. The default remains `gradient = "numeric"`. A supplied `grad` takes precedence; a supplied `hess` is retained.

```r
fit <- hfit(h, inter_arrival, type, mark,
            lambda_component0 = initial,
            gradient = "analytic")
```

Numeric `mu`, `alpha`, `beta`, and optional `eta` are supported, including shared parameters, `reduced = FALSE`, fixed parameters, unequal decay rates, signed components, and multi-kernel mappings. Explicit initial components are held fixed. Automatic initial components include their parameter derivatives. Singular automatic initialization requires explicit `lambda_component0`.

Function components, state-dependent callbacks, mark densities, residual-distribution models, custom likelihoods, equality constraints, and BHHH are excluded. Unsupported requests fail explicitly. The analytic option requires the `auto` or `cpp` engine.

## Calculation

For a component with post-event state `s`, decay `b`, and next interval `d`, the pre-event state is `s * exp(-b*d)`. Its integrated intensity is `s * W`, where `W = (1 - exp(-b*d))/b`.

The native loop propagates partial derivatives with respect to that component's `alpha`, `beta`, `eta`, and initial value. The score adds the derivative of each event's log intensity and subtracts the derivative of each interval integral. Small `b*d` uses a series for the beta derivative. Component-aware multiplication avoids losing finite results at extreme scales.

Each component evolves independently in this numeric model, so time is proportional to events times components. Score storage is proportional to components; no events-by-parameters array is stored. R combines component scores through the existing parameter map. A fit-local cache reuses the value and score at the same parameter vector.

For multivariate automatic initialization, let `K = alpha/beta` and `v = solve(I-K) %*% mu`. The existing initializer is `C[i,j] = K[i,j] * v[j]`. If `Q` is the native score with respect to `C`, the extra score follows from `w = solve(t(I-K), colSums(Q*K))` and `G[i,j] = (Q[i,j] + w[i]) * v[j]`. Add `w`, `G/beta`, and `-G*K/beta` to the `mu`, `alpha`, and `beta` scores. The univariate initializer retains its existing factor of one half. The initializer does not use `eta`, matching the existing likelihood.

## Accuracy

`tests/native-gradient.R` compares the native score with central differences of an independent full-history likelihood. It checks two difference steps and directional derivatives, shared and unshared parameters, overlapping/reordered mappings, marks, signed components, automatic initialization, zero intervals, no events, and extreme rates. It also checks a known Poisson optimum and information, fixed parameters, restored specifications, cache invalidation, supplied derivatives, and optional final-Hessian stabilization.

The tests passed with maxLik 1.5-2.2 and 1.6-10. The full package check exercised all 12 test files, examples, and vignette rebuilding. The existing default path was compared with the pre-change installation at `ff33fe7`: the checked likelihoods, fitted coefficients, optimizer scores, Hessians, and return codes were identical. See [baseline comparisons](gradient-results/baseline-comparison.csv).

The fit benchmark uses seven configurations at 2,000 and 10,000 events, with identical data, starts, constraints, and initial-state rules for both methods. All 28 compared fits returned convergence code 0. Across the 14 pairs, the largest log-likelihood difference was **5.47e-09**, and the largest absolute coefficient difference was **1.23e-05**. The latter occurred in the free-decay multi-kernel fit. These are numerical agreement checks, not a proof of a global optimum for every Hawkes model.

## Whole-fit timings

Measured on Windows with R 4.6.0 and maxLik 1.6-10. Both methods use the compiled likelihood and the default final Hessian. Timings include fit setup, optimization, and final curvature. Three alternating rounds are measured; short fits are batched to exceed 0.2 seconds per batch. The table reports medians. Each data set also has an initial row.

| Model, 10,000 events | Numeric (s) | Analytic (s) | Speedup |
| --- | ---: | ---: | ---: |
| Univariate, fixed initial state | 0.2800 | 0.0875 | 3.20x |
| Univariate, automatic initial state | 0.2800 | 0.1100 | 2.55x |
| Bivariate, shared parameters | 1.0000 | 0.2200 | 4.55x |
| Bivariate, automatic initial state | 1.0200 | 0.2300 | 4.43x |
| Bivariate, separate entries, fixed decay | 2.4200 | 0.3300 | 7.33x |
| Multi-kernel, fixed decay | 2.5500 | 0.5600 | 4.55x |
| Multi-kernel, estimated decay | 0.6600 | 0.1150 | 5.74x |

The observed speedup across all 14 pairs was **2.55-9.92x**. The gain includes fewer likelihood evaluations, one-pass value/score computation, caching, and reuse of the numeric model setup. These measurements are not a speed guarantee for other data or model sizes. Derivative-free methods such as Nelder-Mead do not use the score to choose optimization steps.

## Standard errors

A supplied analytic score changes how maxLik obtains its numerical final Hessian. Therefore, identical estimates do not imply identical default standard errors. Some default numerical fits in this benchmark had nonfinite standard errors; all analytic fits had finite ones.

To separate curvature error from optimizer differences, the benchmark also requests `final_hessian = "stable"` for both methods. Their standard errors differed by at most **3.57e-06** in relative terms. An additional reference differentiates the analytic score at three scaled steps. At the smallest step, the analytic default standard errors differed from this reference by at most **2.72e-07**. See [curvature comparisons](gradient-results/hessian-reference.csv).

These checks establish numerical consistency for the measured examples. They do not establish identifiability or confidence-interval coverage. Weak excitation, nearby decay rates, short samples, and boundaries still require statistical judgment.

## State-dependent extension

The next adapter should accept Jacobians for parameter-only functions and map their component scores to optimization parameters. For callbacks that depend only on observed `N`, `Nc`, or marks, hold those paths fixed when differentiating.

Callbacks that depend on `lambda_component` need a further step. If an event update is `s_next = F(s, theta)`, propagate `d s_next / d theta = F_s * d s / d theta + F_theta`. A reset such as `-l2 + xi * level/mid` must include the derivative of `l2`. The current independent-component score cannot represent such cross-component feedback by itself. Add state and parameter Jacobian contracts, then verify the resulting recursion against the existing R callback path. Parameter-dependent branch thresholds need separate differentiability checks.

## Constraints and reproduction

Unconstrained fitting and supported inequality-constrained optimization retain maxLik's interface. Equality constraints are excluded because maxLik 1.6-10's SUMT penalty gradient is inconsistent with its penalty when `eqB` is nonzero. For example, with `A=1`, `B=-1`, the penalty is `(theta-1)^2`, while the implemented gradient is `2*(theta+1)`. Numerical differentiation avoids that supplied-gradient path. Use the default numeric option for equality constraints until this dependency issue is resolved. BHHH requires event scores, whereas this implementation returns the total score.

Install the current source into a separate library, select it through `R_LIBS`, and run:

```sh
Rscript tests/native-gradient.R
Rscript tools/benchmark-gradient.R tools/gradient-results 2001,10001 3 all
```

Raw [timings](gradient-results/timings.csv), [fit comparisons](gradient-results/comparisons.csv), [coefficients and standard errors](gradient-results/estimates.csv), and [session information](gradient-results/session-info.txt) are saved alongside this report. Benchmark data use fixed seeds. Keep CPU-heavy work separate from timing runs.
