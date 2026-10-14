# Analytic gradients for observed-state impacts

`hfit()` now supports `gradient = "analytic"` for impacts that depend on observed counts, marks, event types, and times. Numeric differentiation remains the default. This extends the [parameter-only component adapter](function-gradient-results.md).

## Usage and scope

```r
count_impact <- function(param = c(a = 0.2), n, N, ...) {
  matrix(param["a"] * (1 + 0.2 * (N[n, 1] %% 2)))
}
count_jac <- function(param, n, N, ...) {
  out <- matrix(0, 1, length(param), dimnames = list(NULL, names(param)))
  out[, "a"] <- 1 + 0.2 * (N[n, 1] %% 2)
  out
}
h <- new("hspec", mu = 0.3,
  alpha = function(param = numeric(0)) matrix(0), beta = 1.2,
  impact = count_impact)
set.seed(1)
path <- hsim(h, size = 501, lambda_component0 = 0.1)
fit <- hfit(h, path$inter_arrival, lambda_component0 = 0.1,
  gradient = "analytic", jacobian = list(impact = count_jac))
```

The impact receives its local `param` vector as before. In contrast, `jacobian$impact` receives the **full named optimization vector**. Return a finite numeric matrix with one row per impact entry in column-major order and one named column per fitted parameter. Include zero columns for parameters with no effect. Column order may differ from the optimization vector.

The Jacobian must include the total impact derivative, including dependencies through the resolved `mu`, `alpha`, and `beta` matrices. Contributions for shared parameter names belong in the same column. Parameter-only core component Jacobians continue to use their local parameter vectors. A constant matrix can replace the Jacobian callback. An impact without local parameters still needs a Jacobian; it is zero only if the impact has no dependence on fitted parameters.

The package checks dimensions, names, and finite values. It does not derive or verify the user's derivative formulas. Callbacks must be deterministic and differentiable in parameters, return fixed-size plain numeric matrices, and keep captured inputs fixed. Branches based on observed counts are supported; a parameter-dependent branch requires a valid derivative at the point being evaluated.

Supported context is `n`, `N`, `Nc`, `mark`, `type`, `inter_arrival`, and resolved numeric `mu`, `alpha`, and `beta`. Declare needed `N` and `Nc` arguments explicitly in the impact, or supply those paths. Ignored `...` remains supported. Accessing latent `lambda`, `lambda_component`, or `lambda_component_n`, directly or through `...`, produces an error. Forcing all of `...` also forces these unavailable values. State-dependent baselines and latent-state feedback remain outside this option, as do mark densities, residual distributions, custom likelihoods, equality constraints, and BHHH.

The callback contract preserves the existing event convention:

- `N[n, ]` and `Nc[n, ]` include the current event. Row `n - 1` contains prior counts.
- The event's intensity and log contribution are evaluated before adding its impact.
- `type_col_map` selects the kernel columns updated by that event.
- The last event ends the likelihood window, so its jump has zero likelihood derivative.
- Automatic initialization retains the existing `mu`/`alpha`/`beta` rule. Additional impact does not enter it. Use an explicit `lambda_component0` when a different initial state is needed.

Default fitting and standalone `logLik()` retain their existing R path. Saved fits preserve the original callback, Jacobians, parameter map, and initial-state settings. A supplied `grad` still takes precedence.

## Calculation

R evaluates observed-state impacts into a components-by-observations matrix. Rcpp runs the intensity recursion and existing component scores while holding this impact matrix fixed. A reverse pass computes the score of each post-event impact. R then contracts those scores with `jacobian$impact` and adds the result to the existing mapped component gradient.

For component `(i, j)`, let `A_n` be the derivative of the likelihood with respect to its state just after event `n`. With `d_n = exp(-beta_ij * dt_n)` and `w_n = (1 - d_n) / beta_ij`, the backward recurrence is:

```text
A_last = 0
impact_score_ijn = A_n              if column j is selected by type[n]
impact_score_ijn = 0                otherwise
A_(n-1) = d_n * A_n - w_n + 1{i = type[n]} * d_n / lambda_type[n]
```

The integral weight uses `expm1` near zero. The reverse pass covers all kernel components, including unequal rates and overlapping mappings. Impact values and scores require O(observations * dimensions * kernels) storage. There is no observations-by-components-by-parameters tensor.

At a fixed parameter vector, the objective and gradient share the native result. Impact Jacobian callbacks run only when the optimizer requests a gradient, then their contracted result is cached. This avoids evaluating R Jacobians for line-search proposals that need only the likelihood.

## Accuracy and compatibility

`tests/observed-impact-gradient.R` checks the native likelihood against an independent full-history calculation and compares derivatives with central differences at two step sizes. It covers unequal rates, signed components and impacts, zero intervals, marks, nonsquare and overlapping mappings, explicit and automatic initial states, nonzero initial counts, shared parameter names, nonlinear core functions, and direct and indirect impact dependencies.

Additional checks cover callback errors, inactive and final-event scores, cache reuse, input and RNG preservation, no events, fixed optimization parameters, supplied gradients, and saved/restored fits. A one-parameter concave model has an independently calculated score root and information matrix; its estimated amplitude and stabilized variance agree with those references.

The benchmark uses three models at 500 and 2,000 events, each with three fitting modes. All 18 compared fits returned convergence code 0. Across the six analytic/default pairs:

- Maximum absolute log-likelihood difference: **2.68e-07**.
- Maximum absolute coefficient difference: **2.17e-05**.
- The original R likelihood at each analytic estimate matched its reported value in these runs.

Compiled likelihood with numeric derivatives produced the same estimates and maximum likelihoods as the R implementation in all six cases. The default path was also run with the previous installation at `6a0b869`. All six fits had identical estimates, maximum likelihoods, optimizer gradients, Hessians, return codes, and iteration counts. See the [baseline comparison](observed-impact-results/baseline-comparison.csv).

The full `R CMD check --no-manual --as-cran` with maxLik 1.6-10 finished with **Status: OK**, including all 14 test files, examples, and vignette rebuilding. The observed-impact, parameter-function, and numeric-gradient test suites also pass with maxLik 1.5-2.2. The usage example above was executed successfully. Updated website pages were rebuilt, and 111 local links and images were checked.

## Timings

Measured on Windows with R 4.6.0 and maxLik 1.6-10. Each data set has an initial row in addition to the stated events. All modes use the same starts, BFGS settings, positive parameter constraints, initial-state rules, and default final-Hessian option. Three alternating rounds follow warm-up; short fits are batched to exceed 0.2 seconds.

The table gives median elapsed seconds for optimization, including maxLik's final Hessian calculation. Objective setup, fit metadata, and the separate stabilized-Hessian comparison are excluded. The compiled/numeric mode is an internal benchmark adapter, not a new public fitting option. Public analytic fits were checked against the timed analytic objective, and default fits against the timed R objective.

| Model | Events | R + numeric (s) | Rcpp + numeric (s) | Rcpp + analytic (s) | Total speedup | Analytic gain within Rcpp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Univariate counts | 500 | 4.39 | 1.75 | 1.07 | 4.10x | 1.64x |
| Bivariate counts and marks | 500 | 6.67 | 3.11 | 1.14 | 5.85x | 2.73x |
| Multi-kernel counts | 500 | 1.71 | 0.67 | 0.55 | 3.11x | 1.22x |
| Univariate counts | 2,000 | 7.25 | 2.81 | 1.75 | 4.14x | 1.61x |
| Bivariate counts and marks | 2,000 | 13.50 | 6.61 | 2.81 | 4.80x | 2.35x |
| Multi-kernel counts | 2,000 | 6.41 | 2.65 | 1.86 | 3.45x | 1.42x |

Total improvement over the default R path was **3.11-5.85x**. Moving likelihood evaluation to Rcpp alone gave **2.04-2.58x**. Within the compiled path, adding analytic gradients gave another **1.22-2.73x**. These measurements isolate the two contributions more clearly than a single before/after ratio. R callback evaluation remains a material cost; the ratios depend on callback complexity, model, data, and optimizer behavior.

## Standard errors

Both public fitting modes were also compared with `final_hessian = "stable"` outside the timed runs. The 500-event univariate and bivariate cases were close to the imposed parameter boundary. Both modes reported `constraint_near_boundary` and retained unavailable standard errors. The other four cases had finite stabilized standard errors with maximum relative difference **6.24e-07**.

The new score does not resolve parameter identification or boundary inference. It preserves those diagnostics instead of manufacturing finite standard errors.

## Reproduction

Install the package into a separate library, select it through `R_LIBS`, and run from the repository root:

```text
Rscript tests/observed-impact-gradient.R
Rscript tools/benchmark-observed-impact.R tools/observed-impact-results 501,2001 3
```

For the previous-version comparison, install `6a0b869` into another library and select that library before running the same fixtures in snapshot mode:

```text
Rscript tools/benchmark-observed-impact.R /temporary/baseline 501,2001 3 all snapshot
```

Compare the saved `default-fits.rds` fields with the current snapshot. Timing runs also save this RDS file, but it is not tracked. The measured artifacts are [comparisons](observed-impact-results/comparisons.csv), [individual timings](observed-impact-results/timings.csv), [estimates and default standard errors](observed-impact-results/estimates.csv), and [session information](observed-impact-results/session-info.txt).
