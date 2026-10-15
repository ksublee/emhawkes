# Analytic gradients for parameter-only functions

`hfit()` now accepts exact component Jacobians for parameter-only `mu`, `alpha`, `beta`, and `eta` functions. Numerical differentiation remains the default.

```r
alpha_fun <- function(param = c(a = 0.2, b = 0.1)) {
  matrix(c(param["a"], 0, param["b"], param["a"]), 2)
}
alpha_jac <- cbind(a = c(1, 0, 0, 1), b = c(0, 0, 1, 0))
h <- new("hspec", mu = c(0.1, 0.1), alpha = alpha_fun,
         beta = matrix(0.9, 2, 2))
fit <- hfit(h, inter_arrival, type, gradient = "analytic",
            jacobian = list(alpha = alpha_jac))
```

For nonlinear transformations, supply a `function(param)` that returns the Jacobian at the current local parameter vector. For example:

```r
mu_fun <- function(param = c(log_mu = -1)) matrix(exp(param))
mu_jac <- function(param) {
  matrix(exp(param), dimnames = list(NULL, "log_mu"))
}
```

Rows follow `as.vector()` of the component matrix, in column-major order. Columns must have the local `param` names; their order can differ. A fixed entry has derivative zero. A component with `param = numeric(0)` requires no Jacobian. The saved fit metadata retains supplied Jacobians, and `fitted_hspec()` preserves the original component functions and their fitted parameter defaults.

Functions must have the single argument `param`, be deterministic and differentiable in it, and return plain numeric matrices of constant dimensions. Captured constants must remain fixed. Missing or malformed Jacobians produce an error. The package validates their structure and finite values, but does not derive or verify the derivative formulas. Supplying `grad` continues to override the built-in gradient.

## Calculation and scope

The existing Rcpp kernel returns component scores. For each function component `c(theta)`, the adapter applies `t(J_c) %*% score_c`. It then adds contributions at the indices from the saved parameter map. A name shared by several functions, or by a function and a numeric component, receives all contributions. Initial-state derivatives are included before this contraction. Fixed initial components remain constant.

Component functions and Jacobians run in R; likelihood and component scores run in Rcpp. One cached evaluation supplies both the objective and gradient at the same parameter vector. Default fitting and standalone `logLik()` retain their existing R path for function components. This change does not add automatic differentiation or compile R function bodies.

This report covers parameter-only components. Observed-state impacts are now supported separately; see the [observed-impact report](observed-impact-results.md). Observed baselines and current-component impacts are now supported with additional Jacobians; see the [state-gradient report](state-gradient-results.md). Full intensity histories, mark densities, residual-distribution models, custom likelihoods, equality constraints, and BHHH remain outside this option. A singular automatic initializer requires an explicit initial state. Finite matrices that overflow multivariate automatic initialization are rejected as invalid optimization proposals in the analytic path.

## Accuracy

`tests/function-gradient.R` compares gradients against central differences of the existing R likelihood at two step sizes and multiple parameter vectors, including directional derivatives. It covers nonlinear transformations, reordered Jacobian columns, fixed zeros, shared names across components, numeric/function mixtures, marks, overlapping multi-kernel mappings, explicit and automatic initial states, zero intervals, and no events. It also checks invalid inputs, cache reuse, a known transformed Poisson optimum and variance, fixed optimization parameters, optional stabilized curvature, supplied gradients, and saved/restored fits.

The full `R CMD check --no-manual --as-cran` finished with `Status: OK`, including all 13 test files, examples, and vignette rebuilding. The updated website pages were rebuilt and 111 local links and images were checked.

Both gradient test files pass with maxLik 1.5-2.2 and 1.6-10. The fit benchmark compares four models at 1,000 and 5,000 events. All 16 fits returned convergence code 0. Across the eight pairs:

- Maximum absolute log-likelihood difference: **5.00e-11**.
- Maximum absolute coefficient difference: **1.14e-06**.
- Evaluating the original R likelihood at each analytic estimate gave the same reported log-likelihood in these runs.

The default path was also run using the previous installation at `08ae157`. All eight fits had identical estimates, maximum log-likelihood, optimizer gradients, Hessians, return codes, and iteration counts. See the [baseline comparison](function-gradient-results/baseline-comparison.csv).

## Whole-fit timings

Measured on Windows with R 4.6.0 and maxLik 1.6-10. Each data set has one initial row in addition to its events. Both modes use the same data, starting values, BFGS settings, constraints, initial-state rules, and default final-Hessian option. Three alternating rounds are measured after warm-up; short fits are batched to exceed 0.2 seconds. The table reports median elapsed time per complete fit, including setup and final curvature.

| Model, 5,000 events | Default numeric (s) | Analytic (s) | Speedup |
| --- | ---: | ---: | ---: |
| Univariate, log parameters, automatic initial state | 3.9100 | 0.0425 | 92.00x |
| Bivariate, fixed zero excitation entry | 7.8600 | 0.1250 | 62.88x |
| Bivariate, shared name in baseline and decay, automatic initial state | 5.0800 | 0.1400 | 36.29x |
| Multi-kernel, nonlinear amplitudes, fixed decay | 3.8900 | 0.0550 | 70.73x |

The observed range across all eight pairs was **31.48-92.00x**. These are combined gains from moving function-model likelihood evaluation from R to Rcpp, using analytic derivatives, and reusing model setup and cached results. They do not isolate the gradient's contribution or guarantee the same speedup for other models and data.

The shared-name model uses common bounds `-8 <= theta <= 2` on its transformed parameters. Without these bounds, the existing numeric fit at 5,000 events proposed a scale parameter near 599 and overflowed automatic initialization. Both methods and the previous-version comparison use the same bounds. The other models retain their stated unconstrained or positive-baseline settings. These bounds do not impose general Hawkes stationarity.

## Standard errors

Default Hessians can differ because maxLik differentiates either likelihood values or the supplied score. One default numeric standard error was nonfinite; all analytic ones were finite in this benchmark. For a common comparison, both fits were also run with `final_hessian = "stable"`, outside the timing measurements. All resulting standard errors were finite and differed by at most **1.68e-05** in relative terms.

This establishes numerical consistency for these examples. It does not establish identifiability, a global optimum, or confidence-interval coverage.

## Reproduction

Install the package into a separate library, select it through `R_LIBS`, and run from the repository root:

```sh
Rscript tests/function-gradient.R
Rscript tests/native-gradient.R
Rscript tools/benchmark-function-gradient.R tools/function-gradient-results 1001,5001 3
```

Raw [timings](function-gradient-results/timings.csv), [fit comparisons](function-gradient-results/comparisons.csv), [coefficients and standard errors](function-gradient-results/estimates.csv), and [session information](function-gradient-results/session-info.txt) accompany this report. The script uses fixed seeds and checks convergence and agreement before recording each comparison. Keep CPU-heavy work separate from timing runs.

To reproduce the previous-version comparison, install `08ae157` into another library and run the same benchmark script with a separate output directory and the final arguments `1001,5001 3 all snapshot`. Compare its `default-fits.rds` with the new run's file; both store the same six fit fields.
