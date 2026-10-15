# Inspect Fit Computation and Curvature Diagnostics

Return the recorded calculation methods and final Hessian checks from
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

## Usage

``` r
fit_diagnostics(fit)
```

## Arguments

- fit:

  A result from
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

## Value

A list with recorded methods in `computation`, calculation status,
reason, free parameters, scale, and, when available, difference-step
checks, curvature changes, recalculated score, scaled condition number,
and Newton decrement.

## Details

`computation` records the engine, optimizer gradient source and stored
final Hessian source without re-evaluating the model. Its `engine` list
contains `requested`, `used` and `fallback_reason`. `used` lists engines
actually invoked in first-use order, including any partial C++ pass
before fallback. It can contain both `"cpp"` and `"r"`. A custom
likelihood uses `"custom"` and has no requested built-in engine.
Recording includes final curvature evaluations.

`computation$gradient$used` is `"numeric"`, `"analytic"`, `"user"` or
`"attribute"`. It describes the returned optimizer score, separately
from a recalculated score in stable curvature checks.
`computation$hessian$used` is `"maxLik_numeric"`,
`"scaled_central_Richardson"`, `"user"`, `"attribute"`, `"BHHH"`,
`"disabled"` or `"unavailable"`. The Hessian list also retains the
requested `final_hessian` choice and any `finalHessian` override. Its
`reason` explains a skipped stabilization or unavailable result. If both
a derivative function and likelihood derivative attributes were
supplied, the source is `NA` with a reason because precedence depends on
maxLik.

The saved record is a snapshot; later objective evaluations or option
changes do not update it. Older fits with model metadata but no
computation record return `recorded = FALSE` and `NA` for unrecorded
methods. Any previously saved gradient selection is retained. Their
engine is not inferred from current options.

`status = "ok"` means the numerical curvature checks passed, not that
the optimizer found a global maximum or that Wald intervals have correct
coverage. The numerical calculation checks information and
standard-error stability across steps. The stored optimizer gradient and
coefficients are not changed. A Newton decrement above 0.1 triggers a
score warning. Its square divided by two estimates the local quadratic
improvement in log-likelihood. An information condition number above
`1e8` triggers a weak-curvature warning; eigenvalues at or below `1e-10`
times the largest absolute eigenvalue do not support numerical
inversion. These thresholds are numerical diagnostics. With the default
`final_hessian = "maxLik"`, these additional checks are not performed
and the status is `"not_checked"`. Request `final_hessian = "stable"` to
run them. Earlier fits without these checks also have status
`"not_checked"`.

## See also

[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md),
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)

## Examples

``` r
h <- new("hspec", mu = .3, alpha = 0, beta = 1)
fit <- hfit(h, c(0, .2, .3, .5, 1), lambda_component0 = 0,
            fixed = c(FALSE, TRUE, TRUE), final_hessian = "stable")
fit_diagnostics(fit)$status
#> [1] "ok"
fit_diagnostics(fit)$computation
#> $recorded
#> [1] TRUE
#> 
#> $engine
#> $engine$requested
#> [1] "auto"
#> 
#> $engine$used
#> [1] "cpp"
#> 
#> $engine$fallback_reason
#> NULL
#> 
#> 
#> $gradient
#> $gradient$requested
#> [1] "numeric"
#> 
#> $gradient$used
#> [1] "numeric"
#> 
#> $gradient$reason
#> NULL
#> 
#> 
#> $hessian
#> $hessian$requested
#> [1] "stable"
#> 
#> $hessian$finalHessian
#> NULL
#> 
#> $hessian$used
#> [1] "scaled_central_Richardson"
#> 
#> $hessian$reason
#> NULL
#> 
#> 
```
