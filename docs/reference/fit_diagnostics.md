# Inspect Final Curvature Diagnostics

Return the final Hessian checks requested through
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

A list describing the calculation status, reason, free parameters,
scale, and, when available, difference-step checks, curvature changes,
recalculated score, scaled condition number, and Newton decrement.

## Details

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
```
