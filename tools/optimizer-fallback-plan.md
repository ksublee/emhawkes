# Optional optimizer providers

`hfit()` now supports an explicit provider. The existing maxLik path is the
unchanged default; the stats path does not call the maxLik API and works when
maxLik is not installed.

```r
hfit(model, dt, optimizer = "maxLik")
hfit(model, dt, optimizer = "stats", method = "BFGS",
     control = list(maxit = 500, reltol = 1e-10))
```

## Implemented scope

- Reuse the existing R/C++ likelihood, analytic score, parameter sharing,
  callback binding and initial-state normalization. `emhawkes.engine` retains
  its existing meaning.
- Use `stats::optim()` for BFGS, CG, Nelder-Mead and L-BFGS-B. Bounds require
  L-BFGS-B, an admissible start and finite likelihood evaluations throughout.
- Use `stats::constrOptim()` for linear inequalities, with BFGS or Nelder-Mead.
  Map `A theta + B` to `ui = A_free` and
  `ci = -(B + A_fixed theta_fixed)`. Require a strictly interior start for
  movable rows; check fixed-only rows separately. Both constraint spellings
  select the same constraints.
- Preserve fixed parameters by logical mask, positions or names. Explicit
  derivatives take precedence over likelihood derivative attributes.
- Provide a separate `emhawkes_stats_fit` class with coefficient, log-likelihood,
  covariance, summary and print methods, plus existing reconstruction and
  diagnostics. Evaluation counts are retained as counts, not iterations.
- Default stats fits to stable final curvature. `final_hessian = "optimizer"`
  requests raw final curvature; `finalHessian = FALSE` disables it. Boundary
  estimates have unavailable Wald covariance even with raw final curvature.
- Make maxLik an optional dependency. Legacy maxLik fits still require maxLik
  and its methods; `fit_diagnostics()` and reconstruction load it on demand.
  Standalone use of a saved legacy object's maxLik methods requires loading
  maxLik. Saved stats fits need neither its namespace nor its methods.

## Deliberate limits

Equality constraints, Newton/BHHH/SANN methods, combined box and linear
constraints, and arbitrary maxLik controls are not supported by stats. Invalid
requests fail explicitly. Use stats controls such as `maxit`, not maxLik's
`iterlim`; vector-valued `parscale` and `ndeps` describe free parameters.
No automatic provider retry or dropping of constraints occurs.

Changing providers does not guarantee identical estimates: stopping rules and
finite differences can differ. This work supplies an independent execution path,
not evidence for improved finite-sample Wald coverage or faster fitting.

## Validation

Package tests cover exact Poisson estimates and information, observed-event
baselines, fixed/shared parameters, derivative precedence, constraint reduction,
boundaries, all-fixed fits, unsupported inputs, nonconvergence, result methods
and reconstruction. A separate R process finds only emhawkes, Rcpp and R's own
packages; it tests fitting, saved-result round trips and clear guidance when
maxLik is requested but unavailable.

`tools/stats-optimizer-validation.R` checks 48 multivariate fits across seeds,
kernel mappings, reduced/full parameter maps, numeric/analytic scores and R/C++
likelihoods. It includes mark effects, compares coefficients/likelihood/covariance
with maxLik, and independently replays fitted likelihoods with both engines.
All replay errors were zero. Poisoning maxLik API entry points did not affect
stats fitting, uncertainty or reconstruction. Sixteen old-versus-new maxLik
comparisons matched exactly, including numeric fit fields, classes and final RNG
state. Full package check results are recorded in `pre-release-audit.md`.

## References

- [R: optim](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html)
- [R: constrOptim](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/constrOptim.html)
