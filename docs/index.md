# emhawkes

`emhawkes` specifies, simulates, and fits exponential Hawkes models,
including marked, state-dependent, and multi-kernel models.

[`library`](https://rdrr.io/r/base/library.html)`(`[`emhawkes`](https://github.com/ksublee/emhawkes)`)`` `` ``h`` ``<-`` ``new``(``"hspec"``,`` `` mu ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.2``, ``0.2``)``,`` `` alpha ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0.5``, ``0.9``, ``0.9``, ``0.5``)``, ``2``, byrow ``=`` ``TRUE``)``,`` `` beta ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``2.25``, ``2``, ``2``)`` ``)`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1107``)`` ``path`` ``<-`` `[`hsim`](https://ksublee.github.io/emhawkes/reference/hsim.md)`(``h``, size ``=`` ``1000``)`` ``fit`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h``, ``path``$``inter_arrival``, ``path``$``type``)`` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit``)`

Use `fitted_hspec(fit)` to restore the estimated model and
`fit_initial_state(fit)` to retrieve the initial-state arguments:

`estimated`` ``<-`` `[`fitted_hspec`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)`(``fit``)`` ``initial`` ``<-`` `[`fit_initial_state`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)`(``fit``)`` ``inferred`` ``<-`` `[`do.call`](https://rdrr.io/r/base/do.call.html)`(``logLik``, `[`c`](https://rdrr.io/r/base/c.html)`(`[`list`](https://rdrr.io/r/base/list.html)`(``object ``=`` ``estimated``,`` `` inter_arrival ``=`` ``path``$``inter_arrival``, type ``=`` ``path``$``type``, mark ``=`` ``path``$``mark``,`` `` infer ``=`` ``TRUE``)``, ``initial``)``)`

The result remains a `maxLik` object, so
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) keep their existing
behavior. The saved parameter map preserves shared parameters and fixed
constants in component functions. Simulation generators are retained,
and their named parameters are updated when they match estimated
coefficient names. Automatic initial intensity remains `NULL` and is
recalculated for the restored model. Initial counts are retained. A
custom objective is not incorporated into the restored model.

Save the fit with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to
retain this information between sessions. Keep the observed data
separately and supply them again for inference. The stored optimization
objective can retain data in its function environment. The saved initial
state describes the start of the observation window. Pass it to
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) to
start a new path from the estimated model. Continuing after the last
observed event requires that event’s state. Fits created without the new
metadata need to be refitted before using the reconstruction helpers.

A successful optimizer return code does not establish reliable standard
errors. Check the score, curvature, numerical difference steps, and
multiple starting values. The [fit diagnostic
report](https://github.com/ksublee/emhawkes/blob/main/tools/fit-diagnostics.md)
reproduces the `NaN` standard errors in the extended examples and
separates Hessian cancellation from the nonsmooth `tzexp` likelihood.

By default, `reduced = TRUE` treats equal starting values **within each
parameter matrix** as one parameter. This example estimates four
parameters: one for `mu`, two for `alpha`, and one for `beta`.
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md)
explains this rule once before optimization.

Use `reduced = FALSE` to estimate the ten matrix entries separately,
even when their starting values are equal. Function components share
parameters by name under either setting. A zero in a numeric matrix is
also a starting value; use a component function to fix selected entries
at zero.

The first data row is the initial state. `inter_arrival`, `type`, and
`mark` start with zero, followed by the observed events. `hsim(size)`
includes this initial row. `N0` and `Nc0` specify initial counts and
mark totals; `Nc0` must be named in
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

The likelihood ends at the last supplied event. `N0` and `Nc0` accept a
scalar or one value per event type. Supplied `N` and `Nc` matrices must
agree with the event data and initial counts. Explicit `type_col_map`
requires initial kernel components for simulation and the built-in
likelihood.

`logLik(h, ..., infer = TRUE)` returns a list with the likelihood,
intensity paths, and residual vectors. Use it in place of the deprecated
[`infer_lambda()`](https://ksublee.github.io/emhawkes/reference/infer_lambda.md)
and
[`residual_process()`](https://ksublee.github.io/emhawkes/reference/residual_process.md)
helpers. Simulation returns an `hreal` object.

Numeric Hawkes specifications use a compiled Rcpp engine for likelihood
evaluation and inference, including multi-kernel models, unequal decay
rates, and numeric mark effects. Function-valued parameters,
user-defined impact or mark densities, residual-distribution models, and
numeric inputs with custom classes use the R implementation. Simulation
continues to use R.

[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) uses
this engine for its built-in likelihood. Gradient and Hessian
calculations remain numerical unless functions are supplied through
`grad` and `hess`. The Rcpp engine does not change the parameter-sharing
rules.

For comparison or troubleshooting, set `options(emhawkes.engine = "r")`
to use the R calculation, or `"auto"` to restore automatic selection.
The development setting `"cpp"` requires a supported numeric
specification and reports unsupported features. Installing from source
requires a C++ compiler and Rcpp; Windows source installations require
Rtools.

See the [basic
vignette](https://ksublee.github.io/emhawkes/articles/example.html) and
[extended
examples](https://ksublee.github.io/emhawkes/articles/extended_example.html).

The website is generated with `Rscript tools/build-docs.R` after
installing the current package and `pkgdown`. The script includes the
extended article from `vignettes_ext`, which is excluded from the CRAN
package. The [GitHub Markdown
example](https://github.com/ksublee/emhawkes/blob/main/github_document/example.md)
uses the same extended source.
