# Computation

This guide covers calculation choices for the built-in likelihood. Start
with [Get
started](https://ksublee.github.io/emhawkes/articles/example.html). For
model interpretation and standard errors, see [Fitting and
inference](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html).

## Working example

[`library`](https://rdrr.io/r/base/library.html)`(`[`emhawkes`](https://github.com/ksublee/emhawkes)`)`` `` ``h`` ``<-`` ``new``(``"hspec"``,`` `` mu ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.2``, ``0.2``)``,`` `` alpha ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0.5``, ``0.9``, ``0.9``, ``0.5``)``, ``2``, byrow ``=`` ``TRUE``)``,`` `` beta ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``2.25``, ``2``, ``2``)`` ``)`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1107``)`` ``path`` ``<-`` `[`hsim`](https://ksublee.github.io/emhawkes/reference/hsim.md)`(``h``, size ``=`` ``1000``)`` ``fit`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h``, ``path``$``inter_arrival``, ``path``$``type``)`` ``#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).`` ``#> Use reduced = FALSE to estimate matrix entries separately.`` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit``)`` ``#> --------------------------------------------`` ``#> Maximum Likelihood estimation`` ``#> BFGS maximization, 35 iterations`` ``#> Return code 0: successful convergence `` ``#> Log-Likelihood: -970.1408 `` ``#> 4 free parameters`` ``#> Estimates:`` ``#> Estimate Std. error t value Pr(> t) `` ``#> mu1 0.19095 0.01635 11.68 < 2e-16 ***`` ``#> alpha1.1 0.48217 0.07373 6.54 6.15e-11 ***`` ``#> alpha1.2 0.98625 0.09484 10.40 < 2e-16 ***`` ``#> beta1.1 2.07987 0.16728 12.43 < 2e-16 ***`` ``#> ---`` ``#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`` ``#> --------------------------------------------`

## Likelihood engines

With automatic engine selection, numeric Hawkes specifications use a
compiled Rcpp engine for likelihood evaluation and inference, including
multi-kernel models, unequal decay rates, and numeric mark effects.
Standalone likelihood and inference also support parameter-only
component functions returning plain numeric matrices. Observed baselines
and observed or sequential state-dependent impacts also use C++.
Baselines cannot depend on latent intensities. Unsupported function
signatures use R. Mark densities, residual-distribution models, and
numeric inputs with custom classes use R. Simulation continues to use R.
The fitting adapters below also select C++ automatically for supported
callback models.

The `emhawkes.engine` option controls the built-in likelihood:

| Setting | Role |
|----|----|
| `"auto"` (default) | Use C++ for supported models and R otherwise. Supported fitting adapters receive this acceleration automatically. |
| `"r"` | Force the R likelihood for comparison or troubleshooting. |
| `"cpp"` | Require C++ support; report an error for unsupported models. |

Setting `"cpp"` is unnecessary to accelerate models already supported by
`"auto"`. Inspect `fit_diagnostics(fit)$computation$engine` for the
engines actually used and any automatic fallback reason. Engine
selection does not enable analytic gradients: request
`gradient = "analytic"` explicitly; unsupported analytic gradients
produce an error. A custom `mylogLik` bypasses built-in engine
selection. Installing from source requires a C++ compiler and Rcpp;
Windows source installations require Rtools.

`hspec2`` ``<-`` ``h`` ``inter_arrival2`` ``<-`` ``path``$``inter_arrival`` ``type2`` ``<-`` ``path``$``type`

The following example compares standalone numeric likelihood
calculations.

This function restores the previous option after each call.

`engine_loglik`` ``<-`` ``function``(``engine``)`` ``{`` `` ``previous`` ``<-`` `[`options`](https://rdrr.io/r/base/options.html)`(``emhawkes.engine ``=`` ``engine``)`` `` `[`on.exit`](https://rdrr.io/r/base/on.exit.html)`(`[`options`](https://rdrr.io/r/base/options.html)`(``previous``)``)`` `` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(`[`logLik`](https://ksublee.github.io/emhawkes/reference/logLik.md)`(``hspec2``, inter_arrival ``=`` ``inter_arrival2``, type ``=`` ``type2``)``)`` ``}`` `[`c`](https://rdrr.io/r/base/c.html)`(``R ``=`` ``engine_loglik``(``"r"``)``, automatic ``=`` ``engine_loglik``(``"auto"``)``)`` ``#> R automatic `` ``#> -974.2809 -974.2809`

## Recorded calculation methods

The fit also records its calculation methods:

`calculation`` ``<-`` `[`fit_diagnostics`](https://ksublee.github.io/emhawkes/reference/fit_diagnostics.md)`(``fit``)``$``computation`` ``calculation``$``engine`` ``# Requested engine, actual engines and auto fallback reason.`` ``#> $requested`` ``#> [1] "auto"`` ``#> `` ``#> $used`` ``#> [1] "cpp"`` ``#> `` ``#> $fallback_reason`` ``#> NULL`` ``calculation``$``gradient`` ``# Numeric, analytic or user-supplied score.`` ``#> $requested`` ``#> [1] "numeric"`` ``#> `` ``#> $used`` ``#> [1] "numeric"`` ``#> `` ``#> $reason`` ``#> NULL`` ``calculation``$``hessian`` ``# Requested and applied final-Hessian methods.`` ``#> $requested`` ``#> [1] "maxLik"`` ``#> `` ``#> $finalHessian`` ``#> NULL`` ``#> `` ``#> $used`` ``#> [1] "maxLik_numeric"`` ``#> `` ``#> $reason`` ``#> NULL`

An auto fit can record both `"cpp"` and `"r"` if it switches during
estimation. The record includes final curvature evaluations and survives
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Older fits without
this record return `recorded = FALSE` and `NA` for unrecorded methods.
The existing curvature checks remain available at the top level.

## Analytic gradients

[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md) uses
this engine for its built-in likelihood. Numerical differentiation
remains the default. For numeric Hawkes models, request the Rcpp
analytic gradient:

`fit`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h``, ``path``$``inter_arrival``, ``path``$``type``, gradient ``=`` ``"analytic"``)`` ``#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).`` ``#> Use reduced = FALSE to estimate matrix entries separately.`

This option supports shared and fixed parameters, multi-kernel mappings,
numeric mark effects, and derivatives of automatic initial states. It
computes the likelihood and score with native recursions. Observed
baselines and current-state impacts require the Jacobians described
below. Full intensity histories, custom likelihoods, mark densities,
residual distributions, equality constraints, and BHHH are not yet
supported. Singular automatic initialization requires explicit
`lambda_component0`. A supplied `grad` takes precedence. A supplied
`hess` is retained; otherwise `maxLik` can differentiate the analytic
gradient for curvature. The optimizer’s path and last digits may differ
from numerical differentiation. The [gradient
report](https://github.com/ksublee/emhawkes/blob/main/tools/gradient-results.md)
describes accuracy, full-fit timings, and the state-dependent extension.

## Component functions and Jacobians

Parameter-only `mu`, `alpha`, `beta`, and `eta` functions can also use
the analytic option when supplied with their exact Jacobians. For
example, this function fixes one excitation entry at zero:

`alpha_fun`` ``<-`` ``function``(``param`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``a ``=`` ``0.2``, b ``=`` ``0.1``)``)`` ``{`` `` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``param``[``"a"``]``, ``0``, ``param``[``"b"``]``, ``param``[``"a"``]``)``, ``2``)`` ``}`` ``alpha_jac`` ``<-`` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``a ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``0``, ``0``, ``1``)``, b ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``, ``1``, ``0``)``)`` ``h_fun`` ``<-`` ``new``(``"hspec"``, mu ``=`` ``h``@``mu``, alpha ``=`` ``alpha_fun``, beta ``=`` ``h``@``beta``)`` ``fit_fun`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h_fun``, ``path``$``inter_arrival``, ``path``$``type``,`` `` gradient ``=`` ``"analytic"``, jacobian ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``alpha ``=`` ``alpha_jac``)``)`` ``#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, beta).`` ``#> Use reduced = FALSE to estimate matrix entries separately.`

Jacobian rows follow [`as.vector()`](https://rdrr.io/r/base/vector.html)
of the returned matrix; columns are named with its `param` names. For
nonlinear transformations, supply a `function(param)` that returns the
Jacobian at that parameter vector. Components and Jacobians must be
deterministic functions of `param` with fixed matrix dimensions. The
package does not derive or verify the formulas automatically. A
component with `param = numeric(0)` needs no Jacobian. Shared parameter
names receive the sum of all component contributions. Missing or
malformed Jacobians produce an error. The analytic fit evaluates
component functions in R and the likelihood and scores in Rcpp. With the
default `"auto"` engine, supported function-valued models also use C++
with numerical differentiation, without Jacobians:

`fit_numeric`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h_fun``, ``path``$``inter_arrival``, ``path``$``type``)`` ``#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, beta).`` ``#> Use reduced = FALSE to estimate matrix entries separately.`

The functions must return plain numeric matrices of constant dimensions
and be deterministic in `param`. Observed-state `impact` functions also
use this numerical path without Jacobians, including count- and
mark-dependent impacts and the Hawkes flocking example. Numerical C++
fitting also supports observed baselines and sequential impacts that
explicitly request current or historical intensity states, including the
bid-ask example. Mark densities and residual distributions remain
unsupported. Supplied `grad` and `hess` are retained. Automatic
selection uses R for unsupported specifications; forcing `"cpp"` reports
unsupported features. Callback errors propagate. Callbacks must be
deterministic with fixed captures; detecting an unsupported state read
through `...` can require a retry in R. Standalone
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md) and
`infer = TRUE` also use C++ for parameter-only component functions,
observed baselines, and observed or sequential state-dependent impacts.
The default automatic selection therefore accelerates replay of restored
models with these components. Each callback keeps its own defaults;
standalone replay does not apply parameter-sharing rules from fitting.

## Observed impacts and baselines

Observed-state `impact` functions can use the analytic option with
`jacobian = list(impact = impact_jac)`. They may depend on observed `N`,
`Nc`, marks, types, and times, and on the resolved `mu`, `alpha`, and
`beta` matrices. Counts at row `n` include that event. The event’s
impact affects subsequent intensities in the kernel columns selected by
its type.

For `jacobian$impact`, `param` is the **full named optimization
vector**. Return a matrix with one row per impact entry in column-major
order and one named column per fitted parameter, including zero columns.
Include derivatives through any core matrices used by the impact. The
impact function keeps receiving its local parameters. See
[`?hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md) for an
example and the callback contract. Current-component feedback can use
the analytic option with additional derivatives below. Default numerical
fitting uses C++ for supported callbacks. Declare `lambda_component_n`
to receive the current pre-jump component matrix. Declare `lambda` or
`lambda_component` only when full pre-jump histories are needed.
Observed baselines may use counts and marks, return a dimensions-by-one
matrix, and remain constant between events. They cannot depend on latent
intensities. See the [bid-ask
benchmark](https://github.com/ksublee/emhawkes/blob/main/tools/bid-ask-numeric-results.md)
for the state convention, compatibility checks and timings.

For an observed baseline, analytic fitting requires `jacobian$mu` with
one row per event dimension and one named column per **full fitted
parameter**, including zero columns. Its callback receives observations
and the full named `param`. Include dependencies through resolved core
components. The baseline itself continues to receive its local
parameters.

## Current-state Jacobians

An impact that declares `lambda_component_n` can use:

`fit`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h``, ``inter_arrival``, ``type``, ``mark``, lambda_component0 ``=`` ``initial``,`` `` gradient ``=`` ``"analytic"``, jacobian ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``mu ``=`` ``mu_jac``, alpha ``=`` ``alpha_jac``,`` `` impact ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``param ``=`` ``impact_param_jac``, state ``=`` ``impact_state_jac``)``)``)`

Here `impact$param` differentiates the impact with respect to all fitted
parameters **holding the current state fixed**. It includes direct
dependencies through core components. `impact$state` differentiates
impact entries with respect to current component entries. Both state
axes follow [`as.vector()`](https://rdrr.io/r/base/vector.html)
column-major order; the matrix is components-by-components. The package
propagates derivatives through decay, selected jumps, and automatic
initialization. Arbitrary reads of full `lambda` or `lambda_component`
histories remain unsupported for analytic fitting. All callbacks must be
deterministic, and supplied derivatives must exist at evaluated
parameters. The [Bid-ask analytic-gradient
section](#bid-ask-analytic-gradient) contains the complete Jacobians and
fit. See the [state-gradient
report](https://github.com/ksublee/emhawkes/blob/main/tools/state-gradient-results.md)
for accuracy checks and timings. Numerical differentiation remains the
default.

## Time and memory

When only the current state is needed, use `lambda_component_n`. Full
histories need storage proportional to events times dimensions times
kernels. If a callback retains these histories, snapshot-preserving
copies can make cumulative copying grow quadratically with event count.
Peak live memory and cumulative allocation are different costs. See the
[scaling
measurements](https://github.com/ksublee/emhawkes/blob/main/tools/scaling-results.md)
for time, allocation and memory comparisons.

See the [count, mark and flocking
benchmark](https://github.com/ksublee/emhawkes/blob/main/tools/impact-numeric-results.md)
for complete-fit timings and accuracy checks.

## Bid-ask analytic gradient

First run the [Bid-ask price
model](https://ksublee.github.io/emhawkes/articles/extended_example.html#bid-ask-price-model).
The code below continues in that R session, using `h_ba`, `hr_ba`,
`mle_ba`, `initial_ask_price`, `initial_bid_price`, `initial_level`, and
`initial_mid_price`.

The current-component version also supports an optional analytic
gradient. The baseline Jacobian uses all fitted parameter names; the
parameter-only alpha Jacobian uses its local names. For impact, supply
separate derivatives with respect to parameters (holding the current
state fixed) and the current state. All matrix entries use
[`as.vector()`](https://rdrr.io/r/base/vector.html) column-major order.
In particular, the state Jacobian differentiates the pre-jump row sums
subtracted by impact. Counts and mark totals are observed data and have
zero parameter derivative. The baseline uses counts before the event,
while impact uses counts after it.

`# Exact derivatives for the current-component Bid-ask example.`` ``bid_ask_jacobians`` ``<-`` ``function``(``initial_ask_price``, ``initial_bid_price``, ``initial_level``,`` `` ``initial_mid_price`` ``=`` ``(``initial_ask_price`` ``+`` ``initial_bid_price``)`` ``/`` ``2``)`` ``{`` `` `[`force`](https://rdrr.io/r/base/force.html)`(``initial_ask_price``)``; `[`force`](https://rdrr.io/r/base/force.html)`(``initial_bid_price``)``; `[`force`](https://rdrr.io/r/base/force.html)`(``initial_level``)``; `[`force`](https://rdrr.io/r/base/force.html)`(``initial_mid_price``)`` `` ``market_state`` ``<-`` ``function``(``n``, ``Nc``, ``after`` ``=`` ``FALSE``)`` ``{`` `` ``if`` ``(``n`` ``==`` ``1L``)`` `[`return`](https://rdrr.io/r/base/function.html)`(`[`list`](https://rdrr.io/r/base/list.html)`(``level ``=`` ``initial_level``, mid ``=`` ``initial_mid_price``)``)`` `` ``k`` ``<-`` ``if`` ``(``after``)`` ``n`` ``else`` ``n`` ``-`` ``1L`` `` ``ask`` ``<-`` ``initial_ask_price`` ``+`` ``Nc``[``k``, ``1``]`` ``-`` ``Nc``[``k``, ``2``]`` `` ``bid`` ``<-`` ``initial_bid_price`` ``+`` ``Nc``[``k``, ``3``]`` ``-`` ``Nc``[``k``, ``4``]`` `` `[`list`](https://rdrr.io/r/base/list.html)`(``level ``=`` ``initial_level`` ``+`` ``Nc``[``k``, ``1``]`` ``-`` ``Nc``[``k``, ``2``]`` ``-`` ``(``Nc``[``k``, ``3``]`` ``-`` ``Nc``[``k``, ``4``]``)``,`` `` mid ``=`` ``(``ask`` ``+`` ``bid``)`` ``/`` ``2``)`` `` ``}`` `` ``mu_jac`` ``<-`` ``function``(``param``, ``n``, ``Nc``, ``...``)`` ``{`` `` ``s`` ``<-`` ``market_state``(``n``, ``Nc``)`` `` ``J`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``4``, `[`length`](https://rdrr.io/r/base/length.html)`(``param``)``, dimnames ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``NULL``, `[`names`](https://rdrr.io/r/base/names.html)`(``param``)``)``)`` `` ``J``[`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``4``)``, ``'mu1'``]`` ``<-`` ``1`` `` ``if`` ``(``s``$``level`` ``>`` ``0``)`` ``J``[`[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``3``)``, ``'zeta1'``]`` ``<-`` ``s``$``level`` ``/`` ``s``$``mid`` `` ``J`` `` ``}`` `` ``alpha_jac`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``16``, ``5``, dimnames ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``NULL``,`` `` `[`c`](https://rdrr.io/r/base/c.html)`(``'alpha_s1'``, ``'alpha_m'``, ``'alpha_s2'``, ``'alpha_w1'``, ``'alpha_w2'``)``)``)`` `` ``alpha_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``16``)``, ``'alpha_s1'``]`` ``<-`` ``1`` `` ``alpha_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``5``, ``12``)``, ``'alpha_m'``]`` ``<-`` ``1`` `` ``alpha_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``9``, ``8``)``, ``'alpha_s2'``]`` ``<-`` ``1`` `` ``alpha_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``15``)``, ``'alpha_w1'``]`` ``<-`` ``1`` `` ``alpha_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``3``, ``14``)``, ``'alpha_w2'``]`` ``<-`` ``1`` `` ``impact_param_jac`` ``<-`` ``function``(``param``, ``n``, ``Nc``, ``...``)`` ``{`` `` ``s`` ``<-`` ``market_state``(``n``, ``Nc``, after ``=`` ``TRUE``)`` `` ``J`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``16``, `[`length`](https://rdrr.io/r/base/length.html)`(``param``)``, dimnames ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``NULL``, `[`names`](https://rdrr.io/r/base/names.html)`(``param``)``)``)`` `` ``J``[`[`c`](https://rdrr.io/r/base/c.html)`(``6``, ``7``, ``10``, ``11``)``, ``'xi'``]`` ``<-`` ``s``$``level`` ``/`` ``s``$``mid`` `` ``J`` `` ``}`` `` ``impact_state_jac`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``16``, ``16``)`` `` ``impact_state_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``6``, ``10``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``6``, ``10``, ``14``)``]`` ``<-`` ``-``1`` `` ``impact_state_jac``[`[`c`](https://rdrr.io/r/base/c.html)`(``7``, ``11``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``3``, ``7``, ``11``, ``15``)``]`` ``<-`` ``-``1`` `` `[`list`](https://rdrr.io/r/base/list.html)`(``mu ``=`` ``mu_jac``, alpha ``=`` ``alpha_jac``,`` `` impact ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``param ``=`` ``impact_param_jac``, state ``=`` ``impact_state_jac``)``)`` ``}`` `` ``ba_jac`` ``<-`` ``bid_ask_jacobians``(``initial_ask_price``, ``initial_bid_price``, ``initial_level``,`` `` ``initial_mid_price``)`

`mle_ba_analytic`` ``<-`` `[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md)`(``h_ba``, inter_arrival ``=`` ``hr_ba``$``inter_arrival``,`` `` type ``=`` ``hr_ba``$``type``, mark ``=`` ``hr_ba``$``mark``,`` `` lambda_component0 ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``1``, ``4``, ``4``)``, gradient ``=`` ``"analytic"``,`` `` jacobian ``=`` ``ba_jac``, final_hessian ``=`` ``"stable"``)`` `[`stopifnot`](https://rdrr.io/r/base/stopifnot.html)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(``mle_ba_analytic``$``maximum`` ``-`` ``mle_ba``$``maximum``)`` ``<`` ``2e-6``)`` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``numeric ``=`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(``mle_ba``)``, analytic ``=`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(``mle_ba_analytic``)``)`` `[`fit_diagnostics`](https://ksublee.github.io/emhawkes/reference/fit_diagnostics.md)`(``mle_ba_analytic``)``$``computation`

Numerical differentiation remains the default. The optimization path and
last digits may differ. Analytic fitting does not support arbitrary
reads of the full `lambda` or `lambda_component` histories. Branches in
this example depend only on observations. Parameter-dependent branch
thresholds need separate differentiability checks. Faster gradients do
not resolve weak identification or guarantee reliable standard errors.

## Documentation sources

The website is generated with `Rscript tools/build-docs.R` after
installing the current package and `pkgdown`. The script includes every
Rmd article in `vignettes_ext`, which is excluded from the CRAN package.
The basic vignette in `vignettes` is included in the package and
published on the website. The [GitHub Markdown
example](https://github.com/ksublee/emhawkes/blob/main/github_document/example.md)
uses the same Extended examples source. Edit the Rmd sources and
rebuild; do not edit the generated pages in `docs` independently.
