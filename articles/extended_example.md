# Extended examples

These examples extend the basic univariate and bivariate workflow in
[Get started](https://ksublee.github.io/emhawkes/articles/example.md).
Choose a topic below. Each article loads emhawkes and defines its own
model specifications and simulated data; the articles can be run
independently.

| Topic | Examples | Main features |
|:---|:---|:---|
| [Multi-kernel models](https://ksublee.github.io/emhawkes/articles/multi-kernel-models.md) | Multi-kernel and synchronized intensity | Multiple decay scales, `type_col_map`, shared parameters |
| [Marked models](https://ksublee.github.io/emhawkes/articles/marked-models.md) | Linear impact | `impact`, `eta`, `rmark`, `dmark` |
| [State-dependent baseline](https://ksublee.github.io/emhawkes/articles/state-dependent-baseline.md) | Fosset spread and Hawkes/M/infinity queue | State from event counts, function-valued `mu` |
| [State-dependent impact](https://ksublee.github.io/emhawkes/articles/state-dependent-impact.md) | Hawkes flocking and bid-ask prices | State-dependent excitation, `lambda_component_n` |

The baseline examples are a useful starting point for state-dependent
models. The impact examples then introduce event-dependent shocks and
adjustments to accumulated intensity components.

## Multi-kernel models

### Multi-kernel model

[Define and fit a model with multiple exponential
kernels](https://ksublee.github.io/emhawkes/articles/multi-kernel-models.html#multi-kernel-model).
Use `type_col_map` to associate kernel columns with event types.

### Synchronized intensity model

[Share parameters across component
functions](https://ksublee.github.io/emhawkes/articles/multi-kernel-models.html#synchronized-intensity-model).
This two-kernel model uses named parameters shared by `mu` and `beta`.

## Marked models

### Linear impact model

[Add mark-dependent excitation and estimate the mark
distribution](https://ksublee.github.io/emhawkes/articles/marked-models.html#linear-impact-model).
The examples cover a general `impact` function and the linear `eta`
interface.

## State-dependent baseline

### Fosset’s state-dependent spread model

[Simulate and fit spread widening and
narrowing](https://ksublee.github.io/emhawkes/articles/state-dependent-baseline.html#fosset-spread-model).
The narrowing baseline becomes zero at the minimum spread.

### A Hawkes/M/infinity queue: website visits and departures

[Model visitor arrivals and
departures](https://ksublee.github.io/emhawkes/articles/state-dependent-baseline.html#hawkes-infinite-server).
Hawkes arrivals increase the population, and each visitor has an
exponential departure time.

## State-dependent impact

### Hawkes flocking model

[Make excitation depend on relative counting-process
levels](https://ksublee.github.io/emhawkes/articles/state-dependent-impact.html#hawkes-flocking-model).
The `impact` function reads the observed event counts.

### Bid-ask price model

[Combine state-dependent baseline and component
adjustments](https://ksublee.github.io/emhawkes/articles/state-dependent-impact.html#bid-ask-price-model).
The model uses marked price changes and the current `lambda_component_n`
matrix.

## Related guides

- [Get started](https://ksublee.github.io/emhawkes/articles/example.md):
  basic univariate and bivariate models.

- [Fitting and
  inference](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html#parameter-settings):
  parameter sharing, observations and initial states.

- [Computation](https://ksublee.github.io/emhawkes/articles/computation.html#calculation-engines):
  R/C++ engines and optional analytic derivatives.

- [Residual
  diagnostics](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html#residual-process):
  checking a fitted Hawkes model.

## Point processes with flexible residuals

FRPP-specific functionality is planned to move to the separate [frpp
package](https://ksublee.github.io/frpp/). Use its [FRPP
examples](https://ksublee.github.io/frpp/articles/frpp-examples.html)
for `fspec()`, `fsim()` and `ffit()`. The full examples previously shown
here, including the original emhawkes fitting code, are preserved in the
[compatibility
vignette](https://ksublee.github.io/frpp/articles/emhawkes-compatibility.html).

Simple custom-residual models can still use
[`set_residual()`](https://ksublee.github.io/emhawkes/reference/set_residual.md),
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) and
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md) in
emhawkes. These functions remain available during the transition; this
documentation change does not remove or deprecate them. For example:

``` r

library(emhawkes)
h_residual <- set_residual(
  new("hspec", mu = 0.5, alpha = 1, beta = 1.2),
  param = c(a = 0.7, ell = 1.6),
  dresidual = dtzexp, presidual = ptzexp,
  qresidual = qtzexp, rresidual = rtzexp
)
set.seed(1107)
path_residual <- hsim(h_residual, size = 100, lambda_component0 = 0)
logLik(h_residual, inter_arrival = path_residual$inter_arrival,
       lambda_component0 = 0)
#> loglikelihood 
#>      63.13951
```

Hawkes goodness-of-fit residual diagnostics remain part of emhawkes; see
[Fitting and
inference](https://ksublee.github.io/emhawkes/articles/fitting-and-inference.html#residual-process).
