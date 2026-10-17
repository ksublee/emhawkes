# Compute a Residual Process (Deprecated)

Legacy helper that integrates a supplied intensity path for one event
type. Use
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
with `infer = TRUE` to compute intensities and residuals together. For a
correctly specified Hawkes model, residual intervals have a
unit-exponential distribution.

## Usage

``` r
residual_process(
  component,
  inter_arrival,
  type,
  rambda_component,
  mu,
  beta,
  dimens = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  lambda_component0 = NULL,
  N0 = NULL,
  ...
)
```

## Arguments

- component:

  The component of type to get the residual process.

- inter_arrival:

  Inter-arrival times of events. This includes inter-arrival for events
  that occur in all dimensions. Start with zero.

- type:

  A vector of types distinguished by numbers, 1, 2, 3, and so on. Start
  with zero.

- rambda_component:

  Right continuous version of lambda process.

- mu:

  Numeric value or matrix or function. If numeric, automatically
  converted to matrix.

- beta:

  Numeric value or matrix or function. If numeric, automatically
  converted to matrix, exponential decay.

- dimens:

  Dimension of the model. If omitted, set to be the length of `mu`.

- mark:

  A vector of realized mark (jump) sizes. Start with zero.

- N:

  A matrix of counting processes.

- Nc:

  A matrix of counting processes weighted by mark.

- lambda_component0:

  Initial intensity components, excluding the baseline `mu`. Must have
  the same matrix dimensions as the model's kernel parameters.

- N0:

  The initial value of N

- ...:

  Further arguments passed to or from other methods.

## Value

A numeric vector of residual intervals for the selected event type.

## Details

This function is deprecated. The replacement
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
returns residuals in `res_process1`, ..., `res_processD` and handles
types with fewer than two events by returning empty vectors. For fitted
parameters, first use
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)
and preserve initialization with
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md).

## See also

[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md),
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)

## Examples

``` r

mu <- c(0.1, 0.1)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=1000)
#> Initial intensity components: auto -> mean (unconditional component means).
inferred <- logLik(h, inter_arrival = res$inter_arrival,
                   type = res$type, infer = TRUE)
#> Initial intensity components: auto -> mean (unconditional component means).
rp <- inferred$res_process1
```
