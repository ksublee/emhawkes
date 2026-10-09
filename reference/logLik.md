# Compute the Log-Likelihood Function

Calculates the log-likelihood for the Hawkes model.

## Usage

``` r
# S4 method for class 'hspec'
logLik(
  object,
  inter_arrival,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  N0 = NULL,
  Nc0 = NULL,
  lambda_component0 = NULL,
  infer = FALSE,
  ...
)
```

## Arguments

- object:

  An
  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
  object containing parameter values for computing the log-likelihood.

- inter_arrival:

  A vector of inter-arrival times for events across all dimensions,
  starting with zero.

- type:

  A vector indicating the dimensions, represented by numbers (1, 2, 3,
  etc.), starting with zero.

- mark:

  A vector of mark (jump) sizes, starting with zero.

- N:

  A matrix representing counting processes.

- Nc:

  A matrix of counting processes weighted by mark sizes.

- N0:

  A matrix of initial values for `N`.

- Nc0:

  A matrix of initial values for `Nc`.

- lambda_component0:

  Initial values for the lambda component \\\lambda\_{ij}\\. Can be a
  numeric value or a matrix. Must have the same number of rows and
  columns as `alpha` or `beta` in `object`.

- infer:

  Logical

- ...:

  Additional arguments passed to or from other methods.

## Details

The first row represents the initial state and is not an event.
Counting-process rows contain values after the corresponding event.
Supplied `N` and `Nc` matrices must agree with `type`, `mark`, and any
supplied initial counts. Event order is preserved, including zero
intervals. Component functions receive the arguments they declare;
functions with `...` receive the available event context. Intensity
histories are `NULL` before initialization. Baseline functions use
`n - 1` to access the preceding event.

## See also

[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md),
[`hfit,hspec-method`](https://ksublee.github.io/emhawkes/reference/hfit.md)
