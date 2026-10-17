# Compute Hawkes volatility

This function computes Hawkes volatility for a constant bivariate Hawkes
model with one kernel per event type and a common decay rate within each
row.

## Usage

``` r
hvol(
  object,
  horizon = 1,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  dependence = FALSE,
  lambda_component0 = NULL,
  ...
)

# S4 method for class 'hspec'
hvol(
  object,
  horizon = 1,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  dependence = FALSE,
  lambda_component0 = NULL,
  ...
)
```

## Arguments

- object:

  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)

- horizon:

  Time horizon used to scale the stationary long-run variance rate.

- inter_arrival:

  Inter-arrival times of events which includes inter-arrival for events
  that occur in all dimensions. Start with zero.

- type:

  A vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on.
  Start with zero.

- mark:

  A vector of mark (jump) sizes. Start with zero. If omitted, all
  observed events have unit marks. At least one event of each type is
  required.

- dependence:

  Dependence between mark and previous sigma-algebra.

- lambda_component0:

  A matrix of initial intensity components, excluding the baseline `mu`.

- ...:

  Further arguments passed to or from other methods.

## Details

Numeric `mu`, `alpha`, and `beta` are required. A missing `eta` is
treated as zero. Custom impact functions, residual distributions, and
noncanonical kernel mappings are not supported by this moment formula.
The result is the square root of `horizon` times the stationary long-run
variance rate of the difference between the two marked counting
processes.
