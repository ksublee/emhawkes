# Compute Hawkes volatility

This function computes Hawkes volatility. Only works for bi-variate
Hawkes process.

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

  Time horizon for volatility.

- inter_arrival:

  Inter-arrival times of events which includes inter-arrival for events
  that occur in all dimensions. Start with zero.

- type:

  A vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on.
  Start with zero.

- mark:

  A vector of mark (jump) sizes. Start with zero.

- dependence:

  Dependence between mark and previous sigma-algebra.

- lambda_component0:

  A matrix of the starting values of lambda component.

- ...:

  Further arguments passed to or from other methods.
