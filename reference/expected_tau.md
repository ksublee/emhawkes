# Expected Inter-Arrival Time

Computes the conditional expected time until the next event.

## Usage

``` r
expected_tau(
  object,
  rambda_component,
  type = 1,
  mu = NULL,
  beta = NULL,
  tol = .Machine$double.eps^0.25,
  max_upper = Inf,
  subdivisions = 400L
)
```

## Arguments

- object:

  An object of class `hspec`.

- rambda_component:

  Finite nonnegative scalar: the current total excitation for the
  selected process dimension (sum its kernel components first).

- type:

  Process dimension index (default is 1).

- mu:

  Optional mu value (overrides object@mu if provided).

- beta:

  Optional beta value (overrides object@beta if provided).

- tol:

  Relative tolerance for numerical integration.

- max_upper:

  Upper integration limit.

- subdivisions:

  Number of subdivisions for numerical integration.

## Value

Expected value of next inter-arrival time.

## Details

This calculation assumes a constant baseline and a common decay rate
within the selected row. With zero baseline, the result is `Inf` if
there is positive probability of no further event. A custom residual law
then requires `presidual` to determine this probability. A finite
`max_upper` instead returns the truncated first moment, excluding
waiting times above that limit.
