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

  Rambda component.

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
