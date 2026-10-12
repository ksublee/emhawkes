# Retrieve Initial-State Arguments from a Hawkes Fit

Retrieve the normalized initial-state arguments saved by
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

## Usage

``` r
fit_initial_state(fit)
```

## Arguments

- fit:

  A result from the current
  [`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md).

## Value

A list with `lambda_component0`, `N0`, and `Nc0`, suitable for passing
to [`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
or [`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md) with
[`do.call()`](https://rdrr.io/r/base/do.call.html).

## Details

A supplied initial intensity component is stored as a matrix. If it was
omitted, `lambda_component0` remains `NULL`, preserving automatic
initialization with the restored model's parameters. It is not a fixed
matrix evaluated at the original starting values.

Initial counts are normalized per dimension and include values inferred
from supplied `N` and `Nc` paths. With a custom likelihood and no event
observations, the count arguments are kept as supplied. These arguments
describe the initial state of the fitted observation window, not the
state after its last event.

## See also

[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md),
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md),
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md)
