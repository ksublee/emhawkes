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

Supplied initial intensity components exclude the baseline `mu` and are
stored as a matrix. If the argument was omitted, the selected mean rule
is stored as `"mean"`, or the selected zero state as a matrix. Mean
components are recomputed with the restored model's parameters; the
automatic selection is not repeated. This is not a matrix evaluated at
the original starting values. Older saved fits without an
initialization-policy record return a matrix evaluated under their
original initialization rule at the fitted parameters. Their retained
objective supplies any observation context needed by callbacks.

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
