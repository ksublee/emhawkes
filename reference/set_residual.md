# Set Residual Distribution Functions for Hawkes Model Specification

Attach residual density, distribution, quantile, and random-generation
functions to an `hspec` object, with named parameter defaults for
simulation and estimation.

## Usage

``` r
set_residual(
  object,
  param,
  dresidual = NULL,
  presidual = NULL,
  qresidual = NULL,
  rresidual = NULL,
  ...
)

# S4 method for class 'hspec'
set_residual(
  object,
  param,
  dresidual = NULL,
  presidual = NULL,
  qresidual = NULL,
  rresidual = NULL
)
```

## Arguments

- object:

  An object of class `hspec` (Hawkes model specification)

- param:

  A named numeric vector of residual-distribution parameters. Names
  match the declared formal arguments of the supplied distribution
  functions. For estimation, the values supply starting points rather
  than fixed constants.

- dresidual:

  Density function of the residual distribution (optional)

- presidual:

  Cumulative distribution function (CDF) of the residual distribution
  (optional)

- qresidual:

  Quantile function of the residual distribution (optional)

- rresidual:

  Random generation function of the residual distribution (optional)

- ...:

  Additional arguments for future extensions

## Value

An updated `hspec` object with residual functions set

## Details

The updated specification represents a model with the supplied residual
distribution. Functions left `NULL` retain their existing slots. Supply
consistent density and CDF functions for likelihood evaluation, a
quantile function for residual inference, and a random generator for
simulation.

The `param` argument in these functions defaults to the parameters
provided during setup and is used for estimation. Only parameter names
declared by each supplied distribution function are forwarded. Matching
names across the wrapped functions share one estimated parameter.
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)
updates the residual generator when its named defaults match estimated
coefficient names. Density wrappers accept `log`; CDF and quantile
wrappers accept `lower.tail` and `log.p`. These options are forwarded
when the supplied distribution supports them. User-defined distributions
should provide log-density and upper-tail log-probability options to
retain accuracy in extreme tails.

## See also

[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md),
[`hfit()`](https://ksublee.github.io/emhawkes/reference/hfit.md),
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`dtzexp()`](https://ksublee.github.io/emhawkes/reference/tzexp.md)

## Examples

``` r
# Create basic Hawkes specification
hspec_obj <- new("hspec",
                 mu = matrix(0.1, nrow = 1),
                 alpha = matrix(0.5, nrow = 1),
                 beta = matrix(1.0, nrow = 1))

# Set residual distribution parameters
params <- c(a = 0.5, ell = 1.0)

# Apply residual functions
hspec_obj <- set_residual(
  hspec_obj,
  param = params,
  dresidual = dtzexp,
  presidual = ptzexp,
  qresidual = qtzexp,
  rresidual = rtzexp
)

# Check resulting functions
hspec_obj@dresidual
#> function (x, param = c(a = 0.5, ell = 1), log = FALSE) 
#> {
#>     residual_density(dresidual, x, log = log, params = distribution_parameters(dresidual, 
#>         param))
#> }
#> <environment: 0x5637fc730ca0>
hspec_obj@rresidual
#> function (size, param = c(a = 0.5, ell = 1)) 
#> {
#>     do.call(rresidual, c(list(size), distribution_parameters(rresidual, 
#>         param)))
#> }
#> <environment: 0x5637fc730ca0>
hspec_obj@dresidual(c(0.2, 1), log = TRUE)
#> [1] -0.1932906 -1.0007753
hspec_obj@presidual(c(0.2, 1), lower.tail = FALSE, log.p = TRUE)
#> [1] -0.1977121 -1.0007753
```
