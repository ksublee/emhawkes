# Register and Restore Model Parameters

Prepare model-independent parameter coordinates and component
reconstruction.

## Usage

``` r
parameter_plan(components, prefixes = names(components),
               named = character(), reduced = TRUE)
```

## Arguments

- components:

  Named list of numeric matrices/scalars, functions or NULL.

- prefixes:

  Matrix name prefixes, a character vector matching components.

- named:

  Names of explicit named parameter-vector components. Their names are
  retained; equal starting values do not imply sharing.

- reduced:

  Share equal starting values within each matrix?

## Value

A locked environment with `start`, the full starting vector, and
`restore(param)`, returning a named list of matrices, functions and
vectors.

## Details

Components are visited in list order. Matrix naming and sharing match
[`hfit`](https://ksublee.github.io/emhawkes/reference/hfit.md). Function
parameters come from their named numeric `param` defaults, evaluated in
the original environment; function bodies are not called. Names shared
across components identify one parameter. Conflicting starts warn and
use the first value. Functions retain their environments and binding
semantics. No model, optimizer, initialization policy or random
generator is supplied by this API. External function environments remain
references; evaluating user parameter defaults can have side effects.

Restore accepts the full vector in starting order, or exactly its names
in any order. Nonfinite proposals are retained for the objective to
handle. Fixed/free selectors and constraints use this full coordinate
order and are handled by the optimizer layer. Create a new plan when
definitions change.

## Examples

``` r
plan <- parameter_plan(list(mu = matrix(c(.2, .2), 2),
  residual = c(shape = 2, rate = 2)), named = "residual")
plan$start
#>   mu1 shape  rate 
#>   0.2   2.0   2.0 
plan$restore(plan$start)
#> $mu
#>      [,1]
#> [1,]  0.2
#> [2,]  0.2
#> 
#> $residual
#> shape  rate 
#>     2     2 
#> 
```
