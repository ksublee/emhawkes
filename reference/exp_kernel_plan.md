# Prepare Reusable Numeric Exponential Kernel Calculations

Validate and capture a numeric exponential kernel and baseline once for
repeated deterministic calculations.

## Usage

``` r
exp_kernel_plan(kernel, mu)
```

## Arguments

- kernel:

  An exp_kernel_spec or convertible numeric hspec.

- mu:

  Finite nonnegative baseline, scalar or one value per event type.

## Value

A locked environment with two functions:
`step(component, dt, type = NULL, mark = 1, integrate = TRUE)` and
`wait(component, residual)`. Their results match
[`exp_kernel_step`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md)
and
[`exp_kernel_wait`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md).

## Details

Later edits to the original inputs do not change the plan. Every call
still checks the current component and changing arguments. Components
are passed explicitly; no current state, histories or random generator
are stored or advanced. Create a new plan when kernel parameters or
baseline change. This experimental numeric interface requires equal beta
values within each row. The likelihood engine option and existing event
loops are unchanged.

## Examples

``` r
kernel <- exp_kernel_spec(.2, 1)
plan <- exp_kernel_plan(kernel, mu = .3)
component <- matrix(.1)
dt <- plan$wait(component, residual = .7)
next_state <- plan$step(component, dt, type = 1)
next_state$after
#>           [,1]
#> [1,] 0.2129613
```
