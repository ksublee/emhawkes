# Shared exponential kernel calculations

## Purpose and availability

The GitHub development version 1.0.1.9000 exports
[`exp_kernel_spec()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md),
[`as_exp_kernel()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md)
and
[`exp_kernel_replay()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_replay.md).
These experimental functions provide a shared calculation boundary
inside emhawkes. The separate [frpp
package](https://ksublee.github.io/frpp/) uses it for its numeric
`state` likelihood. The dependency is frpp -\> emhawkes.

The API supplies kernel states and integrals. The calling model owns its
baseline, initialization policy, residual distribution, likelihood and
random draws. Existing Hawkes event loops and `emhawkes.engine`
selection are unchanged. Publication on GitHub does not update an
already installed CRAN version.

## Specification and mapping

Numeric alpha and beta are D by K matrices; eta is optional with the
same shape. Beta must be positive and equal within each row. Rows are
event types and columns are kernels. `type_col_map[[i]]` selects the
columns updated by type i; entries can overlap across event types.
Signed components and jumps are retained. These shape rules alone do not
guarantee a valid simulator or stationarity. Function-valued kernel
parameters and impact callbacks are outside this first interface;
conversion rejects them without executing them.

``` r

h <- new("hspec", mu = c(.3, .4),
  alpha = matrix(c(.04, .03, .02, .05, .01, .02, .03, .01), 2, 4),
  beta = matrix(rep(c(1.2, 1.7), 4), 2),
  type_col_map = list(c(3L, 1L), c(2L, 4L)))
kernel <- as_exp_kernel(h)
```

[`as_exp_kernel()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md)
reads kernel slots only. Baselines, mark generators, mark densities and
residual laws are not copied or evaluated. The specification contains no
evolving state, so replaying it twice does not accumulate events.

## Replay contract

Supply observed intervals, types, optional marks and explicit initial
components. The first row is the initial state: its interval is zero and
it causes no event update. Subsequent zero intervals update events in
input order. Replay adds no terminal observation interval and consumes
no random numbers.

``` r

dt <- c(0, .2, .4, 0, .3)
type <- c(0, 1, 2, 1, 2)
initial <- matrix(.05, 2, 4)
x <- exp_kernel_replay(kernel, dt, type, lambda_component0 = initial)
x$kernel_integral
#>            [,1]       [,2]
#> [1,] 0.00000000 0.00000000
#> [2,] 0.03556202 0.03390937
#> [3,] 0.06586329 0.05582602
#> [4,] 0.00000000 0.00000000
#> [5,] 0.05751446 0.04875125
x$kernel_before
#>           [,1]       [,2]
#> [1,] 0.2000000 0.20000000
#> [2,] 0.1573256 0.14235406
#> [3,] 0.1282896 0.09744984
#> [4,] 0.1782896 0.15744984
#> [5,] 0.1592723 0.12457271
x$final_component
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.06215525 0.05093330 0.02827364 0.06791007
#> [2,] 0.03796829 0.09085156 0.02892112 0.02683174
```

For M input rows, the ordinary output is:

| Output | Shape | Meaning |
|----|----|----|
| `kernel_integral` | M by D | Kernel-only integral over each interval; initial row is zero |
| `kernel_before` | M by D | Pre-event component sums; initial row is the initial row sums |
| `final_component` | D by K | State after the final event, or unchanged initial state for M = 1 |

A scalar initial value expands over all components. A matrix must match
D by K; a matching vector follows R’s column-major matrix order.
Automatic mean/zero policies and count paths belong to the calling
model.

With `history = TRUE`, `component_before` and `component_after` are
additionally returned as M by (D times K). Within each row, components
are flattened row by row. Histories are allocated only when requested.
Ordinary output still stores M by D interval integrals and sums; it is
not constant-memory streaming.

``` r

full <- exp_kernel_replay(kernel, dt, type, lambda_component0 = initial,
  history = TRUE)
stopifnot(identical(x, full[names(x)]))
```

## Using the result in another model

For a numeric baseline mu, the caller forms pre-event drivers and
integrated drivers as follows:

``` r

mu <- as.vector(h@mu)
driver_before <- sweep(x$kernel_before, 2, mu, "+")
integrated_driver <- x$kernel_integral + outer(dt, mu)
driver_before
#>           [,1]      [,2]
#> [1,] 0.5000000 0.6000000
#> [2,] 0.4573256 0.5423541
#> [3,] 0.4282896 0.4974498
#> [4,] 0.4782896 0.5574498
#> [5,] 0.4592723 0.5245727
integrated_driver
#>            [,1]      [,2]
#> [1,] 0.00000000 0.0000000
#> [2,] 0.09556202 0.1139094
#> [3,] 0.18586329 0.2158260
#> [4,] 0.00000000 0.0000000
#> [5,] 0.14751446 0.1687513
```

FRPP combines these with residual density and survival terms. A
non-exponential residual law changes the conditional hazard: the kernel
driver by itself is not automatically the FRPP conditional intensity.
Hawkes-specific stationary means, gradients and Hessian adjustments are
not supplied by this API.

Use the exported functions rather than `emhawkes:::`. The implementation
is an R reference path, independent of the likelihood engine option.
Internal advance/update helpers and C++ implementation details are not
part of this public contract. Event-wise calculations below support
independent FRPP event generation; callback kernels remain outside this
numeric interface.

## Event-wise calculations

[`exp_kernel_wait()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md)
inverts supplied integrated-driver targets using a numeric baseline and
current kernel components. It draws no residuals and selects no event.
[`exp_kernel_step()`](https://ksublee.github.io/emhawkes/reference/exp_kernel_step.md)
advances the components and optionally applies one mapped event; it also
returns component integrals unless `integrate = FALSE`. Neither function
replays a growing path prefix.

``` r

component <- initial
candidates <- exp_kernel_wait(kernel, component, mu, residual = c(.7, .9))
event <- which.min(candidates)
next_state <- exp_kernel_step(kernel, component, min(candidates), type = event)
next_state$before
#>             [,1]        [,2]        [,3]        [,4]
#> [1,] 0.005503118 0.005503118 0.005503118 0.005503118
#> [2,] 0.002194278 0.002194278 0.002194278 0.002194278
next_state$after
#>            [,1]        [,2]       [,3]        [,4]
#> [1,] 0.04550312 0.005503118 0.01550312 0.005503118
#> [2,] 0.03219428 0.002194278 0.02219428 0.002194278
next_state$integral
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.03708074 0.03708074 0.03708074 0.03708074
#> [2,] 0.02812101 0.02812101 0.02812101 0.02812101
```

A zero interval still updates the selected event. `type = NULL` performs
only decay and integration. Signed components are retained. Waiting-time
inversion requires admissible decreasing row sums; a finite total driver
integral can produce an infinite waiting time when its baseline is zero.

For repeated calculations with fixed parameters, prepare a reusable
plan:

``` r

plan <- exp_kernel_plan(kernel, mu)
stopifnot(identical(candidates, plan$wait(component, c(.7, .9))))
stopifnot(identical(next_state, plan$step(component, min(candidates), event)))
```

The plan captures a validated kernel and baseline once. Its `wait()` and
[`step()`](https://rdrr.io/r/stats/step.html) functions still check the
changing components and arguments. Components remain explicit inputs:
the plan stores no event state, history or RNG. It is reusable across
independent paths; create a new plan when parameters change. Editing the
original specification does not change an existing plan. The plan’s
function bindings are locked. Older event-wise functions remain
available with their existing checked input contracts.

FRPP owns its residual-generation order, earliest-event selection, count
paths, output histories and RNG handling. The event-wise API receives
explicit targets and contains no model-specific residual law. Legacy
Hawkes and FRPP event loops are unchanged. This is an R reference
interface, not a speed improvement claim.
