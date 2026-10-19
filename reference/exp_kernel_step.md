# Event-Wise Numeric Exponential Kernel Calculations

Advance kernel components, apply an optional event update, or invert
supplied integrated-driver targets. These experimental deterministic
interfaces generate no random numbers and do not run a simulation loop.

## Usage

``` r
exp_kernel_step(kernel, component, dt, type = NULL, mark = 1,
                integrate = TRUE)
exp_kernel_wait(kernel, component, mu, residual)
```

## Arguments

- kernel:

  An `exp_kernel_spec` or convertible numeric `hspec`. Decay rates must
  be equal within each row.

- component:

  Finite current kernel components, excluding baseline. A scalar expands
  to all components; a matrix matches the kernel shape, and a vector of
  matching length uses R's column-major order.

- dt:

  A finite nonnegative scalar interval.

- type:

  An event-type row index whose mapped columns are updated, or `NULL` to
  advance without applying an event.

- mark:

  A finite scalar event mark. The default is one.

- integrate:

  Logical; also return the component-wise interval integrals?

- mu:

  Finite nonnegative baseline, scalar or one value per event type.

- residual:

  Nonnegative integrated-driver targets, one per event type. Positive
  infinity is allowed. Targets are supplied by the calling model.

## Details

The specification contains no evolving state. Step returns the pre-event
components after decay and the post-event components after mapped alpha
addition followed by eta times (mark minus one). Signed components are
not clipped. Baseline contributions are excluded from step integrals. A
zero interval still updates when `type` is supplied.

Wait inverts, independently for each row i, \$\$A_i(t) = \mu_i t + C_i
(1-\exp(-\beta_i t))/\beta_i,\$\$ where C is the row sum of the current
components. Kernel row sums must be nonnegative and decreasing; the
existing simulation cancellation guard is retained for signed
components. A zero target yields zero. With zero baseline, a target at
or above the finite total integral yields infinity. No event is selected
and no residual distribution is evaluated by this function. Summation
and root-solving order match the existing equal-rate FRPP arrival
calculation; component integrals are not separately summed during
inversion.

The caller owns residual generation, event selection, initial-state
policy, counts, output recording and callback behavior. A simulator must
check its current state before drawing new residuals and must not
restart another backend after consuming draws. These functions do not
change the existing Hawkes or legacy FRPP paths, or the meaning of
`emhawkes.engine`.

## Value

Step returns a list of `before` and `after` component matrices, and
`integral` when requested. Wait returns a numeric candidate waiting-time
vector in event-type row order. Inputs are not modified.

## See also

[`exp_kernel_spec`](https://ksublee.github.io/emhawkes/reference/exp_kernel_spec.md),
[`exp_kernel_replay`](https://ksublee.github.io/emhawkes/reference/exp_kernel_replay.md)

## Examples

``` r
kernel <- exp_kernel_spec(.2, 1)
wait <- exp_kernel_wait(kernel, component = .1, mu = .3, residual = .7)
exp_kernel_step(kernel, component = .1, dt = wait, type = 1)
#> $before
#>            [,1]
#> [1,] 0.01296127
#> 
#> $after
#>           [,1]
#> [1,] 0.2129613
#> 
#> $integral
#>            [,1]
#> [1,] 0.08703873
#> 
```
