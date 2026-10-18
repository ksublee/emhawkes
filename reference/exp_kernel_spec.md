# Numeric Exponential Kernel Specifications

Create a baseline-independent numeric kernel specification or extract
one from an existing specification. These interfaces are experimental.

## Usage

``` r
exp_kernel_spec(alpha, beta, eta = NULL, type_col_map = NULL)
as_exp_kernel(object)
```

## Arguments

- alpha, beta:

  Finite numeric matrices of equal positive dimensions, or scalars
  representing one by one matrices. Decay rates must be positive and
  equal within each row. Rows identify event types; columns identify
  kernels.

- eta:

  Optional finite numeric matrix of the same shape.

- type_col_map:

  One list entry per event type, containing distinct integer kernel
  columns. Entries can overlap. When omitted, the number of columns must
  equal the number of event types and event type i updates column i.

- object:

  An `hspec` with numeric kernel parameters and no impact callback, or
  an existing `exp_kernel_spec`.

## Details

Construction and conversion execute no user functions and use no random
numbers. Function-valued kernel parameters and impact callbacks are
rejected. Baselines, mark generators, mark densities and residual laws
are not copied. Signed components and jumps are permitted; construction
does not establish that the corresponding model can be simulated.
Execution state is separate from the specification. Existing Hawkes
calculation paths are unchanged.

## Value

An S3 `exp_kernel_spec` containing `alpha`, `beta`, `eta`, `dimens`,
`kernel_cols` and the resolved `type_col_map`.

## See also

[`exp_kernel_replay`](https://ksublee.github.io/emhawkes/reference/exp_kernel_replay.md)

## Examples

``` r
h <- new("hspec", mu = .3, alpha = .2, beta = 1)
kernel <- as_exp_kernel(h)
exp_kernel_replay(kernel, c(0, .25, .5), lambda_component0 = .1)
#> $kernel_integral
#>            [,1]
#> [1,] 0.00000000
#> [2,] 0.02211992
#> [3,] 0.10933729
#> 
#> $kernel_before
#>            [,1]
#> [1,] 0.10000000
#> [2,] 0.07788008
#> [3,] 0.16854279
#> 
#> $final_component
#>           [,1]
#> [1,] 0.3685428
#> 
```
