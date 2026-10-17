# Simulate an Exponential Hawkes Process

Simulate a univariate or multivariate path from an
[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
specification, including marked and multi-kernel models.

## Usage

``` r
hsim(
  object,
  size = 100,
  lambda_component0 = "auto",
  N0 = NULL,
  Nc0 = NULL,
  verbose = FALSE,
  ...
)

# S4 method for class 'hspec'
hsim(
  object,
  size = 100,
  lambda_component0 = "auto",
  N0 = NULL,
  Nc0 = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md).
  S4 object that specifies the parameter values.

- size:

  Positive integer giving the number of output rows, including the
  initial state. At most `size - 1` events are simulated.

- lambda_component0:

  Initial intensity components \\\lambda\_{ij}(0)\\, excluding the
  baseline `mu`. A finite scalar is expanded to every kernel component.
  A matrix must match the shape of `alpha` and `beta`. The default
  `"auto"` (also `NULL`) uses unconditional component means for
  supported stable constant numeric linear models, and zero components
  otherwise. `"mean"` requires supported finite means; `"zero"` starts
  from baseline intensities. Selection is reported once, and zero
  initialization recommends supplying components explicitly. Other size
  adjustments are deprecated and produce a warning. Mean initialization
  excludes nonzero `eta`, impact callbacks and custom residual
  distributions. Function-valued core parameters use zero in `"auto"`.

- N0:

  Initial event counts: a finite scalar or one value per event type. The
  default is zero.

- Nc0:

  Initial mark totals: a finite scalar or one value per event type. The
  default is zero.

- verbose:

  Logical. Controls additional progress output; initialization selection
  is reported independently. Default is `FALSE`.

- ...:

  Further arguments passed to or from other methods.

## Value

An S3 [`hreal`](https://ksublee.github.io/emhawkes/reference/hreal.md)
object containing the simulated path, counting processes, and intensity
components. If no further event can occur, the result contains only the
generated rows and a warning reports the early termination. The initial
row is included.

## Details

The first row is the initial state. Its `inter_arrival`, `type`, and
`mark` entries are zero; its `N` and `Nc` rows contain `N0` and `Nc0`.
Event types identify intensity rows. After a type `i` event,
`type_col_map[[i]]` selects the kernel columns to update; without a
mapping, column `i` is updated.

Simulation requires finite, nonnegative baselines. Signed kernel
components are allowed when their sum is nonnegative and decreasing
between events. For equal decay rates within a row, this requires a
nonnegative row sum. For unequal rates, the simulator uses a sufficient
condition: group components by decay rate, multiply each group sum by
its rate, and check that cumulative sums in increasing rate order are
nonnegative.

Simulation automatically uses an Rcpp event loop for plain numeric `mu`,
`alpha`, and `beta`, including functions with exactly one `param`
argument resolved to such values, and positive baselines in every
dimension. This includes multivariate and mapped multi-kernel models
with equal or unequal decay rates within a row. Equal-rate models use
direct arrival sampling; unequal-rate models retain the R algorithm's
thinning and rejected draws. Equal-rate models also support numeric
`eta`, `rmark` functions requesting no event context, and `impact`
functions requesting only `lambda_component_n` or no event context.
Other optional arguments keep their defaults. Marks and impact
components must be real numeric values. User functions still run in R;
RNG state is synchronized before and after every call, including errors.
Other event callbacks, `...`, event-dependent baselines, custom residual
distributions, zero baselines, and callback models with unequal rates
use R. This selection is independent of `emhawkes.engine`, which
controls likelihood calculations. Both paths use R's RNG in the same
order. Arrival calculations use the current component state. Callbacks
receive the arguments they request; callbacks with `...` continue to
receive the full supported context. Requesting full histories can
increase copying costs for long paths.

To simulate an estimated model, use
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md)
and pass
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)
to `hsim()` with [`do.call()`](https://rdrr.io/r/base/do.call.html).
These arguments describe the start of the fitted observation window.
Continuing after its last event requires the state at that event.

## See also

[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md),
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)

## Examples

``` r

# example 1

mu <- 1; alpha <- 1; beta <- 2
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
hsim(h, size=100)
#> Initial intensity components: auto -> mean (unconditional component means).
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> An object of class "hspec" of 1-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,]    1
#> 
#> Slot alpha: 
#>      [,1]
#> [1,]    1
#> 
#> Slot beta: 
#>      [,1]
#> [1,]    2
#> 
#> Realized path :
#>       arrival N1 mark lambda1 lambda11
#>  [1,]  0.0000  0    0   2.000   1.0000
#>  [2,]  0.2301  1    1   1.631   0.6312
#>  [3,]  0.9198  2    1   1.411   0.4106
#>  [4,]  1.0301  3    1   2.131   1.1313
#>  [5,]  1.0333  4    1   3.118   2.1178
#>  [6,]  1.0601  5    1   3.955   2.9549
#>  [7,]  1.4452  6    1   2.831   1.8307
#>  [8,]  1.5987  7    1   3.083   2.0826
#>  [9,]  2.8213  8    1   1.267   0.2673
#> [10,]  2.8898  9    1   2.105   1.1049
#> [11,]  3.6405 10    1   1.469   0.4690
#> [12,]  3.7504 11    1   2.179   1.1792
#> [13,]  3.7659 12    1   3.113   2.1128
#> [14,]  4.1201 13    1   2.533   1.5327
#> [15,]  4.3839 14    1   2.494   1.4943
#> [16,]  4.5004 15    1   2.976   1.9759
#> [17,]  4.5294 16    1   3.808   2.8083
#> [18,]  4.8439 17    1   3.030   2.0303
#> [19,]  4.8767 18    1   3.838   2.8381
#> [20,]  5.4507 19    1   2.218   1.2177
#> ... with 80 more rows 
#> -------------------------------------------------------


# example 2
mu <- matrix(c(0.1, 0.1), nrow=2)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=100)
#> Initial intensity components: auto -> mean (unconditional component means).
print(res)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> An object of class "hspec" of 2-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,]  0.1
#> [2,]  0.1
#> 
#> Slot alpha: 
#>      [,1] [,2]
#> [1,]  0.2  0.1
#> [2,]  0.1  0.2
#> 
#> Slot beta: 
#>      [,1] [,2]
#> [1,]  0.9  0.9
#> [2,]  0.9  0.9
#> 
#> Realized path :
#>       arrival N1 N2 mark lambda1 lambda2  lambda11  lambda12  lambda21
#>  [1,]   0.000  0  0    0  0.1500  0.1500 3.333e-02 1.667e-02 1.667e-02
#>  [2,]   2.637  1  0    1  0.1047  0.1047 3.107e-03 1.553e-03 1.553e-03
#>  [3,]  17.416  2  0    1  0.1000  0.1000 3.397e-07 2.598e-09 1.698e-07
#>  [4,]  28.979  2  1    1  0.1000  0.1000 6.044e-06 7.851e-14 3.022e-06
#>  [5,]  31.051  2  2    1  0.1155  0.1310 9.367e-07 1.550e-02 4.684e-07
#>  [6,]  34.295  2  3    1  0.1062  0.1125 5.052e-08 6.229e-03 2.526e-08
#>  [7,]  37.816  3  3    1  0.1045  0.1089 2.126e-09 4.469e-03 1.063e-09
#>  [8,]  39.752  3  4    1  0.1358  0.1191 3.501e-02 7.825e-04 1.751e-02
#>  [9,]  40.742  3  5    1  0.1557  0.1899 1.436e-02 4.134e-02 7.181e-03
#> [10,]  40.811  4  5    1  0.2464  0.3725 1.350e-02 1.329e-01 6.751e-03
#> [11,]  51.164  4  6    1  0.1000  0.1000 1.918e-05 1.193e-05 9.589e-06
#> [12,]  53.009  4  7    1  0.1190  0.1380 3.642e-06 1.899e-02 1.821e-06
#> [13,]  53.134  4  8    1  0.2064  0.3127 3.255e-06 1.064e-01 1.628e-06
#> [14,]  56.288  4  9    1  0.1121  0.1241 1.904e-07 1.207e-02 9.521e-08
#> [15,]  63.628  5  9    1  0.1002  0.1003 2.575e-10 1.516e-04 1.288e-10
#> [16,]  66.953  5 10    1  0.1100  0.1050 1.003e-02 7.605e-06 5.017e-03
#> [17,]  68.188  5 11    1  0.1362  0.1675 3.301e-03 3.290e-02 1.651e-03
#> [18,]  68.941  5 12    1  0.1692  0.2359 1.677e-03 6.752e-02 8.386e-04
#> [19,]  68.968  5 13    1  0.2651  0.4278 1.637e-03 1.635e-01 8.184e-04
#> [20,]  69.698  5 14    1  0.2373  0.3734 8.480e-04 1.365e-01 4.240e-04
#>        lambda22
#>  [1,] 3.333e-02
#>  [2,] 3.107e-03
#>  [3,] 5.196e-09
#>  [4,] 1.570e-13
#>  [5,] 3.099e-02
#>  [6,] 1.246e-02
#>  [7,] 8.939e-03
#>  [8,] 1.565e-03
#>  [9,] 8.268e-02
#> [10,] 2.657e-01
#> [11,] 2.387e-05
#> [12,] 3.799e-02
#> [13,] 2.127e-01
#> [14,] 2.414e-02
#> [15,] 3.031e-04
#> [16,] 1.521e-05
#> [17,] 6.580e-02
#> [18,] 1.350e-01
#> [19,] 3.270e-01
#> [20,] 2.730e-01
#> ... with 80 more rows 
#> -------------------------------------------------------
```
