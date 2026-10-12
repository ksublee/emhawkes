# Simulate an Exponential Hawkes Process

Simulate a univariate or multivariate path from an
[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
specification, including marked and multi-kernel models.

## Usage

``` r
hsim(
  object,
  size = 100,
  lambda_component0 = NULL,
  N0 = NULL,
  Nc0 = NULL,
  verbose = FALSE,
  ...
)

# S4 method for class 'hspec'
hsim(
  object,
  size = 100,
  lambda_component0 = NULL,
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

  Initial values for the lambda component \\\lambda\_{ij}\\. A finite
  scalar is expanded to every kernel component. A matrix must match the
  shape of `alpha` and `beta`. Supply this argument when `type_col_map`
  is set. Other size adjustments are deprecated and produce a warning.

- N0:

  Initial event counts: a finite scalar or one value per event type. The
  default is zero.

- Nc0:

  Initial mark totals: a finite scalar or one value per event type. The
  default is zero.

- verbose:

  Logical. If `TRUE`, print progress messages during the simulation.
  Default is `FALSE`.

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

Simulation uses R. To simulate an estimated model, use
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
#>  [1,]  0.0000  0    0   1.500   0.5000
#>  [2,]  0.6689  1    1   1.131   0.1312
#>  [3,]  1.3586  2    1   1.285   0.2848
#>  [4,]  1.4812  3    1   2.005   1.0054
#>  [5,]  1.4846  4    1   2.992   1.9920
#>  [6,]  1.5126  5    1   3.829   2.8290
#>  [7,]  1.9171  6    1   2.705   1.7049
#>  [8,]  2.0790  7    1   2.957   1.9567
#>  [9,]  3.3016  8    1   1.256   0.2564
#> [10,]  3.3701  9    1   2.095   1.0954
#> [11,]  4.1209 10    1   1.467   0.4669
#> [12,]  4.2309 11    1   2.177   1.1771
#> [13,]  4.2464 12    1   3.111   2.1107
#> [14,]  4.6010 13    1   2.531   1.5307
#> [15,]  4.8651 14    1   2.492   1.4923
#> [16,]  4.9816 15    1   2.974   1.9743
#> [17,]  5.0106 16    1   3.807   2.8067
#> [18,]  5.3251 17    1   3.029   2.0295
#> [19,]  5.3578 18    1   3.837   2.8372
#> [20,]  5.9318 19    1   2.217   1.2174
#> ... with 80 more rows 
#> -------------------------------------------------------


# example 2
mu <- matrix(c(0.1, 0.1), nrow=2)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=100)
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
