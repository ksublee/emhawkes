# Infer Intensity Paths (Deprecated)

Deprecated wrapper for
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
with `infer = TRUE`. Use
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
directly to compute intensity paths and residual processes from a
specified model.

## Usage

``` r
infer_lambda(
  object,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  lambda_component0 = NULL,
  N0 = NULL,
  Nc0 = NULL,
  ...
)

# S4 method for class 'hspec'
infer_lambda(
  object,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  lambda_component0 = NULL,
  N0 = NULL,
  Nc0 = NULL,
  ...
)
```

## Arguments

- object:

  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md).
  This object includes the parameter values.

- inter_arrival:

  inter-arrival times of events. This includes inter-arrival for events
  that occur in all dimensions. Start with zero.

- type:

  a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on.
  Start with zero.

- mark:

  a vector of mark (jump) sizes. Start with zero.

- N:

  Optional cumulative event-count matrix consistent with `type`.

- Nc:

  Optional cumulative mark-total matrix consistent with `type` and
  `mark`.

- lambda_component0:

  Initial intensity components \\\lambda\_{ij}(0)\\, excluding the
  baseline `mu`. A scalar is expanded to every kernel component. A
  matrix must match the shape of `alpha` and `beta`.

- N0:

  Initial event counts: a scalar or one value per event type.

- Nc0:

  Initial mark totals: a scalar or one value per event type.

- ...:

  further arguments passed to or from other methods.

## Value

The inference list returned by `logLik(object, infer = TRUE)`. It
contains the likelihood, intensity matrices, and residual vectors.

## See also

[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md),
[`fitted_hspec()`](https://ksublee.github.io/emhawkes/reference/fitted_hspec.md),
[`fit_initial_state()`](https://ksublee.github.io/emhawkes/reference/fit_initial_state.md)

## Examples

``` r
mu <- c(0.1, 0.1)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=100)
#> Initial intensity components: auto -> mean (unconditional component means).
summary(res)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>       arrival N1 N2 mark lambda1 lambda2
#>  [1,]  0.0000  0  0    0 0.15000 0.15000
#>  [2,]  4.5840  1  0    1 0.10081 0.10081
#>  [3,]  5.9065  1  1    1 0.16107 0.13066
#>  [4,]  7.6876  2  1    1 0.13242 0.14643
#>  [5,]  7.9299  2  2    1 0.28688 0.21774
#>  [6,] 17.5999  3  2    1 0.10005 0.10005
#>  [7,] 21.6491  4  2    1 0.10523 0.10262
#>  [8,] 21.9249  5  2    1 0.26012 0.18006
#>  [9,] 30.9103  6  2    1 0.10011 0.10006
#> [10,] 33.2391  6  3    1 0.12461 0.11230
#> [11,] 33.5206  7  3    1 0.19671 0.26478
#> [12,] 33.6337  8  3    1 0.36801 0.33917
#> [13,] 33.6810  9  3    1 0.54848 0.42502
#> [14,] 33.8535 10  3    1 0.65523 0.46389
#> [15,] 34.0128 11  3    1 0.75440 0.50196
#> [16,] 34.1304 11  4    1 0.86854 0.55152
#> [17,] 35.1472 12  4    1 0.44785 0.36093
#> [18,] 35.5445 12  5    1 0.48312 0.35241
#> [19,] 35.7998 12  6    1 0.48395 0.45954
#> [20,] 39.4787 13  6    1 0.11765 0.12041
#> ... with 80 more rows 
#> -------------------------------------------------------
inferred <- logLik(h, inter_arrival = res$inter_arrival,
                   type = res$type, infer = TRUE)
#> Initial intensity components: auto -> mean (unconditional component means).
head(inferred$lambda)
#>        lambda1   lambda2
#> [1,] 0.1500000 0.1500000
#> [2,] 0.1008077 0.1008077
#> [3,] 0.1610746 0.1306601
#> [4,] 0.1324218 0.1464282
#> [5,] 0.2868801 0.2177365
#> [6,] 0.1000477 0.1000528
```
