# Infer lambda process with given Hawkes model and realized path

This method compute the inferred lambda process and returns it as
`hreal` form. If we have realized path of Hawkes process and its
parameter value, then we can compute the inferred lambda processes.
Similarly with other method such as `hfit`, the input arguments are
`inter_arrival`, `type`, `mark`, or equivalently, `N` and `Nc`.

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

  Hawkes process. If not provided, then generate using inter_arrival and
  type.

- Nc:

  mark accumulated Hawkes process. If not provided, then generate using
  inter_arrival, type and mark.

- lambda_component0:

  Initial values for the lambda component \\\lambda\_{ij}\\. Can be a
  numeric value or a matrix. Must have the same number of rows and
  columns as `alpha` or `beta` in `object`.

- N0:

  the initial values of N.

- Nc0:

  the initial values of Nc.

- ...:

  further arguments passed to or from other methods.

## Value

[`hreal`](https://ksublee.github.io/emhawkes/reference/hreal.md)
S3-object, with inferred intensity.

## Examples

``` r
mu <- c(0.1, 0.1)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=100)
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
res2 <- infer_lambda(h, res$inter_arrival, res$type)
#> Warning: The 'infer_lambda' function is deprecated. Please use 'logLik' with 'infer = TRUE' instead.
#> The initial values for intensity processes are not provided. Internally determined initial values are used.
summary(res2)
#>                  Length Class  Mode   
#> loglikelihood      1    -none- numeric
#> lambda           200    -none- numeric
#> lambda_component 400    -none- numeric
#> rambda           200    -none- numeric
#> rambda_component 400    -none- numeric
#> res_process1      42    -none- numeric
#> res_process2      55    -none- numeric
```
