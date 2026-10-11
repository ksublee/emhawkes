# Perform Maximum Likelihood Estimation

Estimate exponential Hawkes model parameters by maximum likelihood. The
supplied
[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
object defines the model and the starting values for optimization with
[`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html).

## Usage

``` r
hfit(
  object,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  lambda_component0 = NULL,
  N0 = NULL,
  mylogLik = NULL,
  reduced = TRUE,
  grad = NULL,
  hess = NULL,
  constraints = NULL,
  method = "BFGS",
  verbose = FALSE,
  ...,
  Nc0 = NULL
)

# S4 method for class 'hspec'
hfit(
  object,
  inter_arrival = NULL,
  type = NULL,
  mark = NULL,
  N = NULL,
  Nc = NULL,
  lambda_component0 = NULL,
  N0 = NULL,
  mylogLik = NULL,
  reduced = TRUE,
  grad = NULL,
  hess = NULL,
  constraints = NULL,
  method = "BFGS",
  verbose = FALSE,
  ...,
  Nc0 = NULL
)
```

## Arguments

- object:

  An
  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
  object containing the parameter values.

- inter_arrival:

  A vector of inter-arrival times for events across all dimensions,
  starting with zero.

- type:

  A vector indicating the dimensions, represented by numbers like 1, 2,
  3, etc., starting with zero.

- mark:

  A vector of mark (jump) sizes, starting with zero.

- N:

  A matrix representing counting processes.

- Nc:

  A matrix of counting processes weighted by mark sizes.

- lambda_component0:

  Initial values for the lambda component \\\lambda\_{ij}\\. Can be a
  numeric value or a matrix. Must have the same number of rows and
  columns as `alpha` or `beta` in `object`.

- N0:

  Initial values for the counting processes matrix `N`.

- mylogLik:

  A user-defined log-likelihood function, which must accept an `object`
  argument consistent with `object`.

- reduced:

  Logical; if `TRUE`, equal entries within each parameter matrix share
  one parameter. If `FALSE`, each matrix entry is estimated separately.
  Function components share parameters by name. Conflicting starting
  values for a shared name produce a warning; the first value is used.

- grad:

  A gradient matrix for the likelihood function. Refer to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) for more
  details.

- hess:

  A Hessian matrix for the likelihood function. Refer to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) for more
  details.

- constraints:

  constraints matrices. Refer to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) for more
  details.

- method:

  The optimization method to be used. Refer to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) for more
  details.

- verbose:

  Logical; if `TRUE`, prints the progress of the estimation process.

- ...:

  Additional parameters for optimization. Refer to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) for more
  details.

- Nc0:

  Initial values for the mark-weighted counting processes. This argument
  must be named.

## Value

[`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) object

## Details

Observations are checked once before optimization. The first row is the
initial state, and supplied counting processes contain values after each
event. Functions with a `param` argument must provide a default named
numeric vector; `numeric(0)` and functions without `param` are allowed
for fixed components. A custom `mylogLik` receives its declared
arguments from `object`, the observation arguments, `lambda_component0`,
`param`, and `reduced`. Optional defaults are preserved; `...` receives
the remaining context. Additional inputs should be captured in the
function's environment or given defaults.

When equal matrix entries are grouped, a message lists the affected
components once before optimization. Use `suppressMessages(hfit(...))`
to hide it. Setting `reduced = FALSE` separates matrix entries; function
parameters still share values by name. A numeric zero is a starting
value, not a fixed zero. To fix entries at zero, use a component
function with constants at those positions.

## See also

[`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md),
[`hsim,hspec-method`](https://ksublee.github.io/emhawkes/reference/hsim.md)

## Examples

``` r

# example 1
mu <- c(0.1, 0.1)
alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
res <- hsim(h, size=100)
summary(hfit(h, inter_arrival=res$inter_arrival, type=res$type))
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).
#> Use reduced = FALSE to estimate matrix entries separately.
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 29 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -245.0405 
#> 4  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.10357    0.01918   5.400 6.67e-08 ***
#> alpha1.1  0.21600    0.09858   2.191  0.02845 *  
#> alpha1.2  0.29650    0.11754   2.523  0.01165 *  
#> beta1.1   1.11245    0.34199   3.253  0.00114 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------


# example 2
# \donttest{
mu <- matrix(c(0.08, 0.08, 0.05, 0.05), nrow = 4)
alpha <- function(param = c(alpha11 = 0, alpha12 = 0.4, alpha33 = 0.5, alpha34 = 0.3)){
  matrix(c(param["alpha11"], param["alpha12"], 0, 0,
           param["alpha12"], param["alpha11"], 0, 0,
           0, 0, param["alpha33"], param["alpha34"],
           0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = TRUE)
}
beta <- matrix(c(rep(0.6, 8), rep(1.2, 8)), nrow = 4, byrow = TRUE)

impact <- function(param = c(alpha1n=0, alpha1w=0.2, alpha2n=0.001, alpha2w=0.1),
                   n=n, N=N, ...){

  Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
                  0, 0, param['alpha1n'], param['alpha1w'],
                  param['alpha2w'], param['alpha2n'], 0, 0,
                  param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=TRUE)

  ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5

  km <- matrix(c(!ind, !ind, !ind, !ind,
                 ind, ind, ind, ind,
                 ind, ind, ind, ind,
                 !ind, !ind, !ind, !ind), nrow = 4, byrow = TRUE)

  km * Psi
}
h <- new("hspec",
         mu = mu, alpha = alpha, beta = beta, impact = impact)
hr <- hsim(h, size=100)
plot(hr$arrival, hr$N[,'N1'] - hr$N[,'N2'], type='s')
lines(hr$N[,'N3'] - hr$N[,'N4'], type='s', col='red')

fit <- hfit(h, hr$inter_arrival, hr$type)
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, beta).
#> Use reduced = FALSE to estimate matrix entries separately.
summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 146 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -158.4592 
#> 12  free parameters
#> Estimates:
#>         Estimate Std. error t value  Pr(> t)    
#> mu1      0.69224    0.42677   1.622   0.1048    
#> mu3      0.07379    0.03508   2.103   0.0354 *  
#> alpha11 -0.12631    0.06528  -1.935   0.0530 .  
#> alpha12  0.01429    0.04928   0.290   0.7719    
#> alpha33  0.96110    0.21444   4.482 7.40e-06 ***
#> alpha34  0.29568    0.11618   2.545   0.0109 *  
#> beta1.1  0.10274    0.08867   1.159   0.2466    
#> beta3.1  1.33410    0.22745   5.866 4.48e-09 ***
#> alpha1n  0.03456    0.06731   0.513   0.6076    
#> alpha1w  0.04976    0.06841   0.727   0.4670    
#> alpha2n -2.47298    0.27532  -8.982  < 2e-16 ***
#> alpha2w  1.00872    0.11327   8.906  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
# }

# example 3
# \donttest{
mu <- c(0.15, 0.15)
alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=TRUE)
beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2, byrow=TRUE)
rmark <- function(param = c(p=0.65), ...){
  rgeom(1, p=param[1]) + 1
}
impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
  ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
  alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
}
h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
          rmark = rmark,
          impact=impact)
res <- hsim(h1, size=100, lambda_component0 = matrix(rep(0.1,4), nrow=2))

fit <- hfit(h1,
            inter_arrival = res$inter_arrival,
            type = res$type,
            mark = res$mark,
            lambda_component0 = matrix(rep(0.1,4), nrow=2))
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).
#> Use reduced = FALSE to estimate matrix entries separately.
summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 37 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -191.1433 
#> 5  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.10755    0.02344   4.589 4.46e-06 ***
#> alpha1.1  0.62799    0.25102   2.502  0.01236 *  
#> alpha1.2  0.55832    0.21978   2.540  0.01107 *  
#> beta1.1   2.11657    0.73318   2.887  0.00389 ** 
#> eta1      0.10179    0.27560   0.369  0.71187    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
# }
# For more information, please see vignettes.
```
