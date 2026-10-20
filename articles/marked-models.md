# Marked models

This article introduces mark-dependent excitation through `impact` and
`eta`, then estimates a mark distribution using `dmark`.

See the [examples
overview](https://ksublee.github.io/emhawkes/articles/extended_example.md)
for other model topics.

``` r

library(emhawkes)
```

The package supports the following family of extended multivariate
marked Hawkes models:

``` math
 \lambda(t) = \mu + \int_{(-\infty,t)\times E} h(t, u, z)M(du \times dz)
```

where the kernel $`h`$ is

``` math
h(t, u, z) = (\alpha + g(t, z))\Gamma(t),
```

The components are:

- $`\alpha`$ is a constant matrix.

- $`g(t, z)`$ describes additional effects on the intensities and may
  depend on the mark or the history of the underlying processes.

- $`\Gamma(t)`$ is a matrix of exponential decay terms, with
  $`\Gamma_{ij}(t) = e^{-\beta_{ij}t}`$.

- $`M`$ is a vector of random measures on the product of the time and
  mark spaces.

## Linear impact model

In the linear impact model,

``` math
 g(t, z) = \eta (z-1).
```

The `impact` function represents $`g(t, z)`$, the additional effect of a
mark on future intensity. For details, see [Marked Hawkes process
modeling of price dynamics and volatility
estimation](https://www.sciencedirect.com/science/article/abs/pii/S0927539816300810).

Its first argument, `param`, contains the model parameters. Additional
arguments, such as `n` and `mark`, provide information about the model
or realized path. Include `...` to accept other arguments supplied by
the package. The `rmark` function generates marks during simulation.

``` r

mu <- matrix(c(0.15, 0.15), nrow=2)
alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=T)
beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2)
rmark <- function(param = c(p=0.65), ...){
  rgeom(1, p=param[1]) + 1
}

impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
  ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
  ma * matrix( rep(param["eta1"], 4), nrow=2)
}

hi <- new("hspec", mu=mu, alpha=alpha, beta=beta,
          rmark = rmark,
          impact=impact)
hi
#> An object of class "hspec" of 2-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,] 0.15
#> [2,] 0.15
#> 
#> Slot alpha: 
#>      [,1] [,2]
#> [1,] 0.75 0.60
#> [2,] 0.60 0.75
#> 
#> Slot beta: 
#>      [,1] [,2]
#> [1,]  2.6  2.6
#> [2,]  2.6  2.6
#> 
#> Slot impact: 
#> function (param = c(eta1 = 0.2), alpha, n, mark, ...) 
#> {
#>     ma <- matrix(rep(mark[n] - 1, 4), nrow = 2)
#>     ma * matrix(rep(param["eta1"], 4), nrow = 2)
#> }
#> 
#> Slot rmark: 
#> function (param = c(p = 0.65), ...) 
#> {
#>     rgeom(1, p = param[1]) + 1
#> }
```

``` r

set.seed(1107)
res_impact <- hsim(hi, size=1000, lambda_component0 = matrix(rep(0.1,4), nrow=2))
summary(res_impact)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>       arrival N1 N2 mark lambda1 lambda2
#>  [1,]  0.0000  0  0    0 0.35000 0.35000
#>  [2,]  1.3001  0  1    1 0.15681 0.15681
#>  [3,]  2.5585  1  1    2 0.17302 0.17871
#>  [4,]  2.7982  1  2    1 0.67171 0.59433
#>  [5,]  3.0525  2  2    2 0.72904 0.76653
#>  [6,]  3.3929  3  2    1 0.78108 0.73464
#>  [7,]  3.4509  4  2    3 1.33764 1.16872
#>  [8,]  3.5101  5  2    2 2.15429 1.88085
#>  [9,]  3.5837  6  2    1 2.58987 2.24016
#> [10,]  3.6041  7  2    1 3.17455 2.70074
#> [11,]  3.7053  8  2    1 3.05159 2.57205
#> [12,]  3.9277  9  2    1 2.19788 1.84482
#> [13,]  3.9610 10  2    1 2.71584 2.25450
#> [14,]  4.2448 11  2    1 1.73576 1.44339
#> [15,]  4.7277 12  2    4 0.81551 0.68947
#> [16,]  5.1473 12  3    1 0.82698 0.73426
#> [17,]  5.3229 12  4    2 0.95877 0.99505
#> [18,]  5.4219 13  4    1 1.39371 1.53772
#> [19,]  5.4869 14  4    2 1.83388 1.82882
#> [20,]  5.8857 14  5    1 1.08385 1.02887
#> ... with 980 more rows 
#> -------------------------------------------------------
```

``` r

fit <- hfit(hi,
            inter_arrival = res_impact$inter_arrival,
            type = res_impact$type,
            mark = res_impact$mark,
            lambda_component0 = matrix(rep(0.1,4), nrow=2))
#> Registered S3 methods overwritten by 'maxLik':
#>   method         from    
#>   coef.maxLik    emhawkes
#>   logLik.maxLik  emhawkes
#>   print.maxLik   emhawkes
#>   summary.maxLik emhawkes
#>   vcov.maxLik    emhawkes
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).
#> Use reduced = FALSE to estimate matrix entries separately.

summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 38 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1525.501 
#> 5  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1      0.131904   0.008914  14.798  < 2e-16 ***
#> alpha1.1 0.676090   0.080586   8.390  < 2e-16 ***
#> alpha1.2 0.603453   0.075153   8.030 9.77e-16 ***
#> beta1.1  2.221258   0.173268  12.820  < 2e-16 ***
#> eta1     0.116345   0.055318   2.103   0.0354 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

For linear mark effects, specify the `eta` slot directly. The resulting
intensity process is

``` math
 \lambda(t) = \mu + \int_{(-\infty, t)\times E} (\alpha + \eta (z-1)) e^{-\beta(t-u)}  M(du \times dz).
```

``` r


rmark <- function(param = c(p=0.65), ...){
  rgeom(1, p=param[1]) + 1
}

h <-  new("hspec", mu=0.15, alpha=0.7, beta=1.6, eta=0.3,
          rmark = rmark)
h
#> An object of class "hspec" of 1-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,] 0.15
#> 
#> Slot alpha: 
#>      [,1]
#> [1,]  0.7
#> 
#> Slot beta: 
#>      [,1]
#> [1,]  1.6
#> 
#> Slot eta: 
#>      [,1]
#> [1,]  0.3
#> 
#> Slot rmark: 
#> function (param = c(p = 0.65), ...) 
#> {
#>     rgeom(1, p = param[1]) + 1
#> }
```

``` r

set.seed(1107)
res <- hsim(h, size = 1000)
#> Initial intensity components: auto -> zero (no supported unconditional mean formula for this model). Supplying lambda_component0 explicitly is recommended for this model.
summary(res)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>       arrival N1 mark lambda1
#>  [1,]  0.0000  0    0 0.15000
#>  [2,]  1.9559  1    1 0.15000
#>  [3,]  3.2143  2    1 0.24347
#>  [4,]  3.9757  3    1 0.38467
#>  [5,]  4.5737  4    1 0.50898
#>  [6,]  7.7746  5    1 0.15632
#>  [7,] 11.9296  6    2 0.15092
#>  [8,] 13.4276  7    1 0.24108
#>  [9,] 23.3489  8    1 0.15000
#> [10,] 31.9672  9    3 0.15000
#> [11,] 32.0264 10    2 1.33257
#> [12,] 32.1604 11    1 1.91124
#> [13,] 36.9857 12    1 0.15109
#> [14,] 37.1174 13    2 0.71795
#> [15,] 38.4267 14    1 0.34298
#> [16,] 38.9809 15    1 0.51792
#> [17,] 42.9833 16    2 0.15177
#> [18,] 48.4586 17    1 0.15016
#> [19,] 48.7744 18    1 0.57247
#> [20,] 48.7916 19    2 1.24194
#> ... with 980 more rows 
#> -------------------------------------------------------
```

``` r

fit <- hfit(h,
            inter_arrival = res$inter_arrival,
            type = res$type,
            mark = res$mark)
#> Initial intensity components: auto -> zero (no supported unconditional mean formula for this model). Supplying lambda_component0 explicitly is recommended for this model.
summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 35 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1684.224 
#> 4  free parameters
#> Estimates:
#>        Estimate Std. error t value  Pr(> t)    
#> mu1    0.145592   0.008433  17.264  < 2e-16 ***
#> alpha1 0.639384   0.066797   9.572  < 2e-16 ***
#> beta1  1.490181   0.116725  12.767  < 2e-16 ***
#> eta1   0.336223   0.077077   4.362 1.29e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

To estimate the parameters of the mark distribution as well, supply its
density or probability mass function in the `dmark` slot.

``` r

h_md <- h

h_md@dmark <- function(param = c(p = 0.1), n=n, mark=mark, ...){
   dgeom(mark[n] - 1, prob = param["p"])
}

mle_md <- hfit(h_md,
               inter_arrival = res$inter_arrival, type = res$type, mark = res$mark)
#> Initial intensity components: auto -> zero (no supported unconditional mean formula for this model). Supplying lambda_component0 explicitly is recommended for this model.
summary(mle_md)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 45 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -2657.065 
#> 5  free parameters
#> Estimates:
#>        Estimate Std. error t value  Pr(> t)    
#> mu1     0.14559    0.00844  17.251  < 2e-16 ***
#> alpha1  0.63939    0.06858   9.323  < 2e-16 ***
#> beta1   1.49020    0.11944  12.476  < 2e-16 ***
#> eta1    0.33622    0.07583   4.434 9.27e-06 ***
#> p       0.65897    0.01214  54.259  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

The density and generator use the same parameter name, `p`.
Reconstruction updates `rmark` with the value estimated through `dmark`,
so new simulations use the fitted mark distribution. Generator
parameters that do not match estimated coefficient names keep their
original defaults.

``` r

estimated_md <- fitted_hspec(mle_md)
c(density = eval(formals(estimated_md@dmark)$param)[["p"]],
  generator = eval(formals(estimated_md@rmark)$param)[["p"]])
#>   density generator 
#> 0.6589701 0.6589701
```
