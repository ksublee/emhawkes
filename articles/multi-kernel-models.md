# Multi-kernel models and shared parameters

These examples combine several exponential kernels and share parameters
across component functions. Each model defines its own specification and
simulated data.

See the [examples
overview](https://ksublee.github.io/emhawkes/articles/extended_example.md)
for other model topics.

``` r

library(emhawkes)
```

## Multi-kernel model

A multi-kernel Hawkes model requires `type_col_map` in its `hspec`
object. This list maps each event type to the corresponding columns of
the kernel parameter matrices. Consider a bivariate multi-kernel model:

``` math
\lambda_t = \mu + \int_{-\infty}^{t} h(t-u) d N(u)
```

where

``` math
h = \sum_{k=1}^{K} h_k, \quad
h_k (t) = \alpha_k \circ \begin{bmatrix}
e^{-\beta_{k11} t} & e^{-\beta_{k12} t} \\
e^{-\beta_{k21} t} & e^{-\beta_{k22} t}
\end{bmatrix}
```

Here, $`k`$ indexes the kernels, $`\alpha_k`$ is the excitation matrix
for kernel $`k`$, and $`\circ`$ denotes elementwise multiplication.

For a bivariate Hawkes model with two kernels, the intensity processes
are

``` math
 \begin{bmatrix} \lambda_1(t) \\ \lambda_2(t) \end{bmatrix} =
\begin{bmatrix} \mu_1 \\ \mu_2 \end{bmatrix} + \int_{-\infty}^{t} \begin{bmatrix}
\alpha_{111} e^{-\beta_{111}(t-s)} & \alpha_{112} e^{-\beta_{112}(t-s)} \\
\alpha_{121}e^{-\beta_{121}(t-s)} & \alpha_{122}e^{-\beta_{122}(t-s)}
\end{bmatrix} \begin{bmatrix} d N_1(s) \\ dN_2(s) \end{bmatrix}
+ \int_{-\infty}^{t} \begin{bmatrix}
\alpha_{211} e^{-\beta_{211}(t-s)} & \alpha_{212} e^{-\beta_{212}(t-s)} \\
\alpha_{221}e^{-\beta_{221}(t-s)} & \alpha_{222}e^{-\beta_{222}(t-s)}
\end{bmatrix} \begin{bmatrix} d N_1(s) \\ dN_2(s) \end{bmatrix}.
```

The parameter matrices are

``` math
 \alpha = \begin{bmatrix} \alpha_{111} & \alpha_{112} & \alpha_{211} & \alpha_{212} \\ \alpha_{121} & \alpha_{122} & \alpha_{221} & \alpha_{222}
\end{bmatrix}, \quad
 \beta = \begin{bmatrix} \beta_{111} & \beta_{112} & \beta_{211} & \beta_{212} \\ \beta_{121} & \beta_{122} & \beta_{221} & \beta_{222}
\end{bmatrix} \quad
```

The mapping must specify which columns correspond to each counting
process $`N_i`$.

``` r

mu <- matrix(c(0.02, 0.02), nrow=2)

beta_1 <- matrix(rep(10, 4), nrow=2)
beta_2 <- matrix(rep(1, 4), nrow=2)
beta  <- cbind(beta_1, beta_2)

alpha_1 <- matrix(c(3, 2,
                    2, 3), nrow=2, byrow=TRUE)
alpha_2 <- matrix(c(0.3, 0.2,
                    0.2, 0.3), nrow=2, byrow=TRUE)
alpha <- cbind(alpha_1, alpha_2)

print(alpha)
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    2  0.3  0.2
#> [2,]    2    3  0.2  0.3
```

The first and third columns of $`\alpha`$ correspond to $`N_1`$, and the
second and fourth columns correspond to $`N_2`$. Thus, define
`type_col_map` as follows:

``` r

type_col_map <- list(c(1,3),  # columns for dN1
                     c(2,4))  # columns for dN2
type_col_map
#> [[1]]
#> [1] 1 3
#> 
#> [[2]]
#> [1] 2 4
```

The entry `type_col_map[[i]]` lists the columns associated with event
type `i`. Rows of the intensity matrix identify event types. After a
type `i` event, `type_col_map[[i]]` selects the kernel columns to
update. Use this mapping to extract the corresponding parameters:

``` r

cat("Part of alpha associated with N1: \n")
#> Part of alpha associated with N1:
alpha[, type_col_map[[1]]]  # associated with N1
#>      [,1] [,2]
#> [1,]    3  0.3
#> [2,]    2  0.2
cat("Part of alpha associated with N2: \n")
#> Part of alpha associated with N2:
alpha[, type_col_map[[2]]]  # associated with N2
#>      [,1] [,2]
#> [1,]    2  0.2
#> [2,]    3  0.3

cat("Part of beta associated with N1: \n")
#> Part of beta associated with N1:
beta[, type_col_map[[1]]]  # associated with N1
#>      [,1] [,2]
#> [1,]   10    1
#> [2,]   10    1
cat("Part of beta associated with N2: \n")
#> Part of beta associated with N2:
beta[, type_col_map[[2]]]  # associated with N2
#>      [,1] [,2]
#> [1,]   10    1
#> [2,]   10    1
```

``` r

h <- new("hspec", mu = mu, alpha = alpha, beta=beta, type_col_map = type_col_map)
h
#> An object of class "hspec" of 2-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,] 0.02
#> [2,] 0.02
#> 
#> Slot alpha: 
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    2  0.3  0.2
#> [2,]    2    3  0.2  0.3
#> 
#> Slot beta: 
#>      [,1] [,2] [,3] [,4]
#> [1,]   10   10    1    1
#> [2,]   10   10    1    1
#> 
#> Slot type_col_map: 
#> [[1]]
#> [1] 1 3
#> 
#> [[2]]
#> [1] 2 4
```

For this multi-kernel model, `lambda_component0` must be supplied for
both simulation and estimation.

``` r

set.seed(620)
lambda_component0_mk <- matrix(seq(1, 1.7, 0.1), nrow = 2)
res_mk <- hsim(h, size = 3000, lambda_component0 = lambda_component0_mk)
res_mk
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> An object of class "hspec" of 2-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,] 0.02
#> [2,] 0.02
#> 
#> Slot alpha: 
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    2  0.3  0.2
#> [2,]    2    3  0.2  0.3
#> 
#> Slot beta: 
#>      [,1] [,2] [,3] [,4]
#> [1,]   10   10    1    1
#> [2,]   10   10    1    1
#> 
#> Slot type_col_map: 
#> [[1]]
#> [1] 1 3
#> 
#> [[2]]
#> [1] 2 4
#> 
#> 
#> Realized path :
#>       arrival N1 N2 mark lambda1 lambda2  lambda11  lambda12 lambda13 lambda14
#>  [1,]  0.0000  0  0    0   5.220   5.620 1.0000000 1.200e+00   1.4000   1.6000
#>  [2,]  0.3339  1  0    1   2.247   2.397 0.0354893 4.259e-02   1.0026   1.1459
#>  [3,]  1.2043  1  1    1   1.046   1.064 0.0005036 7.065e-06   0.5455   0.4799
#>  [4,]  1.2103  2  1    1   3.122   4.182 0.0004743 1.884e+00   0.5423   0.6758
#>  [5,]  1.3704  3  1    1   2.299   2.302 0.6051446 3.799e-01   0.7176   0.5758
#>  [6,]  1.4098  3  2    1   4.238   3.474 2.4298143 2.561e-01   0.9783   0.5535
#>  [7,]  1.5416  3  3    1   2.793   2.894 0.6507825 6.043e-01   0.8575   0.6605
#>  [8,]  1.6639  3  4    1   2.498   2.919 0.1914158 7.660e-01   0.7588   0.7614
#>  [9,]  2.0338  3  5    1   1.282   1.454 0.0047375 6.846e-02   0.5242   0.6642
#> [10,]  2.1098  4  5    1   2.276   2.981 0.0022151 9.671e-01   0.4858   0.8009
#> [11,]  2.3626  5  5    1   1.569   1.623 0.2396920 7.722e-02   0.6103   0.6220
#> [12,]  2.7052  5  6    1   1.216   1.178 0.1053725 2.511e-03   0.6462   0.4416
#> [13,]  2.7692  5  7    1   2.339   2.939 0.0555568 1.056e+00   0.6062   0.6018
#> [14,]  2.9333  6  7    1   1.818   2.272 0.0107638 5.920e-01   0.5144   0.6804
#> [15,]  3.0063  7  7    1   3.146   2.862 1.4508962 2.853e-01   0.7571   0.6325
#> [16,]  3.3222  7  8    1   1.453   1.365 0.1889842 1.211e-02   0.7707   0.4612
#> [17,]  3.3779  7  9    1   2.636   3.241 0.1082998 1.153e+00   0.7290   0.6254
#> [18,]  3.4084  8  9    1   3.931   5.226 0.0798097 2.324e+00   0.7071   0.8006
#> [19,]  3.7871  8 10    1   1.380   1.425 0.0698539 5.270e-02   0.6897   0.5482
#> [20,]  3.8329  9 10    1   2.736   3.505 0.0441702 1.298e+00   0.6588   0.7147
#>        lambda21  lambda22 lambda23 lambda24
#>  [1,] 1.1000000 1.300e+00   1.5000   1.7000
#>  [2,] 0.0390382 4.614e-02   1.0742   1.2175
#>  [3,] 0.0003383 7.654e-06   0.5336   0.5098
#>  [4,] 0.0003186 2.826e+00   0.5304   0.8050
#>  [5,] 0.4034302 5.699e-01   0.6224   0.6859
#>  [6,] 1.6198765 3.841e-01   0.7906   0.6594
#>  [7,] 0.4338551 9.064e-01   0.6930   0.8410
#>  [8,] 0.1276106 1.149e+00   0.6132   1.0095
#>  [9,] 0.0031583 1.027e-01   0.4236   0.9046
#> [10,] 0.0014767 1.451e+00   0.3926   1.1165
#> [11,] 0.1597947 1.158e-01   0.4602   0.8671
#> [12,] 0.0702484 3.767e-03   0.4687   0.6156
#> [13,] 0.0370378 1.584e+00   0.4397   0.8588
#> [14,] 0.0071759 8.881e-01   0.3731   0.9834
#> [15,] 0.9672641 4.280e-01   0.5328   0.9142
#> [16,] 0.1259895 1.817e-02   0.5343   0.6665
#> [17,] 0.0721999 1.730e+00   0.5053   0.9142
#> [18,] 0.0532064 3.485e+00   0.4901   1.1777
#> [19,] 0.0465693 7.905e-02   0.4726   0.8065
#> [20,] 0.0294468 1.947e+00   0.4514   1.0569
#> ... with 2980 more rows 
#> -------------------------------------------------------
```

``` r

fit_mk <- hfit(h, res_mk$inter_arrival, res_mk$type,
               lambda_component0 = lambda_component0_mk,
               final_hessian = "stable")
#> Registered S3 methods overwritten by 'maxLik':
#>   method         from    
#>   coef.maxLik    emhawkes
#>   logLik.maxLik  emhawkes
#>   print.maxLik   emhawkes
#>   summary.maxLik emhawkes
#>   vcov.maxLik    emhawkes
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, alpha, beta).
#> Use reduced = FALSE to estimate matrix entries separately.
summary(fit_mk)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 77 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: 6084.72 
#> 7  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1      0.011162   0.007743   1.442   0.1494    
#> alpha1.1 2.474140   0.350309   7.063 1.63e-12 ***
#> alpha1.2 1.380734   0.298057   4.632 3.61e-06 ***
#> alpha1.3 0.044026   0.097978   0.449   0.6532    
#> alpha1.4 0.231206   0.127207   1.818   0.0691 .  
#> beta1.1  6.202808   1.480638   4.189 2.80e-05 ***
#> beta1.3  0.700741   0.277825   2.522   0.0117 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

This example explicitly requests the checked final Hessian. The original
estimate is retained, but its score calls for tighter optimization and
multiple starting values. Kernel starting order and the fixed initial
state can affect the selected local solution.

``` r

fit_diagnostics(fit_mk)[c("status", "newton_decrement", "warnings")]
#> $status
#> [1] "ok"
#> 
#> $newton_decrement
#> [1] 0.0003106036
#> 
#> $warnings
#> character(0)
```

Reconstruction retains `type_col_map` and the supplied initial intensity
matrix. Each kernel uses its own decay rate in the likelihood and
residual integrals. The inferred likelihood agrees with the fitted
likelihood even though decay rates differ within each row.

``` r

estimated_mk <- fitted_hspec(fit_mk)
initial_mk <- fit_initial_state(fit_mk)
estimated_mk@type_col_map
#> [[1]]
#> [1] 1 3
#> 
#> [[2]]
#> [1] 2 4
initial_mk$lambda_component0
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0  1.2  1.4  1.6
#> [2,]  1.1  1.3  1.5  1.7
inferred_mk <- do.call(logLik, c(list(object = estimated_mk,
  inter_arrival = res_mk$inter_arrival, type = res_mk$type, infer = TRUE), initial_mk))
c(fitted = as.numeric(logLik(fit_mk)), restored = inferred_mk$loglikelihood)
#>   fitted restored 
#>  6084.72  6084.72
all.equal(as.numeric(logLik(fit_mk)), inferred_mk$loglikelihood)
#> [1] TRUE
```

## Synchronized intensity model

This is a two-kernel model with the following reparameterization.

``` math
 \mu = \begin{bmatrix} \theta/(1 - \kappa)/2 + \tilde\theta/(1 + \kappa)/2 \\
\theta/(1 - \kappa)/2 - \tilde\theta/(1 + \kappa)/2 \end{bmatrix}, \quad \theta = (\theta^- + \theta^+)/2,\quad \tilde\theta=(\theta^- -\theta^+)/2
```

``` math
\alpha = \begin{bmatrix} \zeta & \tilde\zeta & \zeta & -\tilde\zeta \\
\zeta & -\tilde\zeta & \zeta & \tilde\zeta
\end{bmatrix}, \quad \zeta = (\eta + \nu) / 2, \quad \tilde \zeta = (\eta - \nu)/ 2
```

``` math
\beta = \begin{bmatrix} \beta_1 & \beta_2 & \beta_1 & \beta_2 \\
\beta_1 & \beta_2 & \beta_1 & \beta_2 \end{bmatrix}, \quad \beta_1 = \beta(1 - \kappa), \quad \beta_2 = \beta(1 + \kappa)
```

To implement this reparameterization, define `mu`, `alpha`, and `beta`
as functions. Each function takes a named parameter vector, `param`, as
its first argument.

``` r

mu <- function(param = c(theta_p = 0.15, theta_n = 0.21, kappa = 0.12)){
  theta    <- (param["theta_n"] + param["theta_p"])/2
  theta_tl <- (param["theta_n"] - param["theta_p"])/2
  matrix(c(theta/2/(1 - param["kappa"]) + theta_tl/2/(1 + param["kappa"]),
           theta/2/(1 - param["kappa"]) - theta_tl/2/(1 + param["kappa"])), nrow=2)
}

alpha <- function(param = c(eta = 5, nu = 3)){
  zeta    <- (param["eta"] + param["nu"])/2
  zeta_tl <- (param["eta"] - param["nu"])/2
  matrix(c(zeta, zeta_tl, zeta, -zeta_tl,
           zeta, -zeta_tl, zeta, zeta_tl), nrow=2, byrow=TRUE)
}

beta <- function(param = c(beta = 12, kappa = 0.12)){
  beta1 <- param["beta"] * (1 - param["kappa"])
  beta2 <- param["beta"] * (1 + param["kappa"])
  matrix(c(beta1, beta2, beta1, beta2,
           beta1, beta2, beta1, beta2), nrow = 2, byrow = TRUE)
}

type_col_map <- list(c(1,2), c(3,4))

h_sy <- new("hspec", mu = mu, alpha = alpha, beta = beta, type_col_map = type_col_map)
h_sy
#> An object of class "hspec" of 2-dimensional Hawkes process
#> 
#> Slot mu: 
#> function (param = c(theta_p = 0.15, theta_n = 0.21, kappa = 0.12)) 
#> {
#>     theta <- (param["theta_n"] + param["theta_p"])/2
#>     theta_tl <- (param["theta_n"] - param["theta_p"])/2
#>     matrix(c(theta/2/(1 - param["kappa"]) + theta_tl/2/(1 + param["kappa"]), 
#>         theta/2/(1 - param["kappa"]) - theta_tl/2/(1 + param["kappa"])), 
#>         nrow = 2)
#> }
#> <bytecode: 0x55e6207ad810>
#> 
#> Slot alpha: 
#> function (param = c(eta = 5, nu = 3)) 
#> {
#>     zeta <- (param["eta"] + param["nu"])/2
#>     zeta_tl <- (param["eta"] - param["nu"])/2
#>     matrix(c(zeta, zeta_tl, zeta, -zeta_tl, zeta, -zeta_tl, zeta, 
#>         zeta_tl), nrow = 2, byrow = TRUE)
#> }
#> <bytecode: 0x55e621041740>
#> 
#> Slot beta: 
#> function (param = c(beta = 12, kappa = 0.12)) 
#> {
#>     beta1 <- param["beta"] * (1 - param["kappa"])
#>     beta2 <- param["beta"] * (1 + param["kappa"])
#>     matrix(c(beta1, beta2, beta1, beta2, beta1, beta2, beta1, 
#>         beta2), nrow = 2, byrow = TRUE)
#> }
#> <bytecode: 0x55e61ca33a80>
#> 
#> Slot type_col_map: 
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3 4
```

``` r

set.seed(1107)
# run simulation
res_sy <- hsim(h_sy, size = 2000, lambda_component0 = matrix(rep(1, 2 * 4), nrow=2))
summary(res_sy)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>       arrival N1 N2 mark lambda1 lambda2
#>  [1,]  0.0000  0  0    0 4.11567 4.08888
#>  [2,]  7.2506  0  1    1 0.11567 0.08888
#>  [3,]  7.4805  0  2    1 0.42304 0.48724
#>  [4,]  7.4915  1  2    1 3.09049 4.86829
#>  [5,] 17.7981  2  2    1 0.11567 0.08888
#>  [6,] 20.0276  2  3    1 0.11567 0.08888
#>  [7,] 20.0484  2  4    1 2.56966 4.05419
#>  [8,] 20.0609  3  4    1 4.94847 7.88899
#>  [9,] 20.1607  4  4    1 3.58697 3.81330
#> [10,] 20.2304  5  4    1 4.09539 3.38445
#> [11,] 20.5057  5  5    1 0.56602 0.47289
#> [12,] 26.4862  6  5    1 0.11567 0.08888
#> [13,] 32.5114  6  6    1 0.11567 0.08888
#> [14,] 32.6227  6  7    1 1.12663 1.54804
#> [15,] 32.6316  6  8    1 3.79509 5.94053
#> [16,] 32.6480  7  8    1 5.81283 9.13210
#> [17,] 32.6921  7  9    1 6.87877 7.59586
#> [18,] 32.7953  8  9    1 3.51752 4.17621
#> [19,] 32.9249  8 10    1 2.20298 1.94570
#> [20,] 34.7648  9 10    1 0.11567 0.08888
#> ... with 1980 more rows 
#> -------------------------------------------------------
```

The parameters to estimate and their initial values are taken from the
default `param` vectors. Parameters with the same name are shared across
functions. For example, `kappa` appears in both `mu` and `beta`, but
only one `kappa` is estimated.

``` r

fit_sy <- hfit(h_sy, inter_arrival=res_sy$inter_arrival,
               type=res_sy$type,
               lambda_component0 = matrix(rep(1, 2 * 4), nrow=2),
               final_hessian = "stable")
summary(fit_sy)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 81 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -260.6273 
#> 6  free parameters
#> Estimates:
#>         Estimate Std. error t value  Pr(> t)    
#> theta_p  0.08291    0.03956   2.096 0.036106 *  
#> theta_n  0.16264    0.03221   5.050 4.42e-07 ***
#> kappa    0.40990    0.15487   2.647 0.008127 ** 
#> eta      5.23254    0.44037  11.882  < 2e-16 ***
#> nu       2.89428    0.39985   7.238 4.54e-13 ***
#> beta    18.01832    4.73854   3.803 0.000143 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

The restored component functions use the estimated defaults. Both
functions retain the same estimated `kappa`. Their parameter-only
interface also allows automatic C++ likelihood and residual replay with
the saved initial state.

``` r

estimated_sy <- fitted_hspec(fit_sy)
c(mu = eval(formals(estimated_sy@mu)$param)[["kappa"]],
  beta = eval(formals(estimated_sy@beta)$param)[["kappa"]])
#>        mu      beta 
#> 0.4099007 0.4099007
```

``` r

replayed_sy <- do.call(logLik, c(list(object = estimated_sy,
  inter_arrival = res_sy$inter_arrival, type = res_sy$type, infer = TRUE),
  fit_initial_state(fit_sy)))
c(fitted = as.numeric(logLik(fit_sy)), replayed = replayed_sy$loglikelihood)
#>    fitted  replayed 
#> -260.6273 -260.6273
```
