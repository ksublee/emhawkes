# State-dependent impact models

These examples use `impact` to make excitation depend on observed event
counts or the current intensity components. Start with the
[state-dependent baseline
examples](https://ksublee.github.io/emhawkes/articles/state-dependent-baseline.md)
for simpler state feedback. Each model below defines its own
specification and simulated data.

See the [examples
overview](https://ksublee.github.io/emhawkes/articles/extended_example.md)
for other model topics.

``` r

library(emhawkes)
```

## Hawkes flocking model

Here, $`g`$ denotes the additional excitation supplied by `impact`. It
depends on the observed state and need not depend on the mark. In the
Hawkes flocking model, the kernel components are
``` math
\alpha = \begin{bmatrix}\alpha_{11} & \alpha_{12} & 0 & 0 \\\alpha_{12}& \alpha_{11} & 0 & 0 \\0 & 0 & \alpha_{33} & \alpha_{34} \\0 & 0 & \alpha_{34} & \alpha_{33} \end{bmatrix},
```
``` math
g = \begin{bmatrix} 0 & 0 & \alpha_{1w} 1_{\{C_1(t) < C_2(t)\}} & \alpha_{1n} 1_{\{C_1(t) < C_2(t)\}} \\
0 & 0 & \alpha_{1n} 1_{\{C_1(t) > C_2(t)\}} & \alpha_{1w}1_{\{C_1(t) > C_2(t)\}} \\
\alpha_{2w} 1_{\{C_2(t) < C_1(t)\}} & \alpha_{2n}1_{\{C_2(t) < C_1(t)\}} & 0 & 0 \\
\alpha_{2n} 1_{\{C_2(t) > C_1(t)\}} & \alpha_{2w}1_{\{C_2(t) > C_1(t)\}} & 0 & 0 \end{bmatrix},
```

where

``` math
C_1(t) = N_1(t) - N_2(t), \quad C_2(t) = N_3(t) - N_4(t).
```
For details, see [Systemic risk in market microstructure of crude oil
and gasoline futures prices: A Hawkes flocking model
approach](https://onlinelibrary.wiley.com/doi/full/10.1002/fut.22048).

The `alpha` slot can contain a function that returns a matrix. In the
following example, the function returns a $`4\times4`$ matrix with some
entries fixed at zero. Only the parameters in `param` are estimated. For
simulation, a fixed matrix and a function returning the same matrix give
the same model.

``` r

mu <- matrix(c(0.02, 0.02, 0.04, 0.04), nrow = 4)


alpha <- function(param = c(alpha11 = 0.2, alpha12 = 0.3, alpha33 = 0.3, alpha34 = 0.4)){
  matrix(c(param["alpha11"], param["alpha12"], 0, 0,
           param["alpha12"], param["alpha11"], 0, 0,
           0, 0, param["alpha33"], param["alpha34"],
           0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = TRUE)
}


beta <- matrix(c(rep(0.7, 8), rep(1.1, 8)), nrow = 4, byrow = TRUE)
```

The `impact` function below implements the model definition. Its `N`
argument contains the counting process values, and `n` identifies the
current event index. Include both arguments because the function uses
them, along with `...` to accept other arguments supplied by the
package. Use `N[n, "N1"]` to read one count directly, rather than
extracting an entire column at every event.

``` r

impact <- function(param = c(alpha1n=0.25, alpha1w=0.1, alpha2n=0.1, alpha2w=0.2),
                   N=N, n=n, ...){

  Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
                  0, 0, param['alpha1n'], param['alpha1w'],
                  param['alpha2w'], param['alpha2n'], 0, 0,
                  param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=TRUE)

  ind <- N[n,"N1"] - N[n,"N2"] > N[n,"N3"] - N[n,"N4"]

  km <- matrix(c(!ind, !ind, !ind, !ind,
                 ind, ind, ind, ind,
                 ind, ind, ind, ind,
                 !ind, !ind, !ind, !ind), nrow = 4, byrow = TRUE)

  km * Psi
}

hspec_fl <- new("hspec",
                mu = mu, alpha = alpha, beta = beta, impact = impact)
hspec_fl
#> An object of class "hspec" of 4-dimensional Hawkes process
#> 
#> Slot mu: 
#>      [,1]
#> [1,] 0.02
#> [2,] 0.02
#> [3,] 0.04
#> [4,] 0.04
#> 
#> Slot alpha: 
#> function (param = c(alpha11 = 0.2, alpha12 = 0.3, alpha33 = 0.3, 
#>     alpha34 = 0.4)) 
#> {
#>     matrix(c(param["alpha11"], param["alpha12"], 0, 0, param["alpha12"], 
#>         param["alpha11"], 0, 0, 0, 0, param["alpha33"], param["alpha34"], 
#>         0, 0, param["alpha34"], param["alpha33"]), nrow = 4, 
#>         byrow = TRUE)
#> }
#> <bytecode: 0x563856cf3448>
#> 
#> Slot beta: 
#>      [,1] [,2] [,3] [,4]
#> [1,]  0.7  0.7  0.7  0.7
#> [2,]  0.7  0.7  0.7  0.7
#> [3,]  1.1  1.1  1.1  1.1
#> [4,]  1.1  1.1  1.1  1.1
#> 
#> Slot impact: 
#> function (param = c(alpha1n = 0.25, alpha1w = 0.1, alpha2n = 0.1, 
#>     alpha2w = 0.2), N = N, n = n, ...) 
#> {
#>     Psi <- matrix(c(0, 0, param["alpha1w"], param["alpha1n"], 
#>         0, 0, param["alpha1n"], param["alpha1w"], param["alpha2w"], 
#>         param["alpha2n"], 0, 0, param["alpha2n"], param["alpha2w"], 
#>         0, 0), nrow = 4, byrow = TRUE)
#>     ind <- N[n, "N1"] - N[n, "N2"] > N[n, "N3"] - N[n, "N4"]
#>     km <- matrix(c(!ind, !ind, !ind, !ind, ind, ind, ind, ind, 
#>         ind, ind, ind, ind, !ind, !ind, !ind, !ind), nrow = 4, 
#>         byrow = TRUE)
#>     km * Psi
#> }
```

``` r

set.seed(1107)
hr_fl <- hsim(hspec_fl, size=2000)
#> Initial intensity components: auto -> zero (no supported unconditional mean formula for this model). Supplying lambda_component0 explicitly is recommended for this model.
summary(hr_fl)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>        arrival N1 N2 N3 N4 mark  lambda1 lambda2  lambda3  lambda4
#>  [1,]  0.00000  0  0  0  0    0 0.020000 0.02000 0.040000 0.040000
#>  [2,]  0.84055  0  0  1  0    1 0.020000 0.02000 0.040000 0.040000
#>  [3,]  7.35095  0  0  2  0    1 0.021049 0.02000 0.040233 0.040310
#>  [4,]  8.75286  1  0  2  0    1 0.057874 0.02000 0.104229 0.125638
#>  [5,] 29.54869  2  0  2  0    1 0.020000 0.02000 0.040000 0.040000
#>  [6,] 30.91017  2  1  2  0    1 0.097114 0.13567 0.040000 0.062366
#>  [7,] 31.08958  3  1  2  0    1 0.352605 0.29841 0.040000 0.222539
#>  [8,] 33.43266  4  1  2  0    1 0.123299 0.13218 0.040000 0.061465
#>  [9,] 33.76215  5  1  2  0    1 0.260826 0.34728 0.179195 0.054939
#> [10,] 34.28701  6  1  2  0    1 0.325287 0.45441 0.230420 0.048387
#> [11,] 35.07381  6  2  2  0    1 0.311303 0.44340 0.204307 0.043530
#> [12,] 35.31892  6  3  2  0    1 0.518076 0.54511 0.241844 0.042695
#> [13,] 35.92431  6  4  2  0    1 0.542398 0.49464 0.195086 0.041385
#> [14,] 36.30708  7  4  2  0    1 0.649096 0.53607 0.141792 0.172180
#> [15,] 36.96777  7  5  2  0    1 0.542096 0.53389 0.185908 0.103906
#> [16,] 37.00203  7  5  2  1    1 0.822615 0.71697 0.180512 0.294145
#> [17,] 37.02489  8  5  2  1    1 0.809878 0.80433 0.567094 0.580389
#> [18,] 37.58468  9  5  2  1    1 0.688962 0.75279 0.432793 0.331931
#> [19,] 39.13120  9  5  2  2    1 0.314339 0.36983 0.148166 0.093268
#> [20,] 40.16087 10  5  2  2    1 0.163160 0.23879 0.203723 0.153817
#> ... with 1980 more rows 
#> -------------------------------------------------------
```

``` r

fit_fl <- hfit(hspec_fl, hr_fl$inter_arrival, hr_fl$type)
#> Registered S3 methods overwritten by 'maxLik':
#>   method         from    
#>   coef.maxLik    emhawkes
#>   logLik.maxLik  emhawkes
#>   print.maxLik   emhawkes
#>   summary.maxLik emhawkes
#>   vcov.maxLik    emhawkes
#> Initial intensity components: auto -> zero (no supported unconditional mean formula for this model). Supplying lambda_component0 explicitly is recommended for this model.
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; mu, beta).
#> Use reduced = FALSE to estimate matrix entries separately.
summary(fit_fl)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 83 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -3050.707 
#> 12  free parameters
#> Estimates:
#>         Estimate Std. error t value  Pr(> t)    
#> mu1     0.028131   0.004309   6.528 6.65e-11 ***
#> mu3     0.033344   0.004265   7.818 5.39e-15 ***
#> alpha11 0.199823   0.028718   6.958 3.45e-12 ***
#> alpha12 0.312352   0.030265  10.320  < 2e-16 ***
#> alpha33 0.283701   0.036342   7.806 5.89e-15 ***
#> alpha34 0.341242   0.038529   8.857  < 2e-16 ***
#> beta1.1 0.734287   0.052146  14.081  < 2e-16 ***
#> beta3.1 0.902526   0.072865  12.386  < 2e-16 ***
#> alpha1n 0.249943   0.044270   5.646 1.64e-08 ***
#> alpha1w 0.107222   0.029765   3.602 0.000315 ***
#> alpha2n 0.084355   0.037596   2.244 0.024852 *  
#> alpha2w 0.173571   0.034026   5.101 3.38e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

This impact depends on observed event counts, so automatic selection
uses the C++ likelihood with numerical differentiation. No Jacobian is
needed; R still evaluates the impact.

``` r

fit_diagnostics(fit_fl)$computation$engine
#> $requested
#> [1] "auto"
#> 
#> $used
#> [1] "cpp"
#> 
#> $fallback_reason
#> NULL
```

The same observed impact also supports automatic C++ intensity and
residual replay. No Jacobian is needed; retain the fit’s initial-state
arguments.

``` r

replayed_fl <- do.call(logLik, c(list(object = fitted_hspec(fit_fl),
  inter_arrival = hr_fl$inter_arrival, type = hr_fl$type, infer = TRUE),
  fit_initial_state(fit_fl)))
c(fitted = as.numeric(logLik(fit_fl)), replayed = replayed_fl$loglikelihood)
#>    fitted  replayed 
#> -3050.707 -3050.707
```

## Bid-ask price model

This model describes bid and ask prices using a system of counting
processes and their conditional intensities:

``` math
N_t = \begin{bmatrix}
N_1(t) \\ N_2(t) \\ N_3(t) \\ N_4(t)
\end{bmatrix},
\quad
\lambda_t =
\begin{bmatrix} \lambda_1(t) \\ \lambda_2(t) \\ \lambda_3(t) \\ \lambda_4(t)
\end{bmatrix}
```

In the unmarked case, the changes in ask and bid prices are
$`N_1(t) - N_2(t)`$ and $`N_3(t) - N_4(t)`$, respectively. The marked
example below uses `Nc` to accumulate price changes. The midprice
$`p(t)`$ is the average of the ask and bid prices, including their
initial levels.

The base intensity process is

``` math
\mu = \begin{bmatrix} \mu_1 \\ \zeta \ell(t-) \\ \zeta \ell(t-) \\ \mu_1 \end{bmatrix}, \quad \ell(t) = \frac{L(t)}{p(t)} 
```

where $`L(t) \in \{0, 1, 2, \ldots\}`$ measures the bid-ask spread above
its minimum level, so $`L(t)=0`$ corresponds to the minimum spread. For
details, see [Modeling bid and ask price dynamics with an extended
Hawkes process and its empirical applications for high-frequency stock
market
data](https://academic.oup.com/jfec/article/21/4/1099/6513935?login=true).

In the `mu` function below, `n` identifies the current event index, and
`Nc` is used to calculate the spread level and midprice.

``` r

# presumed initial bid and ask prices
initial_ask_price <- 1250 #cents
initial_bid_price <- 1150 #cents

initial_level <- round((initial_ask_price - initial_bid_price) - 1)
initial_mid_price <- (initial_bid_price + initial_ask_price) / 2

mu <- function(param = c(mu1 = 0.08, zeta1 = 0.10), n=n, Nc=Nc, ...){

  if(n == 1){

    level <- initial_level
    mid <- initial_mid_price

  } else {

    level <- Nc[n-1,1] - Nc[n-1,2] - (Nc[n-1,3] - Nc[n-1,4]) + initial_level
    ask <- initial_ask_price + (Nc[n-1,1] - Nc[n-1,2])
    bid <- initial_bid_price + (Nc[n-1,3] - Nc[n-1,4])
    mid <- (ask + bid) / 2

  }

  if(level <= 0){
    matrix(c(param["mu1"], 0,
             0, param["mu1"]), nrow = 4)
  } else {
    matrix(c(param["mu1"], param["zeta1"] * level / mid,
             param["zeta1"]*level / mid, param["mu1"]), nrow = 4)

  }

}
```

In addition, the kernel is represented by

``` math
h(t, u) =
\begin{bmatrix} \alpha_{s1} & \alpha_{m} & \alpha_{s2} & 0 \\   \alpha_{w1} & \alpha_{n1}(u) & \alpha_{n1}(u) & \alpha_{w2} \\  \alpha_{w2} & \alpha_{n2}(u) & \alpha_{n2}(u) & \alpha_{w1} \\  0 & \alpha_{s2} & \alpha_{m} & \alpha_{s1} \\   \end{bmatrix},
```

where

``` math
\alpha_{n1}(u) = - \sum_{j=1}^4 \lambda_{2j}(u) + \xi \ell(u), \quad \alpha_{n2}(u) = - \sum_{j=1}^4 \lambda_{3j}(u) + \xi \ell(u),
```

Here, $`\xi \geq 0`$ is a constant, and $`\lambda_{ij}`$ is a component
of $`\lambda_i`$, defined by

``` math
\lambda_{ij}(t) = \int_{-\infty}^t h_{ij}(t, u) d N_j(u).
```

The following code places the constant part of $`h`$ in `alpha` and the
state-dependent part in `impact`. The `lambda_component_n` argument
supplies the current pre-jump matrix of $`\lambda_{ij}`$ values. This
avoids allocating a full intensity history when only the current state
is needed. Individual components can be negative after a reset. The
equal decay rates used here allow simulation as long as the sum of the
components in each row remains nonnegative.

``` r


alpha <- function(param = c(alpha_s1=4, alpha_m=26, alpha_s2=5,
                            alpha_w1=11, alpha_w2=7)){
  matrix(c(param["alpha_s1"], param["alpha_m"], param["alpha_s2"], 0,
           param["alpha_w1"], 0, 0, param["alpha_w2"],
           param["alpha_w2"], 0, 0, param["alpha_w1"],
           0, param["alpha_s2"], param["alpha_m"], param["alpha_s1"]), nrow = 4, byrow = TRUE)
}

impact <- function(param = c(xi = 2.7), n=n, Nc=Nc, lambda_component_n = lambda_component_n, ... ){
  if(n == 1){
    level <-  initial_level
    # mid <- initial_mid_price
  } else {
    level <- Nc[n,1] - Nc[n,2] - (Nc[n,3] - Nc[n,4]) + initial_level
    ask <- initial_ask_price + (Nc[n,1] - Nc[n,2])
    bid <- initial_bid_price + (Nc[n,3] - Nc[n,4])
    mid <- (ask + bid)/2
  }

  lambda_component_matrix <- lambda_component_n

  l2 <- sum(lambda_component_matrix[2,]) # sum of second row
  l3 <- sum(lambda_component_matrix[3,]) # sum of third row

  im <- matrix(c(0, 0, 0, 0,
                 0, -l2 + param["xi"]*level/mid, -l2 + param["xi"]*level/mid, 0,
                 0, -l3 + param["xi"]*level/mid, -l3 + param["xi"]*level/mid, 0,
                 0, 0, 0, 0), nrow = 4, byrow = TRUE)

}

beta <- matrix(rep(50, 16), nrow = 4, byrow=TRUE)

rmark <- function(n=n, Nc=Nc, type, ...){
  if(n == 1){
    level <-  initial_level
  } else {
    level <- Nc[n-1,1] - Nc[n-1,2] - (Nc[n-1,3] - Nc[n-1,4]) + initial_level
  }
  if (type[n] == 2 | type[n] == 3){
    min(level,  rgeom(1, p=0.65) + 1)
  } else {
    rgeom(1, p=0.65) + 1
  }
}


h_ba <- new("hspec", mu = mu, alpha = alpha, beta = beta, impact=impact, rmark = rmark)
h_ba
#> An object of class "hspec" of 4-dimensional Hawkes process
#> 
#> Slot mu: 
#> function (param = c(mu1 = 0.08, zeta1 = 0.1), n = n, Nc = Nc, 
#>     ...) 
#> {
#>     if (n == 1) {
#>         level <- initial_level
#>         mid <- initial_mid_price
#>     }
#>     else {
#>         level <- Nc[n - 1, 1] - Nc[n - 1, 2] - (Nc[n - 1, 3] - 
#>             Nc[n - 1, 4]) + initial_level
#>         ask <- initial_ask_price + (Nc[n - 1, 1] - Nc[n - 1, 
#>             2])
#>         bid <- initial_bid_price + (Nc[n - 1, 3] - Nc[n - 1, 
#>             4])
#>         mid <- (ask + bid)/2
#>     }
#>     if (level <= 0) {
#>         matrix(c(param["mu1"], 0, 0, param["mu1"]), nrow = 4)
#>     }
#>     else {
#>         matrix(c(param["mu1"], param["zeta1"] * level/mid, param["zeta1"] * 
#>             level/mid, param["mu1"]), nrow = 4)
#>     }
#> }
#> 
#> Slot alpha: 
#> function (param = c(alpha_s1 = 4, alpha_m = 26, alpha_s2 = 5, 
#>     alpha_w1 = 11, alpha_w2 = 7)) 
#> {
#>     matrix(c(param["alpha_s1"], param["alpha_m"], param["alpha_s2"], 
#>         0, param["alpha_w1"], 0, 0, param["alpha_w2"], param["alpha_w2"], 
#>         0, 0, param["alpha_w1"], 0, param["alpha_s2"], param["alpha_m"], 
#>         param["alpha_s1"]), nrow = 4, byrow = TRUE)
#> }
#> <bytecode: 0x563859928b40>
#> 
#> Slot beta: 
#>      [,1] [,2] [,3] [,4]
#> [1,]   50   50   50   50
#> [2,]   50   50   50   50
#> [3,]   50   50   50   50
#> [4,]   50   50   50   50
#> 
#> Slot impact: 
#> function (param = c(xi = 2.7), n = n, Nc = Nc, lambda_component_n = lambda_component_n, 
#>     ...) 
#> {
#>     if (n == 1) {
#>         level <- initial_level
#>     }
#>     else {
#>         level <- Nc[n, 1] - Nc[n, 2] - (Nc[n, 3] - Nc[n, 4]) + 
#>             initial_level
#>         ask <- initial_ask_price + (Nc[n, 1] - Nc[n, 2])
#>         bid <- initial_bid_price + (Nc[n, 3] - Nc[n, 4])
#>         mid <- (ask + bid)/2
#>     }
#>     lambda_component_matrix <- lambda_component_n
#>     l2 <- sum(lambda_component_matrix[2, ])
#>     l3 <- sum(lambda_component_matrix[3, ])
#>     im <- matrix(c(0, 0, 0, 0, 0, -l2 + param["xi"] * level/mid, 
#>         -l2 + param["xi"] * level/mid, 0, 0, -l3 + param["xi"] * 
#>             level/mid, -l3 + param["xi"] * level/mid, 0, 0, 0, 
#>         0, 0), nrow = 4, byrow = TRUE)
#> }
#> 
#> Slot rmark: 
#> function (n = n, Nc = Nc, type, ...) 
#> {
#>     if (n == 1) {
#>         level <- initial_level
#>     }
#>     else {
#>         level <- Nc[n - 1, 1] - Nc[n - 1, 2] - (Nc[n - 1, 3] - 
#>             Nc[n - 1, 4]) + initial_level
#>     }
#>     if (type[n] == 2 | type[n] == 3) {
#>         min(level, rgeom(1, p = 0.65) + 1)
#>     }
#>     else {
#>         rgeom(1, p = 0.65) + 1
#>     }
#> }
```

``` r

set.seed(1107)
hr_ba <- hsim(h_ba, size=1000, lambda_component0 = matrix(rep(1, 16), 4))
summary(hr_ba)
#> -------------------------------------------------------
#> Simulation result of exponential (marked) Hawkes model.
#> Realized path :
#>       arrival N1 N2 N3 N4 mark lambda1    lambda2   lambda3 lambda4
#>  [1,]  0.0000  0  0  0  0    0 4.08000  4.0082500 4.0082500  4.0800
#>  [2,]  3.6673  1  0  0  0    2 0.08000  0.0082500 0.0082500  0.0800
#>  [3,]  6.9225  2  0  0  0    1 0.08000  0.0084097 0.0084097  0.0800
#>  [4,] 20.0251  2  0  0  1    1 0.08000  0.0084894 0.0084894  0.0800
#>  [5,] 24.3957  2  0  0  2    3 0.08000  0.0085762 0.0085762  0.0800
#>  [6,] 24.4026  2  0  1  2    1 0.08000  4.9856950 7.8296138  2.9239
#>  [7,] 27.9067  3  0  1  2    1 0.08000  0.0087500 0.0087500  0.0800
#>  [8,] 27.9378  3  1  1  2    1 0.92442  2.3309869 1.4865661  0.0800
#>  [9,] 33.8868  3  1  2  2    2 0.08000  0.0087500 0.0087500  0.0800
#> [10,] 33.8900  3  1  2  3    1 4.33027  0.2054120 0.2054120 22.1814
#> [11,] 34.3764  4  1  2  3    1 0.08000  0.0086631 0.0086631  0.0800
#> [12,] 34.4342  4  2  2  3    1 0.30328  0.6227650 0.3994842  0.0800
#> [13,] 36.5641  5  2  2  3    2 0.08000  0.0086631 0.0086631  0.0800
#> [14,] 36.5658  5  3  2  3    2 3.74917 10.0990280 6.4298623  0.0800
#> [15,] 40.4590  5  3  2  4    1 0.08000  0.0086631 0.0086631  0.0800
#> [16,] 40.4666  5  3  3  4    1 0.08000  4.8059475 7.5472032  2.8213
#> [17,] 43.8835  6  3  3  4    5 0.08000  0.0086631 0.0086631  0.0800
#> [18,] 44.8818  7  3  3  4    2 0.08000  0.0090607 0.0090607  0.0800
#> [19,] 47.6976  7  3  3  5    2 0.08000  0.0092193 0.0092193  0.0800
#> [20,] 51.7138  8  3  3  5    1 0.08000  0.0093932 0.0093932  0.0800
#> ... with 980 more rows 
#> -------------------------------------------------------
```

The following likelihood calculation and fit do not estimate the mark
distribution because no `dmark` function is supplied.

``` r

logLik(h_ba, inter_arrival = hr_ba$inter_arrival, type = hr_ba$type,
       mark = hr_ba$mark, Nc = hr_ba$Nc,
       lambda_component0 = matrix(rep(1, 16), 4))
#> loglikelihood 
#>     -1817.794
```

``` r

mle_ba <- hfit(h_ba, inter_arrival = hr_ba$inter_arrival, type = hr_ba$type,
               mark = hr_ba$mark,
               lambda_component0 = matrix(rep(1, 16), 4),
               final_hessian = "stable")
#> hfit: equal starting values within each matrix share a parameter (reduced = TRUE; beta).
#> Use reduced = FALSE to estimate matrix entries separately.
summary(mle_ba)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 31 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1810.892 
#> 9  free parameters
#> Estimates:
#>           Estimate Std. error t value  Pr(> t)    
#> mu1       0.069541   0.003566  19.501  < 2e-16 ***
#> zeta1     0.101083   0.007967  12.687  < 2e-16 ***
#> alpha_s1  4.656214   0.789575   5.897 3.70e-09 ***
#> alpha_m  23.561095   2.345578  10.045  < 2e-16 ***
#> alpha_s2  5.850205   1.063972   5.498 3.83e-08 ***
#> alpha_w1 10.286378   1.135937   9.055  < 2e-16 ***
#> alpha_w2  7.403771   0.945758   7.828 4.94e-15 ***
#> beta1.1  51.872684   2.874442  18.046  < 2e-16 ***
#> xi        1.764779   0.838421   2.105   0.0353 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

Automatic selection uses the C++ likelihood for this observed baseline
and sequential impact with numerical differentiation. The impact reads
the current state before its jump. No Jacobian is required. Inspect the
saved engine record to confirm which calculation was used and whether
automatic selection fell back to R.

This example needs only the current component matrix, so it uses
`lambda_component_n`. Requesting `lambda_component` also stores the full
history. If a callback retains histories beyond its call, preserving
those snapshots can copy the history repeatedly. Cumulative allocation
can then grow quadratically with event count. This cost is repeated for
each likelihood evaluation during numerical optimization. Peak live
memory can be much smaller because old copies are reclaimed.

``` r

fit_diagnostics(mle_ba)$computation$engine
#> $requested
#> [1] "auto"
#> 
#> $used
#> [1] "cpp"
#> 
#> $fallback_reason
#> NULL
```

Replay the fitted model with the same initialization and observations.
Automatic selection also uses C++ for its likelihood, intensity paths,
and residuals.

``` r

replayed_ba <- do.call(logLik, c(list(object = fitted_hspec(mle_ba),
  inter_arrival = hr_ba$inter_arrival, type = hr_ba$type, mark = hr_ba$mark,
  infer = TRUE), fit_initial_state(mle_ba)))
stopifnot(isTRUE(all.equal(unname(replayed_ba$loglikelihood),
  as.numeric(mle_ba$maximum), tolerance = 1e-8)))
```

For the optional analytic fit, see the [Bid-ask Jacobians and analytic
gradient](https://ksublee.github.io/emhawkes/articles/computation.html#bid-ask-analytic-gradient).
That section continues from this example and gives the complete
derivatives and fit. Numerical differentiation remains the default.
