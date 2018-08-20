Introduction to emhawkes
================

### Univariate Hawkes process

``` r
library(emhawkes)
```

This subsection explaines how to construct, simulate, and estimate a univariate Hawkes model. First, create a `hspec` which defines the Hawkes model. S4 class `hspec` contains slots of model parameters, `mu`, `alpha`, `beta`, `dimens`, `rmark`, and `impact`.

The parameters of the model, `mu`, `alpha`, `beta`, are defined by matrices. The following is an example of a univariate model (without mark).

``` r
set.seed(1107)
mu1 <- 0.3; alpha1 <- 1.2; beta1 <- 1.5
hspec1 <- new("hspec", mu=mu1, alpha=alpha1, beta=beta1)
show(hspec1)
#> An object of class "hspec"
#> Slot "mu":
#>      [,1]
#> [1,]  0.3
#> 
#> Slot "alpha":
#>      [,1]
#> [1,]  1.2
#> 
#> Slot "beta":
#>      [,1]
#> [1,]  1.5
#> 
#> Slot "dimens":
#> [1] 1
#> 
#> Slot "rmark":
#> function(...) 1
#> <environment: 0x000000001a4fccf0>
#> 
#> Slot "impact":
#> NULL
```

The function `hsim` implements simulation where the input arguments are `hspec`, `size` and the initial values of intensity `lambda` and Hawkes processes `N`. If the initial values of `lambda` are omitted, the internally determined initial values are used. The default initial value of `N` is zero.

``` r
res1 <- hsim(hspec1, size=100)
summary(res1)
#> ------------------------------------`------
#> Simulation result of marked Hawkes model.
#> Realized path (with right continuous representation):
#>        arrival N1 lambda1
#>  [1,]  0.00000  0 0.90000
#>  [2,]  0.97794  1 0.43838
#>  [3,]  1.09001  2 1.43128
#>  [4,]  1.28999  3 2.02711
#>  [5,]  1.53225  4 2.33527
#>  [6,]  1.65001  5 3.01139
#>  [7,]  2.51807  6 1.36377
#>  [8,]  2.81710  7 1.74553
#>  [9,]  2.87547  8 2.72378
#> [10,]  3.16415  9 2.65016
#> [11,]  3.51378 10 2.40131
#> [12,]  4.22355 11 1.43843
#> [13,] 16.96752 12 0.30000
#> [14,] 17.71654 13 0.69015
#> [15,] 19.10293 14 0.49874
#> [16,] 24.06354 15 0.30082
#> [17,] 24.09256 16 1.44967
#> [18,] 28.40173 17 0.30366
#> [19,] 28.53743 18 1.28198
#> [20,] 28.56702 19 2.38725
#> ... with 80 more rows 
#> ------------------------------------------
```

The results of `hsim` is an S3 class `hreal` which consists of `hspec`, `inter_arrival`, `arrival`, `type`, `mark`, `N`, `Nc`, `lambda`, `lambda_component`, `rambda`, `rambda\_component`. The components `inter_arrival`, `type`, `mark`, `N`, `Nc`, `lambda`, `lambda_component`, `rambda`, `rambda_component` can be excessed during simulation and estimation. The column names of `N` are `N1`, `N2`, `N3` and so on. Similarly, the column names of `lambda` are `lambda1`, `lambda2`, `lambda3` and so on. The column names of `lambda_component` are `lambda_component11`, `lambda_component12`, `lambda_component13` and so on. `inter_arrival`, `type`, `mark`, `N`, `Nc` start at zero. Using `summary` function, one can print the first 20 elements of `arrival`, `N` and `lambda`.

`hspec` is the model specification, `inter_arrival` is the inter-arrival time of every event, and `arrival` is the cumulative sum of `inter_arrival`. `type` is the type of events, i.e., *i* for *N*<sub>*i*</sub>, and `mark` is a numeric vector which represents additional information for the event. `lambda` represents *λ* which is the left continuous and right limit version. The right continuous version of intensity is `rambda`. `lambda_component` represents *λ*<sub>*i**j*</sub> and `rambda_component` is the right continuous version.

The log-likelihood function is computed by `logLik` method. In this case, the inter-arrival times and `hspec` are inputs of the function.

``` r
logLik(hspec1, inter_arrival = res1$inter_arrival)
#> [1] -28.09111
```

The likelihood estimation is performed using `mhfit` function. The specification of the initial values of the parameters, `mhspec0` is needed. In the following example, `hspec0` is set to be `hspec1`, which is defined previously, for simplicity, but any candidates for the starting value of the numerical procedure can be used.

Note that only `inter_arrival` is needed. (Indeed, for more precise simulation, `lambda0`, the inital value of lambda compoment, should be specified. If not, internally determined initial values are set.)

``` r
mu0 <- 0.5; alpha0 <- 1.0; beta0 <- 1.8
hspec0 <- new("hspec", mu=mu0, alpha=alpha0, beta=beta0)
mle <- hfit(hspec0, inter_arrival = res1$inter_arrival)
#>    mu1 alpha1  beta1 
#>    0.5    1.0    1.8
summary(mle)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 44 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -26.68325 
#> 3  free parameters
#> Estimates:
#>        Estimate Std. error t value  Pr(> t)    
#> mu1     0.29680    0.09286   3.196  0.00139 ** 
#> alpha1  1.75745    0.40859   4.301 1.70e-05 ***
#> beta1   2.14300    0.49441   4.334 1.46e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

### Bivariate Hawkes model

In a bivariate model, the parameters, the slots of `hspec`, are matrices. `mu` is 2-by-1, and `alpha` and `beta` are 2-by-2 matrices. `rmark` is a random number generating function. `lambda0`, 2-by-2 matrix, represents the initial values of `lambda_component`, a set of `lambda11, lambda12, lambda21, lambda22`. The intensity processes are represented by

*λ*<sub>1</sub>(*t*)=*μ*<sub>1</sub> + *λ*<sub>11</sub>(*t*)+*λ*<sub>12</sub>(*t*)

*λ*<sub>2</sub>(*t*)=*μ*<sub>2</sub> + *λ*<sub>21</sub>(*t*)+*λ*<sub>22</sub>(*t*)

*λ*<sub>*i**j*</sub> called lambda components and `lambda0` is the time zero values of *λ*<sub>*i**j*</sub>, i.e., *λ*<sub>*i**j*</sub>(0). `lambda0` can be omitted and then internally determined initial values are used.

``` r
mu2 <- matrix(c(0.2), nrow = 2)
alpha2 <- matrix(c(0.5, 0.9, 0.9, 0.5), nrow = 2, byrow = TRUE)
beta2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow = TRUE)
lambda0 <- matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2, byrow = TRUE)
hspec2 <- new("hspec", mu=mu2, alpha=alpha2, beta=beta2)
print(hspec2)
#> An object of class "hspec"
#> Slot "mu":
#>      [,1]
#> [1,]  0.2
#> [2,]  0.2
#> 
#> Slot "alpha":
#>      [,1] [,2]
#> [1,]  0.5  0.9
#> [2,]  0.9  0.5
#> 
#> Slot "beta":
#>      [,1] [,2]
#> [1,] 2.25 2.25
#> [2,] 2.25 2.25
#> 
#> Slot "dimens":
#> [1] 2
#> 
#> Slot "rmark":
#> function(...) 1
#> <bytecode: 0x000000001a126840>
#> <environment: 0x000000001d55f608>
#> 
#> Slot "impact":
#> NULL
```

To simulate, use function `mhsim`.

``` r
res2 <- hsim(hspec2,  size=1000)
#> Warning in hsim(hspec2, size = 1000): The initial values for intensity processes are not provided. Internally determined initial values are used.
summary(res2)
#> ------------------------------------`------
#> Simulation result of marked Hawkes model.
#> Realized path (with right continuous representation):
#>       arrival N1 N2 lambda1 lambda2
#>  [1,]  0.0000  0  0 0.52941 0.52941
#>  [2,]  2.4745  1  0 0.20126 0.20126
#>  [3,]  2.5307  2  0 0.64180 0.99436
#>  [4,]  2.9229  2  1 0.58968 0.90105
#>  [5,]  5.8320  3  1 0.20185 0.20173
#>  [6,]  7.3825  4  1 0.21533 0.22754
#>  [7,]  7.6700  4  2 0.46986 0.68572
#>  [8,]  8.2456  4  3 0.52043 0.46999
#>  [9,]  8.3886  4  4 1.08460 0.75811
#> [10,]  8.5108  5  4 1.55553 1.00371
#> [11,]  8.6228  6  4 1.64232 1.52431
#> [12,]  8.8592  6  5 1.34104 1.50670
#> [13,]  8.9916  7  5 1.71515 1.54119
#> [14,]  9.0341  8  5 2.03128 2.23670
#> [15,]  9.0367  8  6 2.51809 3.12008
#> [16,]  9.0584  9  6 3.26436 3.45670
#> [17,]  9.6773 10  6 1.08553 1.23269
#> [18,]  9.7228 10  7 1.45073 1.94465
#> [19,]  9.7569 11  7 2.19206 2.27905
#> [20,]  9.8171 11  8 2.37611 2.80136
#> ... with 980 more rows 
#> ------------------------------------------
```

Frome the result, we get a vector of realized `inter_arrival` and `type`. Bivariate model requires `inter_arrival` and `type` for estimation.

``` r
inter_arrival2 <- res2$inter_arrival
type2 <- res2$type
```

Log-likelihood is computed by a function `logLik`.

``` r
logLik(hspec2, inter_arrival = inter_arrival2, type = type2)
#> Warning in .local(object, ...): The initial values for intensity processes are not provided. Internally determined initial values are used.
#> [1] -1063.402
```

A log-likelihood estimation is performed using `hfit`. In the following, the values of parameter slots in `mhspec0`, such as `mu, alpha, beta`, are regarded as a starting point of the numerical optimization. For simplicity, we use `hspec0 <- hspec2`. Since the true parameter values are not known in the actual problem, the initial value should be guessed. The realized `inter_arrival` and `type` are used.

``` r
hspec0 <- hspec2
mle <- hfit(hspec0, inter_arrival = inter_arrival2, type = type2)
#>      mu1 alpha1,1 alpha1,2  beta1,1 
#>     0.20     0.50     0.90     2.25
summary(mle)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 31 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1061.707 
#> 4  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.20316    0.01580  12.854  < 2e-16 ***
#> alpha1,1  0.59956    0.08137   7.369 1.72e-13 ***
#> alpha1,2  0.89903    0.09506   9.457  < 2e-16 ***
#> beta1,1   2.26760    0.20169  11.243  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

### Parameter setting

This subsection explains about the relation between parameter setting and estimation procedure in two-dimensional Hawkes model. The number of parameters to be estimated in the model depends on how we set the parameter slots such as `alpha` and `beta` in `mhspec0`, the sepcification for initial values.. Since the paremeter slot such as `alpha` is a matrix, and the element in the matrix can be the same or different. The number of parameters in the estimation varies depending on whether or not some of the elements in the initial setting are the same or different.

For example, if `alpha[1,1]` and `alpha[1,2]` in `mhspec0` are different, the numerical procedure tries to estimate both parameters of `alpha[1,1]` and `alpha[1,2]`. If `alpha[1,1]` and `alpha[1,2]` are the same in the initial setting, then the estimation procedure considered two parameters are the same in the model and hence only one of them is estimated.

Recall that th example in the previous section is of a symmetric Hawkes model where `alpha` is symmetric.

``` r
print(hspec2)
#> An object of class "hspec"
#> Slot "mu":
#>      [,1]
#> [1,]  0.2
#> [2,]  0.2
#> 
#> Slot "alpha":
#>      [,1] [,2]
#> [1,]  0.5  0.9
#> [2,]  0.9  0.5
#> 
#> Slot "beta":
#>      [,1] [,2]
#> [1,] 2.25 2.25
#> [2,] 2.25 2.25
#> 
#> Slot "dimens":
#> [1] 2
#> 
#> Slot "rmark":
#> function(...) 1
#> <bytecode: 0x000000001a126840>
#> <environment: 0x000000001d55f608>
#> 
#> Slot "impact":
#> NULL
```

``` r
res2 <- hsim(hspec2,  size=1000)
```

In the first example of estimation, `alpha0` is a matrix where the all elements have the same value, 0.75. In this setting, `mhfit` considers that `alpha11 == alpha12 == alpha21 == alpha22` in the model (even though the actual parameters have different values). Similarly for other parmater matrix `mu0` and `beta0`.

``` r
mu0 <- matrix(c(0.15, 0.15), nrow = 2)
alpha0 <- matrix(c(0.75, 0.75, 0.75, 0.75), nrow = 2, byrow=TRUE)
beta0 <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow = 2, byrow=TRUE)

hspec0 <- new("hspec", mu=mu0, alpha=alpha0, beta=beta0)
summary(hfit(hspec0, inter_arrival = res2$inter_arrival, type = res2$type))
#>      mu1 alpha1,1  beta1,1 
#>     0.15     0.75     2.60
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 28 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1324.235 
#> 3  free parameters
#> Estimates:
#>          Estimate Std. error t value Pr(> t)    
#> mu1       0.18792    0.01335   14.08  <2e-16 ***
#> alpha1,1  0.72289    0.06156   11.74  <2e-16 ***
#> beta1,1   2.34526    0.19789   11.85  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

In the second example, `alpha0`'s elements are not same, but symmetric as in the original simulation. We have `alpha11 == alpha22` and `alpha11 == alpha22` in `alpha0` and hence `alpha11` and `alpha12` will be estimated.

``` r
mu0 <- matrix(c(0.15, 0.15), nrow = 2)
alpha0 <- matrix(c(0.75, 0.751, 0.751, 0.75), nrow = 2, byrow=TRUE)
beta0 <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow = 2, byrow=TRUE)

hspec0 <- new("hspec", mu=mu0, alpha=alpha0, beta=beta0)
summary(hfit(hspec0, inter_arrival = res2$inter_arrival, type = res2$type))
#>      mu1 alpha1,1 alpha1,2  beta1,1 
#>    0.150    0.750    0.751    2.600
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 32 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1314.076 
#> 4  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.18883    0.01327  14.233  < 2e-16 ***
#> alpha1,1  0.48009    0.06788   7.073 1.52e-12 ***
#> alpha1,2  0.96746    0.09065  10.673  < 2e-16 ***
#> beta1,1   2.35526    0.18428  12.781  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

By simply setting `reduced = FALSE`, all parameters are estimated (not recommended).

``` r
summary(hfit(hspec2, inter_arrival = res2$inter_arrival, type = res2$type, reduced=FALSE))
#>      mu1      mu2 alpha1,1 alpha2,1 alpha1,2 alpha2,2  beta1,1  beta2,1 
#>     0.20     0.20     0.50     0.90     0.90     0.50     2.25     2.25 
#>  beta1,2  beta2,2 
#>     2.25     2.25
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 56 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1311.888 
#> 10  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.18813    0.01915   9.826  < 2e-16 ***
#> mu2       0.18929    0.01828  10.354  < 2e-16 ***
#> alpha1,1  0.33976    0.08497   3.999 6.37e-05 ***
#> alpha2,1  0.96210    0.14672   6.557 5.48e-11 ***
#> alpha1,2  1.04192    0.14743   7.067 1.58e-12 ***
#> alpha2,2  0.61487    0.12641   4.864 1.15e-06 ***
#> beta1,1   1.67851    0.40245   4.171 3.04e-05 ***
#> beta2,1   2.56714    0.42657   6.018 1.76e-09 ***
#> beta1,2   2.42639    0.37305   6.504 7.81e-11 ***
#> beta2,2   2.74344    0.60919   4.503 6.69e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

### More complicated model

#### Linear impact model

In the following, `impact` represents the impact of mark. It is a function, and the first argument is `param` represents the parameter of the model. `impact` function can have additional arguments such as `alpha`, `n`, `mark`, etc., which represents the dependence. Do not miss `...` as the ellipsis is omitted, an error occurs.

``` r
mu <- matrix(c(0.15, 0.15), nrow=2)
alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=T)
beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2)
rmark <- function(param = c(p=0.65), ...){
  rgeom(1, p=param[1]) + 1
}

impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
  ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
  alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
  #alpha * ma * matrix( c(rep(param["eta1"], 2), rep(param["eta2"], 2)), nrow=2)
}
h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
          rmark = rmark,
          impact=impact)
res <- hsim(h1, size=1000, lambda0 = matrix(rep(0.1,4), nrow=2))

fit <- hfit(h1, 
            inter_arrival = res$inter_arrival,
            type = res$type,
            mark = res$mark,
            lambda0 = matrix(rep(0.1,4), nrow=2))
#>      mu1 alpha1,1 alpha1,2  beta1,1     eta1 
#>     0.15     0.75     0.60     2.60     0.20

summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 46 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1526.811 
#> 5  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1       0.15046    0.00928  16.213  < 2e-16 ***
#> alpha1,1  0.80019    0.10186   7.856 3.96e-15 ***
#> alpha1,2  0.68370    0.08870   7.708 1.28e-14 ***
#> beta1,1   2.76112    0.25024  11.034  < 2e-16 ***
#> eta1      0.17281    0.10292   1.679   0.0931 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

#### Hawkes flocking model

In the basic model, `alpha` is a matrix, but it can be a function as in the following code. The function `alpha` simply return a 4 × 4 matrix but by doing so, we can fix some of elements as specific vales when estimating. When estimating, the optimization only works for the specified parameters in `param`. In the case of simulation, there is no difference whether the parameter set is represented by a matrix or a function.

``` r
mu <- matrix(c(0.02, 0.02, 0.04, 0.04), nrow = 4)


alpha <- function(param = c(alpha11 = 0, alpha12 = 0.3, alpha33 = 0.3, alpha34 = 0.4)){
  matrix(c(param["alpha11"], param["alpha12"], 0, 0,
           param["alpha12"], param["alpha11"], 0, 0,
           0, 0, param["alpha33"], param["alpha34"],
           0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = T)
}


beta <- matrix(c(rep(0.7, 8), rep(1.1, 8)), nrow = 4, byrow = T)

impact <- function(param = c(alpha1n=0, alpha1w=0.1, alpha2n=0.1, alpha2w=0.2),
                   n=n, N=N, ...){
  
  Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
                  0, 0, param['alpha1n'], param['alpha1w'],
                  param['alpha2w'], param['alpha2n'], 0, 0,
                  param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=T)
  
  ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5
  
  km <- matrix(c(!ind, !ind, !ind, !ind,
                 ind, ind, ind, ind,
                 ind, ind, ind, ind,
                 !ind, !ind, !ind, !ind), nrow = 4, byrow = T)
  
  km * Psi
}

h <- new("hspec",
         mu = mu, alpha = alpha, beta = beta, impact = impact)

hr <- hsim(h, size=1000)


fit <- hfit(h, hr$inter_arrival, hr$type)
#>     mu1     mu3 alpha11 alpha12 alpha33 alpha34 beta1,1 beta3,1 alpha1n 
#>    0.02    0.04    0.00    0.30    0.30    0.40    0.70    1.10    0.00 
#> alpha1w alpha2n alpha2w 
#>    0.10    0.10    0.20
summary(fit)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> BFGS maximization, 91 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -2627.193 
#> 12  free parameters
#> Estimates:
#>          Estimate Std. error t value  Pr(> t)    
#> mu1      0.022172   0.002502   8.863  < 2e-16 ***
#> mu3      0.040147   0.003421  11.736  < 2e-16 ***
#> alpha11 -0.045187   0.013858  -3.261  0.00111 ** 
#> alpha12  0.310024   0.044909   6.903 5.08e-12 ***
#> alpha33  0.231770   0.038615   6.002 1.95e-09 ***
#> alpha34  0.480105   0.052190   9.199  < 2e-16 ***
#> beta1,1  0.689674   0.081276   8.486  < 2e-16 ***
#> beta3,1  1.154105   0.103202  11.183  < 2e-16 ***
#> alpha1n -0.006219   0.017935  -0.347  0.72878    
#> alpha1w  0.110979   0.021764   5.099 3.41e-07 ***
#> alpha2n  0.230138   0.075533   3.047  0.00231 ** 
#> alpha2w  0.207419   0.060753   3.414  0.00064 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```
