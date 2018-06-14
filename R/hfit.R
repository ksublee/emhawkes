#' @include hspec.R hmoment.R hllf.R
NULL

#' Perform Maximum Likelihood Estimation
#'
#' Generic function hfit.
#'
#'
#' @param object hspec
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param lambda0 the inital values of lambda component. Must have the same dimensional matrix (n by n) with hspec.
#' @param N0 the initial values of N.
#' @param reduced When TRUE, reduced estimation performed.
#' @param constraint constraint matrix. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param method method for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param grad gradient matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param hess Hessian matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param verbose If true, mle progress is printed
#' @param ... other parameters for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#'
#'
#' @docType methods
#' @rdname hfit
#' @export
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hsim,hspec-method}}
#'
#' @examples
#' #example 1
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' summary(hfit(h, res$inter_arrival, res$type))
#'
#'
#' #example 2
#' mu <- matrix(c(0.08, 0.08, 0.05, 0.05), nrow = 4)
#' alpha <- function(param = c(alpha11 = 0, alpha12 = 0.4, alpha33 = 0.5, alpha34 = 0.3)){
#'   matrix(c(param["alpha11"], param["alpha12"], 0, 0,
#'            param["alpha12"], param["alpha11"], 0, 0,
#'            0, 0, param["alpha33"], param["alpha34"],
#'            0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = TRUE)
#' }
#' beta <- matrix(c(rep(0.6, 8), rep(1.2, 8)), nrow = 4, byrow = TRUE)
#'
#' impact <- function(param = c(alpha1n=0, alpha1w=0.2, alpha2n=0.001, alpha2w=0.1),
#'                    n=n, N=N, ...){
#'
#'   Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
#'                   0, 0, param['alpha1n'], param['alpha1w'],
#'                   param['alpha2w'], param['alpha2n'], 0, 0,
#'                   param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=TRUE)
#'
#'   ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5
#'
#'   km <- matrix(c(!ind, !ind, !ind, !ind,
#'                  ind, ind, ind, ind,
#'                  ind, ind, ind, ind,
#'                  !ind, !ind, !ind, !ind), nrow = 4, byrow = TRUE)
#'
#'   km * Psi
#' }
#' h <- new("hspec",
#'          mu = mu, alpha = alpha, beta = beta, impact = impact)
#' hr <- hsim(h, size=1000)
#' plot(hr$arrival, hr$N[,'N1'] - hr$N[,'N2'], type='s')
#' lines(hr$N[,'N3'] - hr$N[,'N4'], type='s', col='red')
#' fit <- hfit(h, hr$inter_arrival, hr$type)
#' summary(fit)
#'
#'
#' #example 3
#' mu <- c(0.15, 0.15)
#' alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=TRUE)
#' beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2, byrow=TRUE)
#' rmark <- function(param = c(p=0.65), ...){
#'   rgeom(1, p=param[1]) + 1
#' }
#' impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
#'   ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#'   alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#'   #alpha * ma * matrix( c(rep(param["eta1"], 2), rep(param["eta2"], 2)), nrow=2)
#' }
#' h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#'           rmark = rmark,
#'           impact=impact)
#' res <- hsim(h1, size=100, lambda0 = matrix(rep(0.1,4), nrow=2))
#'
#' fit <- hfit(h1,
#'             inter_arrival = res$inter_arrival,
#'             type = res$type,
#'             mark = res$mark,
#'             lambda0 = matrix(rep(0.1,4), nrow=2))
#' summary(fit)
#'
#'
#' #example 4
#' mu <- function(param = c(mu1 = 0.08, eta1 = 0.7), n=n, N=N, ...){
#'   if(n == 1){
#'     level <- N[,"N1"][1] - N[,"N2"][1] - (N[,"N3"][1] - N[,"N4"][1])
#'     matrix(c(param["mu1"], param["eta1"]*level, param["eta1"]*level, param["mu1"]), nrow = 4)
#'   } else {
#'     level <- N[,"N1"][n-1] - N[,"N2"][n-1] - (N[,"N3"][n-1] - N[,"N4"][n-1])
#'     matrix( c(param["mu1"], param["eta1"]*level, param["eta1"]*level, param["mu1"]), nrow = 4)
#'   }
#' }
#' alpha <- function(param = c(alpha11 = 0.6, alpha14=0.7)){
#'   matrix(c(param["alpha11"], 0, 0, param["alpha14"],
#'            0, 0, 0, 0,
#'            0, 0, 0, 0,
#'            param["alpha14"], 0, 0, param["alpha11"]), nrow = 4, byrow = TRUE)
#' }
#
#' beta <- matrix(rep(2.6, 16), nrow=4, byrow=TRUE)
#' h <- new("hspec", mu, alpha, beta)
#' hr <- hsim(h, size=100)
#'
#' fit <- hfit(h, hr$inter_arrival, hr$type)
#' summary(fit)
#'
setGeneric("hfit", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL, lambda0 = NULL, N0 = NULL,
                            reduced = TRUE,
                            grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
                            verbose = FALSE, ...) standardGeneric("hfit"))


#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#'
#' @rdname hfit
setMethod(
  f="hfit",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL, lambda0 = NULL, N0 = NULL,
           reduced = TRUE,
           grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
           verbose = FALSE, ...){


    if(is.null(lambda0)){
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
    }

    # When the mark sizes are not provided, the jumps are all unit jumps.
    # unit <- FALSE
    # if(is.null(mark)) {
    #   mark <- c(0, rep(1, length(inter_arrival)-1))
    #   unit <- TRUE
    # }

    # parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens


    # parameter setting
    if(is.function(mu)){
      pr_mus <- eval(formals(object@mu)[[1]])
    }else{
      pr_mus <- as.param(object@mu, "mu", reduced)
    }

    pr_alphas <- as.param(object@alpha, "alpha", reduced)
    pr_betas <- as.param(object@beta, "beta", reduced)

    len_mu <- length(pr_mus)
    len_alpha <- length(pr_alphas)
    len_beta <- length(pr_betas)

    if(!is.null(impact)){
      pr_impact <- eval(formals(object@impact)[[1]])
    } else {
      pr_impact <- NULL
    }
    starting_point <- c(pr_mus, pr_alphas, pr_betas, pr_impact)
    print(starting_point)

    # loglikelihood function for maxLik

    llh_function <- function(param){


      pr_mus <- param[1:len_mu]
      pr_alphas <- param[(len_mu + 1):(len_mu + len_alpha)]
      pr_betas <- param[(len_mu + len_alpha + 1):(len_mu + len_alpha + len_beta)]
      pr_impact <- param[(len_mu + len_alpha + len_beta + 1):length(param)]


      if (is.function(object@mu)){
        mu0 <- hijack(object@mu, param = pr_mus)
      } else{
        if(reduced){
          mu0 <- matrix(pr_mus[look_up_mtrx(mu, "mu")], nrow=dimens)
        } else {
          mu0 <- matrix(pr_mus, nrow=dimens)
        }
      }
      if (is.function(object@alpha)){
        alpha0 <- hijack(object@alpha, param = pr_alphas)
      } else{
        if(reduced){
          alpha0 <- matrix(pr_alphas[look_up_mtrx(alpha, "alpha")], nrow=dimens)
        } else {
          alpha0 <- matrix(pr_alphas, nrow=dimens)
        }
      }
      if (is.function(object@beta)){
        beta0 <- hijack(object@beta, param = pr_betas)
      } else{
        if(reduced){
          beta0 <- matrix(pr_betas[look_up_mtrx(beta, "beta")], nrow=dimens)
        } else {
          beta0 <- matrix(pr_betas, nrow=dimens)
        }
      }

      #object@impact is user defined impact function
      if (!is.null(object@impact)){
        impact0 <- hijack(object@impact, param = pr_impact)

        hspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0,
                                impact = impact0,
                                rmark = object@rmark)

      } else {
        hspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0,
                                rmark = object@rmark)
      }

      logl <- logLik(hspec0, inter_arrival = inter_arrival, type = type,
                     mark = mark, N0 = N0, lambda0 = lambda0)
      if(verbose){
        cat("Parameters : ", param, "\n")
        cat("Log likelihood : ", logl, "\n")
      }

      logl

    }

    #llh_function(starting_point)
    maxLik::maxLik(logLik=llh_function,
                   start=starting_point, grad, hess, method = method)

  }
)

