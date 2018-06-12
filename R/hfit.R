#' @include hspec.R hmoment.R
NULL

#' Compute the loglikelihood function
#'
#' The loglikelihood of the ground process of the Hawkes model.
#' (The estimation for jump distribution is not provided.)
#'
#' @param object \code{\link{hspec-class}}. The parameter values in the object are used to compute the log-likelihood.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param lambda0 the initial values of lambda component. Must have the same dimensional matrix with \code{hspec}.
#' @param N0 the initial value of N
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hfit,hspec-method}}
#'
#' @export
setMethod(
  f="logLik",
  signature(object="hspec"),
  function(object, inter_arrival, type = NULL, mark = NULL, N0 = NULL, lambda0 = NULL){
    # parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) {
      mark <- rep(1, length(inter_arrival))
    }
    # if dimens == 1 and type is not provided, then all type is 1.
    if(dimens==1 & is.null(type)) {
      type <- rep(1, length(inter_arrival))
    } else if (dimens != 1 & is.null(type)) {
      stop("The argument type should be provided.")
    }

    mu_args <- c()
    impct_args <- c()

    if (is.function(mu) & length(formals(mu)) > 1){
      mu_args <- methods::formalArgs(mu)
    }
    if(!is.null(impact)){
      impct_args <- methods::formalArgs(impact)
    }


    # size is length(inter_arrival) - 1
    size <- length(inter_arrival)

    if("N" %in% impct_args | "N" %in% mu_args){
      N <- matrix(numeric(length = dimens * size), ncol = dimens)
      colnames(N) <- paste0("N", 1:dimens)
      if (!is.null(N0)){
        N[1,] <- N0
      }
    }

    if("Nc" %in% impct_args | "Nc" %in% mu_args){
      Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)
      colnames(Nc)  <- paste0("Nc", 1:dimens)
      if (!is.null(N0)){
        Nc[1,] <- N0
      }
    }

    if (is.function(mu)){
      # mu is represeted by function
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                N = N, Nc = Nc,
                alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      mu0 <- mu
    }

    # default lambda0
    if(is.null(lambda0)) {
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             alpha = alpha, beta = beta)
    }
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))


    if("lambda_component" %in% impct_args | "lambda_component" %in% mu_args){
      lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = dimens^2)
      indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
      colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
    }

    if("lambda" %in% impct_args | "lambda" %in% mu_args){
      lambda   <- matrix(sapply(mu0 + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)
      colnames(lambda) <- paste0("lambda", 1:dimens)
    }

    #if (dimens==1) rowSums_lambda0 <- lambda0
    #else rowSums_lambda0 <- rowSums(lambda0)

    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0
    sum_mu_inter_arrival <- sum(mu0) * inter_arrival[1]

    current_lambda <- lambda0

    for (n in 2:size) {

      if (is.function(mu)){
        # mu is represeted by function
        mu_n <- mu(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                  N = N, Nc = Nc,
                  alpha = alpha, beta = beta)
      } else{
        # mu is a matrix
        mu_n <- mu
      }

      decayed <- exp(-beta * inter_arrival[n])
      decayed_lambda <- current_lambda * decayed

      if("lambda_component" %in% impct_args | "lambda_component" %in% mu_args)
        lambda_component[n, ] <- t(decayed_lambda)
      if("lambda" %in% impct_args | "lambda" %in% mu_args)
        lambda[n, ] <- mu + rowSums(decayed_lambda)

      if("N" %in% impct_args | "N" %in% mu_args){
        N[n, ] <- N[n-1, ]
        N[n, type[n]] <- N[n-1, type[n]] + 1
      }
      if("Nc" %in% impct_args | "Nc" %in% mu_args){
        Nc[n, ] <- Nc[n-1, ]
        Nc[n, type[n]] <- Nc[n-1, type[n]] + mark[n]
      }

      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens^2), nrow = dimens)
      impact_alpha[ , type[n]] <- alpha[ , type[n]]

      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens^2), nrow = dimens)
        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                             mu = mu, alpha = alpha, beta = beta)

        #impact_res <- impact()

        impact_mark[ , type[n]] <- impact_res[ , type[n]]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }

      # sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component +
        sum(current_lambda / beta * ( 1 - decayed ))

      # sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- mu_n + decayed_lambda
      else lambda_lc <- mu_n + rowSums(decayed_lambda)

      #log(lambda_lc[type[n]]) can be NaN, so warning is turned off for a moment
      oldw <- getOption("warn")
      options(warn = -1)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[type[n]])
      options(warn = oldw)


      # sum of mu * inter_arrival
      sum_mu_inter_arrival <- sum_mu_inter_arrival + sum(mu_n) * inter_arrival[n]

      # current_lambda <- matrix(lambda_component[n, ], nrow = dimens, byrow = TRUE)
      current_lambda <- new_lambda  # lambda determined in the previous loop

    }

    # log likelihood for ground process
    # sum_log_lambda - sum(mu_n*sum(inter_arrival)) - sum_integrated_lambda_component
    sum_log_lambda - sum_mu_inter_arrival - sum_integrated_lambda_component

  }
)

#' Perform a maximum likelihood estimation
setGeneric("hfit", function(object, ...) standardGeneric("hfit"))
#'
#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#'
#' @param object hspec, or can be omitted.
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
#' @seealso \code{\link{hspec-class}}, \code{\link{hsim,hspec-method}}
#'
#' @rdname hfit
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
#' @export
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

    pr_impact <- eval(formals(object@impact)[[1]])
    starting_point <- c(pr_mus, pr_alphas, pr_betas, pr_impact)


    # constraint matrix
    # mu, alpha, beta should be larger than zero
    #if (unit) A <- diag(1, nrow = length(starting_point) - pr_impact)
    #else A <- cbind(diag(1, nrow = length(starting_point) - length(unique_etas)), rep(0, length(starting_point) - length(unique_etas)))


    # constraint : sum of alpha < beta
    #A <- rbind(A, c(0, rep(-1, len_alpha), 1, rep(0, len_eta)))
    #B <- rep(0, nrow(A))


    # loglikelihood function for maxLik

    llh_function <- function(param){


      pr_mus <- param[1:len_mu]
      pr_alphas <- param[(len_mu + 1):(len_mu + len_alpha)]
      pr_betas <- param[(len_mu + len_alpha + 1):(len_mu + len_alpha + len_beta)]
      pr_impact <- param[(len_mu + len_alpha + len_beta + 1):length(param)]


      if (is.function(object@mu)){
        mu0 <- hijack(object@mu, param = pr_mus)
      } else{
        mu0 <- matrix(pr_mus[look_up_mtrx(mu, "mu")], nrow=dimens)
      }
      if (is.function(object@alpha)){
        alpha0 <- hijack(object@alpha, param = pr_alphas)
      } else{
        alpha0 <- matrix(pr_alphas[look_up_mtrx(alpha, "alpha")], nrow=dimens)
      }
      if (is.function(object@beta)){
        beta0 <- hijack(object@beta, param = pr_betas)
      } else{
        beta0 <- matrix(pr_betas[look_up_mtrx(beta, "beta")], nrow=dimens)
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

