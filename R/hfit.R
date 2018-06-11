#' @include hspec.R hmoment.R
#'
#' Compute the loglikelihood function
#'
#' This is a generic function.
#' The loglikelihood of the ground process of the Hawkes model.
#' (The estimation for jump distribution is not provided.)
#'
#' @param object \code{\link{hspec-class}}. The parameter values in the object are used to compute the log-likelihood.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param lambda0 The starting values of lambda. Must have the same dimensional matrix with \code{hspec}.
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hfit,hspec-method}}
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


    # default lambda0
    if(is.null(lambda0)) {
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object)
    }

    # size is length(inter_arrival) - 1
    size <- length(inter_arrival)

    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))


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
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[type[n]])

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

setGeneric("hfit", function(object, ...) standardGeneric("hfit"))

#' Perform a maximum likelihood estimation
#'
#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#'
#' @param object mhspec, or can be omitted.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param lambda0 the starting values of lambda. Must have the same dimensional matrix (n by n) with mhspec.
#' @param reduced When TRUE, reduced estimation performed.
#' @param constraint constraint matrix. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param method method for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param grad gradient matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param hess Hessian matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param ... other parameters for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#'
#' @examples
#' @seealso \code{\link{mhspec-class}}, \code{\link{mhsim,mhspec-method}}
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

        mhspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0,
                                impact = impact0,
                                rmark = object@rmark)

      } else {
        mhspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0,
                                rmark = object@rmark)
      }

      logl <- logLik(mhspec0, inter_arrival = inter_arrival, type = type,
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

