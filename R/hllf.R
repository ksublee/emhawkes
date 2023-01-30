#' @include hspec.R hmoment.R
NULL

#' Compute the log-likelihood function
#'
#' The log-likelihood of the ground process of the Hawkes model.
#' (The log-likelihood for mark (jump) distribution is not provided.)
#'
#' @param object \code{\link{hspec-class}}. The parameter values in the object are used to compute the log-likelihood.
#' @param inter_arrival A vector of realized inter-arrival times of events whici includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type A vector of realized dimensions distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark A vector of realized mark (jump) sizes. Start with zero.
#' @param N A matrix of counting processes.
#' @param Nc A matrix of counting processes weighted by mark.
#' @param lambda_component0 The initial values of lambda component. Must have the same dimensional matrix with \code{object}.
#' @param N0 A matrix of initial values of N.
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hfit,hspec-method}}
#'
#' @docType methods
#' @rdname logLik
#'
#' @export
setMethod(
  f="logLik",
  signature(object="hspec"),
  function(object, inter_arrival, type = NULL, mark = NULL, N = NULL, Nc = NULL,
           N0 = NULL, lambda_component0 = NULL, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0.")

      lambda_component0 <- additional_argument[["lambda0"]]

    }



    # parameter setting
    # after parameter setting, functions mu, alpha, beta become matrices


    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta
    impact <- plist$impact
    rmark <- plist$rmark
    dmark <- plist$dmark
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

    # check the argument lists in mu and impact
    mu_args <- c()
    impct_args <- c()

    if (is.function(mu)){
      if(length(formals(mu)) > 1){
        mu_args <- methods::formalArgs(mu)
      }
    }
    if(!is.null(impact)){
      impct_args <- methods::formalArgs(impact)
    }


    # size is length(inter_arrival) - 1
    size <- length(inter_arrival)

    if("N" %in% c(impct_args, mu_args) & is.null(N)){
      N <- matrix(numeric(length = dimens * size), ncol = dimens)
      colnames(N) <- paste0("N", 1:dimens)
      if (!is.null(N0)){
        N[1,] <- N0
      }
      N[cbind(2:size, type[-1])] <- 1
      N <- apply(N, 2, cumsum)
    }

    if("Nc" %in% c(impct_args, mu_args) & is.null(Nc)){
      Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)
      colnames(Nc)  <- paste0("Nc", 1:dimens)
      if (!is.null(N0)){
        Nc[1,] <- N0
      }
      Nc[cbind(2:size, type[-1])] <- mark[-1]
      Nc <- apply(Nc, 2, cumsum)
    }



    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                lambda_component_n = lambda_component_n,
                alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      mu0 <- mu
    }

    # default lambda_component0
    if(length(object@type_col_map) > 0 & is.null(lambda_component0)){
      stop("In this model, please provide lambda_component0.")
    }
    if(is.null(lambda_component0)) {
      if(!("showWarning" %in% names(additional_argument))){

        message("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")

      }

      lambda_component0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu = mu, alpha = alpha, beta = beta)
    }
    rowSums_lambda_component0 <- rowSums(matrix(lambda_component0, nrow=dimens))


    # little bit complicated code
    if("lambda_component" %in% c(impct_args, mu_args)){
      lambda_component <- matrix(sapply(lambda_component0, c, numeric(length = size - 1)), ncol = ncol(beta) * dimens)
      indxM <- matrix(rep(1:ncol(beta), dimens), byrow = TRUE, nrow = dimens)
      #colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
      colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:dimens, FUN = paste0)),
                                                      1:ncol(beta), FUN=paste0)))
    }

    if("lambda" %in% c(impct_args, mu_args)){
      lambda   <- matrix(sapply(as.vector(mu0) + rowSums_lambda_component0, c, numeric(length = size - 1)), ncol = dimens)
      colnames(lambda) <- paste0("lambda", 1:dimens)
    }


    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0
    sum_mu_inter_arrival <- 0
    sum_log_dmark <- 0

    current_lambda_component_without_mu <- lambda_component0

    # currently only piecewise constant mu is available, hope to be updated
    if (is.function(mu)){
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                 lambda_component_n = lambda_component_n,
                 alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }


    is_lambda_necessary <- "lambda" %in% c(impct_args, mu_args)
    is_lambda_component_necessary <- "lambda_component" %in% c(impct_args, mu_args)

    zero_mat_for_alpha <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize
    zero_mat_for_eta <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize
    zero_mat_for_impact <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize


    for (n in 2:size) {

      mu_n <- rmu_n
      type_n <- type[n]
      inter_arrival_n <- inter_arrival[n]

      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival_n)
      decayed_lambda <- lambda_component_n <- current_lambda_component_without_mu * decayed

      ## 1. sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component +
        sum(current_lambda_component_without_mu / beta * ( 1 - decayed ))


      ## 2. sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- mu_n + sum(decayed_lambda)
      else lambda_lc_type_n <- mu_n[type_n] + sum(decayed_lambda[type_n, ])


      #log(lambda_lc[type[n]]) can be NaN, so warning is turned off for a moment
      #oldw <- getOption("warn")

      #options(warn = -1)
      if (dimens == 1) sum_log_lambda <- sum_log_lambda + suppressWarnings(log(lambda_lc))
      else sum_log_lambda <- sum_log_lambda + suppressWarnings(log(lambda_lc_type_n))


      #options(warn = oldw)

      ## 2.1. sum of log mark p.d.f. when jump occurs
      if(!is.null(dmark)){
        sum_log_dmark <- sum_log_dmark +
          log(dmark(mark = mark, n = n, type = type, inter_arrival = inter_arrival,
                    N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                    lambda_component_n = lambda_component_n,
                    mu = mu, alpha = alpha, beta = beta))

      }



      ## 3. sum of mu * inter_arrival
      sum_mu_inter_arrival <- sum_mu_inter_arrival + sum(mu_n) * inter_arrival_n

      if(is_lambda_component_necessary)
        lambda_component[n, ] <- t(decayed_lambda)
      if(is_lambda_necessary)
        lambda[n, ] <- as.vector(mu_n) + rowSums(decayed_lambda)

      # impact

      if( is.null(object@type_col_map) ){
        types <- type_n
      } else if ( length(object@type_col_map) > 0 ) {
        types <- object@type_col_map[[type_n]]
      } else {
        stop("Check the type_col_map argument.")
      }

      # 1. impact by alpha

      impact_alpha <- zero_mat_for_alpha   # matrix
      impact_alpha[ , types] <- alpha[ , types]

      current_lambda_component_without_mu <- decayed_lambda + impact_alpha

      # 2. additional impact by eta

      if(!is.null(eta)) {

        impact_eta <- zero_mat_for_eta   # matrix
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        current_lambda_component_without_mu <- current_lambda_component_without_mu + impact_eta
      }


      # 3. impact by impact function
      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){

        impact_mark <- zero_mat_for_impact
        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                             lambda_component_n = lambda_component_n,
                             mu = mu, alpha = alpha, beta = beta)

        impact_mark[ , types] <- impact_res[ , types]

        # for next step
        current_lambda_component_without_mu <- current_lambda_component_without_mu + impact_mark

      }


      # new mu, i.e., right continuous version of mu
      if (is.function(mu)){
        # mu is represented by a function
        rmu_n <- mu(n = n + 1, mark = mark, type = type, inter_arrival = inter_arrival,
                   N = N, Nc = Nc,
                   alpha = alpha, beta = beta)
      } else{
        # mu is a matrix
        rmu_n <- mu
      }


    }


    # log likelihood for ground process
    # sum_log_lambda - sum(mu_n*sum(inter_arrival)) - sum_integrated_lambda_component
    sum_log_lambda - sum_mu_inter_arrival - sum_integrated_lambda_component + sum_log_dmark

  }
)



