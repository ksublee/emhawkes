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
  function(object, inter_arrival, type = NULL, mark = NULL, N = NULL, Nc = NULL,
           N0 = NULL, lambda0 = NULL){
    # parameter setting
    # after parameter setting, functions mu, alpha, beta become matrices
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

    # default lambda0
    if(length(type_col_map) > 0 & is.null(lambda0)){
      stop("In this model, please provide lambda0.")
    }
    if(is.null(lambda0)) {
      if(!exists("this_flag_represents_binding_env_is_hfit")){
        warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      }
      lambda0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu = mu, alpha = alpha, beta = beta)
    }
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))


    # too complicated code
    if("lambda_component" %in% c(impct_args, mu_args)){
      lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = ncol(beta) * dimens)
      indxM <- matrix(rep(1:ncol(beta), dimens), byrow = TRUE, nrow = dimens)
      #colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
      colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:dimens, FUN = paste0)),
                                                      1:ncol(beta), FUN=paste0)))
    }

    if("lambda" %in% c(impct_args, mu_args)){
      lambda   <- matrix(sapply(as.vector(mu0) + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)
      colnames(lambda) <- paste0("lambda", 1:dimens)
    }


    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0
    sum_mu_inter_arrival <- 0

    current_lambda <- lambda0

    # currently only piecewise constant mu is available, hope to be updated
    if (is.function(mu)){
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                 lambda_component_n = lambda_component_n,
                 alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }

    for (n in 2:size) {

      mu_n <- rmu_n

      decayed <- exp(-beta * inter_arrival[n])
      decayed_lambda <- lambda_component_n <- current_lambda * decayed

      ## 1. sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component +
        sum(current_lambda / beta * ( 1 - decayed ))

      ## 2. sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- mu_n + sum(decayed_lambda)
      else lambda_lc <- mu_n + rowSums(decayed_lambda)

      #log(lambda_lc[type[n]]) can be NaN, so warning is turned off for a moment
      oldw <- getOption("warn")
      options(warn = -1)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[type[n]])
      options(warn = oldw)

      ## 3. sum of mu * inter_arrival
      sum_mu_inter_arrival <- sum_mu_inter_arrival + sum(mu_n) * inter_arrival[n]

      if("lambda_component" %in% c(impct_args, mu_args))
        lambda_component[n, ] <- t(decayed_lambda)
      if("lambda" %in% c(impct_args, mu_args))
        lambda[n, ] <- as.vector(mu_n) + rowSums(decayed_lambda)

      # impact
      # 1. impact by alpha
      impact_alpha <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)

      if( length(object@type_col_map) == 0){
        types <- type[n]
      } else{
        types <- object@type_col_map[[type[n]]]
      }

      impact_alpha[ , types] <- alpha[ , types]

      # 2. impact by mark
      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){

        impact_mark <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)
        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                             lambda_component_n = lambda_component_n,
                             mu = mu, alpha = alpha, beta = beta)

        #impact_res <- impact()

        impact_mark[ , types] <- impact_res[ , types]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }


      # for next step
      # current_lambda <- matrix(lambda_component[n, ], nrow = dimens, byrow = TRUE)
      current_lambda <- new_lambda

      # new mu, i.e., right continuous version of mu
      if (is.function(mu)){
        # mu is represeted by function
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
    sum_log_lambda - sum_mu_inter_arrival - sum_integrated_lambda_component


  }
)



