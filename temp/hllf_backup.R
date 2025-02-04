#' @include hspec.R hmoment.R
NULL

#' Compute the Log-Likelihood Function
#'
#' Calculates the log-likelihood for the Hawkes model.
#'
#' @param object An \code{\link{hspec-class}} object containing parameter values for computing the log-likelihood.
#' @param inter_arrival A vector of inter-arrival times for events across all dimensions, starting with zero.
#' @param type A vector indicating the dimensions, represented by numbers (1, 2, 3, etc.), starting with zero.
#' @param mark A vector of mark (jump) sizes, starting with zero.
#' @param N A matrix representing counting processes.
#' @param Nc A matrix of counting processes weighted by mark sizes.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' Can be a numeric value or a matrix.
#' Must have the same number of rows and columns as \code{alpha} or \code{beta} in \code{object}.
#' @param N0 A matrix of initial values for \code{N}.
#' @param Nc0 A matrix of initial values for \code{Nc}.
#' @param infer
#' @param ... Additional arguments passed to or from other methods.
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
           N0 = NULL, Nc0 = NULL, lambda_component0 = NULL,
           infer = FALSE, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0.")

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    ## Parameter setting and initialization
    # after parameter setting, functions mu, alpha, beta become matrices, if possible
    # Note that mu, alpha, beta, eta can be given functions even they are essentially matrices.
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta

    size <- length(inter_arrival)

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) mark <- rep(1, size)

    if(object@dimens != 1 && is.null(type)) stop("The argument type should be provided.")
    if(object@dimens == 1 && is.null(type)) type <- rep(1, size)

    # check the argument lists in mu and impact
    mu_args <- impct_args <- c()
    if(is.function(mu) && length(formals(mu)) > 1) mu_args <- methods::formalArgs(mu)
    if(!is.null(object@impact)) impct_args <- methods::formalArgs(object@impact)

    # N and Nc are constructed only if they are needed.
    if("N" %in% c(impct_args, mu_args) && is.null(N)) N <- type_to_N(type, object@dimens, N0 = N0)
    if("Nc" %in% c(impct_args, mu_args) && is.null(Nc)) Nc <- type_mark_to_Nc(type, mark, object@dimens, Nc0 = Nc0)


    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                N = N, Nc = Nc,
                lambda = lambda, lambda_component = lambda_component,
                lambda_component_n = lambda_component_n,
                alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      mu0 <- mu
    }

    # default lambda_component0
    if(length(object@type_col_map) > 0 && is.null(lambda_component0)){
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


    is_lambda_necessary <- infer || "lambda" %in% c(impct_args, mu_args)
    is_lambda_component_necessary <- infer || "lambda_component" %in% c(impct_args, mu_args)


    if (is_lambda_component_necessary) {
      # Initialize the lambda_component matrix with zeros
      lambda_component <- matrix(0, nrow = size, ncol = ncol(beta) * object@dimens)

      # Set the initial values for the first row based on lambda_component0
      lambda_component[1, ] <- t(lambda_component0)

      # Generate column names using outer and paste
      colnames(lambda_component) <- paste0("lambda", rep(1:object@dimens, each = ncol(beta)), 1:ncol(beta))
    }

    if(is_lambda_necessary){
      # Initialize the lambda matrix with zeros
      lambda <- matrix(0, nrow = size, ncol = object@dimens)

      # Set the initial values for the first row
      lambda[1, ] <- as.vector(mu0) + rowSums(matrix(lambda_component0, nrow=object@dimens))

      # Set column names
      colnames(lambda) <- paste0("lambda", 1:object@dimens)
    }

    if (infer) {
      # Initialize the rambda_component matrix with zeros
      rambda_component <- matrix(0, nrow = size, ncol = ncol(beta) * object@dimens)

      # Set the initial values for the first row based on lambda_component0
      rambda_component[1, ] <- t(lambda_component0)

      # Generate column names using outer and paste
      colnames(rambda_component) <- paste0("rambda", rep(1:object@dimens, each = ncol(beta)), 1:ncol(beta))

      # Initialize the rambda matrix with zeros
      rambda <- matrix(0, nrow = size, ncol = object@dimens)

      # Set the initial values for the first row
      rambda[1, ] <- as.vector(mu0) + rowSums(matrix(lambda_component0, nrow=object@dimens))

      # Set column names
      colnames(rambda) <- paste0("rambda", 1:object@dimens)

      # for residual
      phis <- matrix(nrow = length(inter_arrival), ncol = object@dimens)
      colnames(phis) <- paste0("phi", 1:object@dimens)


    }


    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0
    sum_mu_inter_arrival <- 0
    sum_log_f_epsilon_phi <- 0
    sum_log_dmark <- 0


    # Only piece-wise constant mu is available
    if (is.function(mu)){
      # for inter-arrival between n=1 and n=2
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                 lambda_component_n = lambda_component_n,
                 alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }


    zero_mat_for_alpha <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize
    zero_mat_for_eta <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize
    zero_mat_for_impact <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize


    rambda_component_n <- lambda_component0


    # handling both basic and discrete models
    dresidual <- if(is.null(object@dresidual)) stats::dexp else object@dresidual
    presidual <- if(is.null(object@presidual)) stats::pexp else object@presidual

    ## main loop
    for (n in 2:size) {

      mu_n <- rmu_n; rambda_component_n_minus_1 <- rambda_component_n

      type_n <- type[n]

      decayed <- exp(-beta * inter_arrival[n])

      # lambda decayed due to time, impact due to mark is not added yet
      lambda_component_n <- rambda_component_n_minus_1 * decayed

      # integrate mu over inter_arrival w.r.t. time
      # mu is matrix or vector
      if ((is.matrix(mu_n) && nrow(mu_n) == object@dimens && ncol(mu_n) == 1) ||
          (is.vector(mu_n) && length(mu_n) == 1)) {

        integrated_mu <- sum(mu_n) * inter_arrival[n]

      } else {

        # When mu is observed as piecewise constant between n-1 and n
        # Rows from 1 to object@dimens are mu values for each dimension.
        # The last row represents the inter-arrival where mu changes between n-1 and n.

        sum_mu_n <- colSums(mu_n[1:object@dimens, , drop = FALSE])
        integrated_mu <- sum(sum_mu_n * mu_n[object@dimens + 1, ])

      }


      if (is.null(object@dresidual) && !infer){

        # basic model
        ## 1. sum of integrated_lambda_component and sum of mu * inter_arrival
        sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi - sum(rambda_component_n / beta * (1 - decayed)) -
          integrated_mu


      } else {


        # discrete model or when needed to infer
        for (k in 1:object@dimens){

          phi_k_wo_mu <- phi_wo_mu(tau = inter_arrival[n],
                                   rowSum_rambda_component_n_miuns_1 = sum(rambda_component_n_minus_1[k, ]),
                                   beta[k, 1])


          #phi_value_k <- phi(tau = inter_arrival[n],
          #                   rowSum_rambda_component_n_miuns_1 = sum(rambda_component_n_minus_1[k, ]),
          #                   mu = mu_n[k], beta[k, 1])

          phi_value_k <- phi_k_wo_mu + integrated_mu

          log_value <- suppressWarnings(if (k == type_n) log(dresidual(x=phi_value_k)) else log(1-presidual(q=phi_value_k)))

          sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi + log_value

          if(infer) phis[n, k] <- phi_value_k

        }

      }



      ## 2. sum of log lambda when jump occurs
      # Calculate sum of log lambda when a jump occurs
      lambda_at_event <- if (object@dimens == 1) {
        mu_n + lambda_component_n
      } else {
        mu_n[type_n] + sum(lambda_component_n[type_n, ])
      }

      # Update the sum of log lambda with a warning suppression for log calculation
      sum_log_lambda <- sum_log_lambda + suppressWarnings(log(lambda_at_event))


      if(is_lambda_component_necessary) lambda_component[n, ] <- t(lambda_component_n)
      if(is_lambda_necessary){
        lambda[n, ] <- as.vector(mu_n[1:object@dimens]) + rowSums(lambda_component_n)
      }

      ## 3. sum of log mark p.d.f. when jump occurs
      if(!is.null(object@dmark)){
        sum_log_dmark <- sum_log_dmark +
          log(object@dmark(mark = mark, n = n, type = type, inter_arrival = inter_arrival,
                           N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                           lambda_component_n = lambda_component_n,
                           mu = mu, alpha = alpha, beta = beta))

      }



      ## impact

      # Determine types based on the type_col_map
      types <- if (is.null(object@type_col_map)) {
        type_n
      } else if (length(object@type_col_map) > 0) {
        object@type_col_map[[type_n]]
      } else {
        stop("Check the type_col_map argument.")
      }

      # 1. impact by alpha

      impact_alpha <- zero_mat_for_alpha   # matrix
      impact_alpha[ , types] <- alpha[ , types]

      rambda_component_n <- lambda_component_n + impact_alpha

      # 2. additional impact by eta

      if(!is.null(eta)) {

        impact_eta <- zero_mat_for_eta   # matrix
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        rambda_component_n <- rambda_component_n + impact_eta
      }


      # 3. impact by impact function
      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(object@impact)){

        impact_mark <- zero_mat_for_impact
        impact_res <- object@impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                                    N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                                    lambda_component_n = lambda_component_n,
                                    mu = mu, alpha = alpha, beta = beta)

        impact_mark[ , types] <- impact_res[ , types]

        # for next step
        rambda_component_n <- rambda_component_n + impact_mark

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

      if(infer){

        rambda_component[n, ] <- t(rambda_component_n)
        rambda[n, ] <- as.vector(rmu_n[1:object@dimens]) + rowSums(rambda_component_n)

      }

    }


    if(!infer) {
      # only log-likelihood
      c(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark)

    } else {

      return_object <- list(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark,
                            lambda = lambda, lambda_component = lambda_component,
                            rambda = rambda, rambda_component = rambda_component)

      if (is.null(object@dresidual) || object@dimens == 1 ){


        for (p in 1:object@dimens){

          np <- which(type == p)
          res_process <- rep(0, length(np) - 1)

          for(k in seq_along(res_process)){
            res_process[k] <- sum(phis[(np[k]+1):np[k+1], p])
          }

        }


        return_object[[paste0("res_process", p)]] <- res_process



      } else if (!is.null(object@dresidual)){

        for (p in 1:object@dimens){

          np <- which(type == p)
          res_process <- rep(0, length(np) - 1)

          for(k in seq_along(res_process)){
            res_process[k] <- phis[(np[k]+1), p]
          }

          return_object[[paste0("min_res_process", p)]] <- res_process

        }

      }

      return_object

    }

  }
)



