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
#' @param infer Logical
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
    # mu, alpha, beta, eta
    list2env(plist, envir = environment())

    size <- length(inter_arrival)

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) mark <- rep(1, size)

    if (is.null(type)) {
      if (object@dimens != 1) stop("The argument type should be provided.")
      type <- rep(1, size)
    }

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

    # ---- Lambda component0 check & setup ----
    if (is.null(lambda_component0)) {
      if (length(object@type_col_map) > 0) {
        stop("In this model, please provide lambda_component0.")
      }

      if (!("showWarning" %in% names(additional_argument))) {
        message("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      }

      lambda_component0 <- get_lambda0(
        object,
        mark = mark,
        type = type,
        inter_arrival = inter_arrival,
        N = N, Nc = Nc,
        mu = mu,
        alpha = alpha,
        beta = beta
      )
    }

    # ---- Determine necessity for lambda/lambda_component ----
    is_lambda_necessary <- infer || "lambda" %in% c(impct_args, mu_args)
    is_lambda_component_necessary <- infer || "lambda_component" %in% c(impct_args, mu_args)

    # ---- lambda_component ----
    if (is_lambda_component_necessary) {
      lambda_component <- initialize_lambda_matrix(
        prefix = "lambda",
        size = size,
        dimens = object@dimens,
        beta_cols = ncol(beta),
        lambda0 = lambda_component0
      )
    }

    # ---- lambda ----
    if (is_lambda_necessary) {
      lambda <- initialize_lambda_total(
        prefix = "lambda",
        size = size,
        dimens = object@dimens,
        mu0 = mu0,
        lambda0 = lambda_component0
      )
    }

    # ---- rambda (inference mode) ----
    if (infer) {
      rambda_component <- initialize_lambda_matrix(
        prefix = "rambda",
        size = size,
        dimens = object@dimens,
        beta_cols = ncol(beta),
        lambda0 = lambda_component0
      )

      rambda <- initialize_lambda_total(
        prefix = "rambda",
        size = size,
        dimens = object@dimens,
        mu0 = mu0,
        lambda0 = lambda_component0
      )

      # phi matrix (phi(tau) = epsilon), so this is related to residual process
      phis <- matrix(nrow = length(inter_arrival), ncol = object@dimens)
      colnames(phis) <- paste0("phi", 1:object@dimens)

      # historical right-continuous version of mus
      rmus <- matrix(nrow = length(inter_arrival), ncol = object@dimens)
      rmus[1, ] <-  if (is.function(mu)) {
        mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
           N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
           lambda_component_n = lambda_component_n,
           alpha = alpha, beta = beta)
      } else {
        mu
      }
    }


    # ---- Initialize accumulators ----
    accumulators <- c("sum_log_lambda", "sum_integrated_lambda_component",
                      "sum_mu_inter_arrival", "sum_log_f_epsilon_phi",
                      "sum_log_dmark")

    for (name in accumulators) assign(name, 0)


    # ---- Evaluate piecewise constant mu (only supported form) ----
    rmu_n <- if (is.function(mu)) {
      mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
         N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
         lambda_component_n = lambda_component_n,
         alpha = alpha, beta = beta)
    } else {
      mu
    }


    # ---- Preallocate zero matrices for alpha, eta, and impact ----
    zero_mat_for_alpha <- init_zero_matrix(object@dimens, ncol(beta))
    zero_mat_for_eta <- zero_mat_for_alpha  # same size
    zero_mat_for_impact <- zero_mat_for_alpha  # same size

    rambda_component_n <- lambda_component0


    # handling both basic and flexible models
    dresidual <- if(is.null(object@dresidual)) stats::dexp else object@dresidual
    presidual <- if(is.null(object@presidual)) stats::pexp else object@presidual
    qresidual <- object@qresidual

    ## main loop
    for (n in 2:size) {

      mu_n <- rmu_n; rambda_component_n_minus_1 <- rambda_component_n

      type_n <- type[n]

      decayed <- exp(-beta * inter_arrival[n])


      # lambda decayed due to time, impact due to mark is not added yet
      lambda_component_n <- rambda_component_n_minus_1 * decayed

      # integrate mu over inter_arrival w.r.t. time
      if ((is.matrix(mu_n) && nrow(mu_n) == object@dimens && ncol(mu_n) == 1) ||
          (is.vector(mu_n) && length(mu_n) == 1)) {

        # col num : dimens
        integrated_mu <- mu_n * inter_arrival[n]

      } else {

        # When mu is observed as piecewise constant between n-1 and n
        # Rows from 1 to object@dimen are mu values for each dimension.
        # The last row represents the inter-arrival where mu changes between n-1 and n.
        sum_mu_n <- colSums(mu_n[1:object@dimens, ])

        # col num : dimens
        integrated_mu <- sum_mu_n * mu_n[object@dimens + 1, ]

      }

      if (is.null(object@dresidual) && !infer){

        # basic model
        ## 1. sum of integrated_lambda_component and sum of mu * inter_arrival
        ## In the log-likelihood calculation, a negative sign is applied, so it's subtracted here.
        sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi - sum(rambda_component_n_minus_1 / beta * (1 - decayed)) - sum(integrated_mu)


      } else {

        # flexible model or when needed to infer
        for (k in 1:object@dimens){

          phi_k_wo_mu <- phi_wo_mu(tau = inter_arrival[n],
                                   rowSum_rambda_component_n_miuns_1 = sum(rambda_component_n_minus_1[k, ]),
                                   beta[k, 1])

          # partial residual for k-th dimension
          phi_value_k <- phi_k_wo_mu + integrated_mu[k]

          log_value <- suppressWarnings(if (k == type_n) log(dresidual( phi_value_k)) else log(1-presidual(phi_value_k)))
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
      if(is_lambda_necessary) lambda[n, ] <- as.vector(mu_n) + rowSums(lambda_component_n)

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

      if(infer){

        rambda_component[n, ] <- t(rambda_component_n)
        rambda[n, ] <- as.vector(rmu_n) + rowSums(rambda_component_n)
        rmus[n, ] <- rmu_n

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


    if(!infer) {
      # only log-likelihood
      c(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark)

    } else {

      return_object <- list(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark,
                            lambda = lambda, lambda_component = lambda_component,
                            rambda = rambda, rambda_component = rambda_component)

      if (is.null(object@dresidual)){

        # basic model

        for (k in 1:object@dimens){

          # the ordered sequence of indices associated with type k
          n_k <- which(type == k)
          res_process <- rep(0, length(n_k) - 1)

          for(j in seq_along(res_process)){
            # n_k[j]+1 : index of previous jump + 1
            # n_k[j+1] : index of current jump
            res_process[j] <- sum(phis[(n_k[j]+1):n_k[j+1], k])
          }

          return_object[[paste0("res_process", k)]] <- res_process

        }

      } else {

        # flexible model
        for (k in 1:object@dimens){

          # the ordered sequence of indices associated with type k
          n_k <- which(type == k)
          res_process <- rep(0, length(n_k) - 1)
          eres_process <- rep(0, length(n_k) - 1)

          for(j in seq_along(res_process)){


            # n_k[j]+1 : index of previous jump + 1
            # n_k[j+1] : index of current jump
            # F(phi(tau)) = F(epsilon)
            F_phi_taus <- presidual(phis[(n_k[j]+1):n_k[j+1], k])

            res_process[j] <- qresidual(1 - prod(1 - F_phi_taus))

            eres_process[j] <- - sum(log(1 - F_phi_taus))

          }

          return_object[[paste0("res_process", k)]] <- res_process
          return_object[[paste0("eres_process", k)]] <- eres_process

        }
      }

      return_object

    }

  }
)

# Helpers: Initialization for lambda matrices
initialize_lambda_matrix <- function(prefix, size, dimens, beta_cols, lambda0) {
  mat <- matrix(0, nrow = size, ncol = dimens * beta_cols)
  mat[1, ] <- t(lambda0)
  colnames(mat) <- paste0(prefix, rep(1:dimens, each = beta_cols), 1:beta_cols)
  return(mat)
}

initialize_lambda_total <- function(prefix, size, dimens, mu0, lambda0) {
  mat <- matrix(0, nrow = size, ncol = dimens)
  mat[1, ] <- as.vector(mu0) + rowSums(matrix(lambda0, nrow = dimens))
  colnames(mat) <- paste0(prefix, 1:dimens)
  return(mat)
}

# ---- Preallocate zero matrices for alpha, eta, and impact ----
init_zero_matrix <- function(dimens, beta_cols) {
  matrix(0, nrow = dimens, ncol = beta_cols)
}

