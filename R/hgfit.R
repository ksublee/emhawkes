#' @include hspec.R hmoment.R utilities.R
NULL

#' Infer Intensity Paths (Deprecated)
#'
#' Deprecated wrapper for [logLik()] with `infer = TRUE`. Use `logLik()` directly
#' to compute intensity paths and residual processes from a specified model.
#'
#' @param object \code{\link{hspec-class}}. This object includes the parameter values.
#' @param inter_arrival inter-arrival times of events. This includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param N Optional cumulative event-count matrix consistent with `type`.
#' @param Nc Optional cumulative mark-total matrix consistent with `type` and `mark`.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' A scalar is expanded to every kernel component. A matrix must match the
#' shape of `alpha` and `beta`.
#' @param N0 Initial event counts: a scalar or one value per event type.
#' @param Nc0 Initial mark totals: a scalar or one value per event type.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The inference list returned by `logLik(object, infer = TRUE)`.
#' It contains the likelihood, intensity matrices, and residual vectors.
#' @seealso [logLik()], [fitted_hspec()], [fit_initial_state()]
#'
#' @docType methods
#' @rdname infer_lambda
#' @export
#'
#' @examples
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' summary(res)
#' inferred <- logLik(h, inter_arrival = res$inter_arrival,
#'                    type = res$type, infer = TRUE)
#' head(inferred$lambda)
#'
setGeneric("infer_lambda", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL,
                            N = NULL, Nc = NULL,
                            lambda_component0 = NULL, N0 = NULL, Nc0 = NULL, ...) standardGeneric("infer_lambda"))

#' @rdname infer_lambda
setMethod(
  f="infer_lambda",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL,
           N = NULL, Nc = NULL,
           lambda_component0 = NULL, N0 = NULL, Nc0 = NULL, ...){

    warning("infer_lambda() is deprecated; use logLik(..., infer = TRUE) instead.", call. = FALSE)


    infered_res <- logLik(object, inter_arrival = inter_arrival, type = type,
                          mark = mark, N = N, Nc = Nc,  lambda_component0 = lambda_component0,
                          N0 = N0, Nc0 = Nc0,
                          infer = TRUE, ...)

    return(infered_res)


    # Similar to logLik function, need to integrate


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

    size <-  length(inter_arrival) #

    # check the argument lists in mu and impact
    mu_args <- impct_args <- c()
    if(is.function(mu) && length(formals(mu)) > 1) mu_args <- methods::formalArgs(mu)
    if(!is.null(object@impact)) impct_args <- methods::formalArgs(object@impact)

    # N and Nc are constructed only if they are needed.
    if("N" %in% c(impct_args, mu_args) && is.null(N)) N <- type_to_N(type, object@dimens, N0 = N0)
    if("Nc" %in% c(impct_args, mu_args) && is.null(Nc)) N <- type_mark_to_Nc(type, mark, object@dimens, Nc0 = Nc0)


    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                N = N, Nc = Nc,
                lambda = lambda, lambda_component = lambda_component,
                alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      mu0 <- mu
    }

    # default lambda_component0
    if(!is.null(object@type_col_map)){
      if(length(object@type_col_map) > 0 & is.null(lambda_component0)){
        stop("In this model, please provide lambda_component0.")
      }
    }


    if(!is.null(lambda_component0)){

      # If the dimensions of model and lambda_component0 do not match, lambda_component0 will be adjusted
      if (object@dimens * ncol(beta) > length(lambda_component0)){
        warning("The size of lambda_component0 does not match to the dimension of the model and is adjusted. \n
                lambda_component0 is now :")
        lambda_component0 <- rep(lambda_component0, object@dimens^2)
      }
      if (object@dimens * ncol(beta) < length(lambda_component0)){
        warning("The size of lambda_component0 does not match to the dimension of the model and is adjusted.\n
                lambda_component0 is now :")
        lambda_component0 <- lambda_component0[1:object@dimens^2]
      }

      lambda_component0 <- as.matrix(lambda_component0, nrow = object@dimens)

    } else {
      # default lambda_component0
      message("The initial values for intensity processes are not provided. Internally determined initial values are used for inference.\n")
      lambda_component0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu = mu, alpha = alpha, beta = beta)
    }


    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda_component0, c, numeric(length = size - 1)), ncol = object@dimens * ncol(beta))
    rowSums_lambda_component0 <- rowSums(matrix(lambda_component0, nrow=object@dimens))

    lambda <- matrix(sapply(mu0 + rowSums_lambda_component0, c, numeric(length = size - 1)), ncol = object@dimens)
    rambda <- lambda
    rambda_component <- lambda_component

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:object@dimens)
    indxM <- matrix(rep(1:ncol(beta), object@dimens), byrow = TRUE, nrow = object@dimens)
    colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:object@dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))

    colnames(rambda) <- paste0("rambda", 1:object@dimens)
    colnames(rambda_component) <- as.vector(t(outer(as.vector(outer("rambda", 1:object@dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))


    current_lambda <- lambda_component0

    # only piecewise constant mu is available, to be updated
    if (is.function(mu)){
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                  N = N, Nc = Nc,
                  lambda = lambda, lambda_component = lambda_component,
                  lambda_component_n = lambda_component_n,
                  alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }

    zero_mat_for_alpha <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize
    zero_mat_for_eta <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize
    zero_mat_for_impact <- matrix(rep(0, object@dimens * ncol(beta)), nrow = object@dimens)  # for fast computation, pre-initialize

    for (n in 2:size) {

      mu_n <- rmu_n
      type_n <- type[n]
      inter_arrival_n <- inter_arrival[n]

      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival_n)
      decayed_lambda <- lambda_component_n <- current_lambda * decayed

      # update lambda
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu_n + rowSums(decayed_lambda)


      ## impact

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

      current_lambda <- decayed_lambda + impact_alpha

      # 2. additional impact by eta

      if(!is.null(eta)) {

        impact_eta <- zero_mat_for_eta   # matrix
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        current_lambda <- current_lambda + impact_eta
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
        current_lambda <- current_lambda + impact_mark

      }


      # new mu, i.e., right continuous version of mu
      if (is.function(mu)){
        # mu is represeted by function
        rmu_n <- mu(n = n + 1, mark = mark, type = type, inter_arrival = inter_arrival,
                    N = N, Nc = Nc,
                    lambda_component_n = lambda_component_n,
                    alpha = alpha, beta = beta)
      } else{
        # mu is a matrix
        rmu_n <- mu
      }


      # update rambda
      # rambda_component = {"rambda11", "rambda12", ..., "rambda21", "rambda22", ...}
      rambda_component[n, ] <- t(current_lambda)
      rambda[n, ] <- rmu_n + rowSums(current_lambda)

      #current_rambda_component <- new_lambda
    }

    realization <- list(object, inter_arrival, cumsum(inter_arrival), type, mark,
                        N, Nc, lambda, lambda_component, rambda, rambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "type", "mark",
                            "N", "Nc", "lambda", "lambda_component", "rambda", "rambda_component")
    class(realization) <- c("hreal")

    realization
  }

)

integrate_rambda_component <- function(inter_arrival, rambda_componet, beta, dimens){

  size <- length(inter_arrival)

  integrated_rambda_component <- matrix(rep(0, nrow(rambda_componet)*ncol(rambda_componet)),
                                            nrow=nrow(rambda_componet))

  current_lambda <- matrix(rambda_componet[1,], nrow = dimens, byrow=TRUE)


  for (n in seq_len(size - 1L) + 1L) {

    decayed <- exp(-beta * inter_arrival[n])
    integrated_rambda_component[n,] <- t(integrate_exp_kernel(current_lambda, beta, inter_arrival[n]))
    new_lambda <- rambda_componet[n,]
    current_lambda <- matrix(new_lambda, nrow=dimens, byrow=TRUE)

  }

  integrated_rambda_component
}


integrate_rambda <- function(inter_arrival, rambda_component, mu, beta, dimens,
                             alpha = NULL,
                             type = NULL, mark = NULL,
                             N = NULL, Nc = NULL, lambda = NULL,
                             lambda_component = NULL, lambda_component_n = NULL,
                             lambda_component0 = NULL, N0 = NULL){

  size <- length(inter_arrival)


  if (is.function(mu)){
    mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
              N = N, Nc = Nc,
              lambda_component_n = lambda_component_n,
              alpha = alpha, beta = beta)
  } else {
    mu0 <- mu
  }

  integrated_rambda <- matrix(rep(0, length(rambda_component)/dimens),
                                        nrow=nrow(rambda_component))

  current_lambda <- matrix(rambda_component[1,], nrow = dimens, byrow=TRUE)
  size <- length(inter_arrival)

  # only piecewise constant mu is available, to be updated
  if (is.function(mu)){
    rmu_n <- mu(n = 2, type = type, inter_arrival = inter_arrival, mark = mark,
                N = N, Nc = Nc,
                lambda = lambda, lambda_component = lambda_component,
                alpha = alpha, beta = beta)
  } else{
    rmu_n <- mu
  }


  for (n in seq_len(size - 1L) + 1L) {

    mu_n <- rmu_n

    decayed <- exp(-beta * inter_arrival[n])

    integrated_rambda[n, ] <- rowSums(integrate_exp_kernel(current_lambda, beta, inter_arrival[n])) +
      mu_n * inter_arrival[n]

    new_lambda <- rambda_component[n,]
    current_lambda <- matrix(new_lambda, nrow=dimens, byrow=TRUE)

    # new mu, i.e., right continuous version of mu
    if (is.function(mu)){
      # mu is represented by function
      rmu_n <- mu(n = n + 1, mark = mark, type = type, inter_arrival = inter_arrival,
                  N = N, Nc = Nc,
                  lambda_component_n = lambda_component_n,
                  alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      rmu_n <- mu
    }
  }

  integrated_rambda

}

#' Compute a Residual Process (Deprecated)
#'
#' Legacy helper that integrates a supplied intensity path for one event type.
#' Use [logLik()] with `infer = TRUE` to compute intensities and residuals together.
#' For a correctly specified Hawkes model, residual intervals have a unit-exponential distribution.
#'
#' @param component The component of type to get the residual process.
#' @param type A vector of types distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param inter_arrival Inter-arrival times of events. This includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param rambda_component Right continuous version of lambda process.
#' @param mu Numeric value or matrix or function. If numeric, automatically converted to matrix.
#' @param beta Numeric value or matrix or function. If numeric, automatically converted to matrix, exponential decay.
#' @param dimens Dimension of the model. If omitted, set to be the length of \code{mu}.
#' @param mark A vector of realized mark (jump) sizes. Start with zero.
#' @param N A matrix of counting processes.
#' @param Nc A matrix of counting processes weighted by mark.
#' @param lambda_component0 The initial values of lambda component. Must have the same dimensional matrix with \code{hspec}.
#' @param N0 The initial value of N
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A numeric vector of residual intervals for the selected event type.
#' @details
#' This function is deprecated. The replacement `logLik()` returns residuals in
#' `res_process1`, ..., `res_processD` and handles types with fewer than two events
#' by returning empty vectors. For fitted parameters, first use [fitted_hspec()]
#' and preserve initialization with [fit_initial_state()].
#' @seealso [logLik()], [fitted_hspec()], [fit_initial_state()]
#'
#' @examples
#'
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=1000)
#' inferred <- logLik(h, inter_arrival = res$inter_arrival,
#'                    type = res$type, infer = TRUE)
#' rp <- inferred$res_process1
#'
#' @export
residual_process <- function(component, inter_arrival, type, rambda_component, mu, beta, dimens=NULL,
                             mark = NULL,
                             N = NULL, Nc = NULL,
                             lambda_component0 = NULL, N0 = NULL, ...){

  warning("residual_process() is deprecated; use logLik(..., infer = TRUE) and its res_process fields instead.", call. = FALSE)


  additional_argument <- list(...)
  if ("lambda0" %in% names(additional_argument)) {

    warning("'lambda0' is deprecated; use 'lambda_component0' instead.", call. = FALSE)

    lambda_component0 <- additional_argument[["lambda0"]]

  }

  if (is.null(dimens)){
    if (is.function(mu)){
      stop("argument dimens is missing.")

    } else{
      dimens <- length(mu)
    }
  }
  integrated_rambda <- integrate_rambda(inter_arrival, rambda_component, mu, beta, dimens,
                                        type = type, mark = mark,
                                        N = N, Nc = Nc,
                                        lambda_component0 = lambda_component0, N0 = N0)

  row_idx <- which(type == component)
  res_process <- rep(0, (length(row_idx)-1))
  for( i in 1:length(res_process)){

    res_process[i] <- sum(integrated_rambda[(row_idx[i] + 1):row_idx[i+1], component])
  }


  res_process

}
