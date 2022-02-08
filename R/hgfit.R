#' @include hspec.R hmoment.R utilities.R
NULL

#' Infer lambda process with given Hawkes model and realized path
#'
#' This method compute the inferred lambda process and returns it as \code{hreal} form.
#' If we have realized path of Hawkes process and its parameter value, then we can compute the inferred lambda processes.
#' Similarly with other method such as \code{hfit}, the input arguments are \code{inter_arrival}, \code{type}, \code{mark},
#' or equivalently, \code{N} and \code{Nc}.
#'
#' @param object \code{\link{hspec-class}}. This object includes the parameter values.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param N Hawkes process. if not provided, then generate using inter_arrival and type.
#' @param Nc mark accumulated Hawkes process. if not provided, then generate using inter_arrival, type and mark.
#' @param lambda0 the initial values of lambda component. Must have the same dimensional matrix (n by n) with hspec.
#' @param N0 the initial values of N.
#'
#' @return hreal S3-object, the Haweks model with inferred intensity, lambda
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
#' res2 <- infer_lambda(h, res$inter_arrival, res$type)
#'
setGeneric("infer_lambda", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL,
                            N = NULL, Nc = NULL,
                            lambda0 = NULL, N0 = NULL) standardGeneric("infer_lambda"))
#'
#' @rdname infer_lambda
#'
setMethod(
  f="infer_lambda",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL,
           N = NULL, Nc = NULL,
           lambda0 = NULL, N0 = NULL){

    # Similar to logLik function, need to integrate


    #parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens
    size <-  length(inter_arrival)

    if(is.null(N) | is.null(Nc)){
      temp <- type_to_N(type, mark, dimens)
      if(is.null(N)) N <- temp[[1]]
      if(is.null(Nc)) Nc <- temp[[2]]
    }

    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                N = N, Nc = Nc,
                lambda_component_n = lambda_component_n,
                alpha = alpha, beta = beta)
    } else {
      mu0 <- mu
    }

    # default lambda0
    if(!is.null(object@type_col_map)){
      if(length(object@type_col_map) > 0 & is.null(lambda0)){
        stop("In this model, please provide lambda0.")
      }
    }


    if(!is.null(lambda0)){

      # If the dimensions of model and lambda0 do not match, lambda0 will be adjusted
      if (dimens * ncol(beta) > length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted. \n
                lambda0 is now :")
        lambda0 <- rep(lambda0, dimens^2)
      }
      if (dimens * ncol(beta) < length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted.\n
                lambda0 is now :")
        lambda0 <- lambda0[1:dimens^2]
      }

      lambda0 <- as.matrix(lambda0, nrow = dimens)

    } else {
      # default lambda0
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu = mu, alpha = alpha, beta = beta)
    }


    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = dimens * ncol(beta))
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    lambda <- matrix(sapply(mu0 + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)

    rambda <- lambda
    rambda_component <- lambda_component

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:ncol(beta), dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))

    colnames(rambda) <- paste0("rambda", 1:dimens)
    colnames(rambda_component) <- as.vector(t(outer(as.vector(outer("rambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))


    current_lambda <- lambda0

    # only piecewise constant mu is available, to be updated
    if (is.function(mu)){
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                  N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                  lambda_component_n = lambda_component_n,
                  alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }

    zero_mat_for_alpha <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize
    zero_mat_for_eta <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize
    zero_mat_for_impact <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)  # for fast computation, pre-initialize

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
      if(!is.null(impact)){

        impact_mark <- zero_mat_for_impact
        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
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


  for (n in 2:size) {

    decayed <- exp(-beta * inter_arrival[n])
    integrated_rambda_component[n,] <- t(current_lambda / beta * ( 1 - decayed ))
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
                             lambda0 = NULL, N0 = NULL){

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

    integrated_rambda[n, ] <- rowSums(current_lambda / beta * ( 1 - decayed )) +
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

#' Compute residual process
#'
#' Using random time change, this function compute the residual process, which is the inter-arrival time of a standard Poisson process.
#' Therefore, the return values should follow the exponential distribution with rate 1, if model and rambda are correctly specified.
#'
#' @param component the component of type to get the residual process
#' @param type a vector of types. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param rambda_component right continuous version of lambda process
#' @param mu numeric value or matrix or function, if numeric, automatically converted to matrix
#' @param beta numeric value or matrix or function, if numeric, automatically converted to matrix, exponential decay
#' @param dimens dimension of the model. if omitted, set to be the length of \code{mu}.
#' @param mark a vector of realized mark (jump) sizes. Start with zero.
#' @param N a matrix of counting processes
#' @param Nc a matrix of cumulated counting processes
#' @param lambda0 the initial values of lambda component. Must have the same dimensional matrix with \code{hspec}.
#' @param N0 the initial value of N
#'
#' @examples
#'
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=1000)
#' rp <- residual_process(1, res$type, res$inter_arrival, res$rambda_component, mu, beta)
#' p <- ppoints(100)
#' q <- quantile(rp,p=p)
#' plot(qexp(p), q, xlab="Theoretical Quantiles",ylab="Sample Quantiles")
#' qqline(q, distribution=qexp,col="blue", lty=2)
#'
#' @export
residual_process <- function(component, type, inter_arrival, rambda_component, mu, beta, dimens=NULL,
                             mark = NULL,
                             N = NULL, Nc = NULL,
                             lambda0 = NULL, N0 = NULL){

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
                                        lambda0 = lambda0, N0 = N0)

  row_idx <- which(type == component)
  res_process <- rep(0, (length(row_idx)-1))
  for( i in 1:length(res_process)){

    res_process[i] <- sum(integrated_rambda[(row_idx[i] + 1):row_idx[i+1], component])
  }


  res_process

}
