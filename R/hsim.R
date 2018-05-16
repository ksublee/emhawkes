#' @include hspec.R hmoment.R

setGeneric("hsim", function(object, ...) standardGeneric("hsim"))

#' Simulate a (marked) Hawkes process
#'
#' The method simulate marked Hawkes processes.
#' The object \code{\link{hspec-class}} contains the parameter values such as \code{mu}, \code{alpha}, \code{beta}.
#' The mark (jump) structure may or may not be included.
#' It returns an object of class '\code{\link{hreal-class}}.
#'
#' @param object \code{\link{hspec-class}}. This object includes the parameter values.
#' @param lambda0 the starting values of lambda (intensity process). numeric or matrix.
#' @param size the number of observations.
#'
#'
#' @export
setMethod(
  f="hsim",
  signature(object = "hspec"),
  definition = function(object, lambda0 = NULL, size = 100){

    # dimension of Hawkes process
    dimens <- length(object@mu)

    if(!is.null(lambda0)){

      # If the dimensions of model and lambda0 do not match, lambda0 will be adjusted
      if (dimens^2 > length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted. \n
                lambda0 is now :")
        lambda0 <- rep(lambda0, dimens^2)
      }
      if (dimens^2 < length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted.\n
                lambda0 is now :")
        lambda0 <- lambda0[1:dimens^2]
      }

      lambda0 <- as.matrix(lambda0, nrow = dimens)

    } else {
      # default lambda0
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object)
    }

    # parameter setting
    mu <- object@mu
    alpha <- object@alpha
    beta <- object@beta
    rmark <- object@rmark
    impact <- object@impact

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = dimens^2)
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    lambda <- matrix(sapply(mu + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)

    rambda <- lambda
    rambda_component <- lambda_component

    N <- matrix(numeric(length = dimens * size), ncol = dimens)
    Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)

    type <- numeric(length = size)
    inter_arrival <- numeric(length = size)
    mark <- numeric(length = size)

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))

    colnames(rambda) <- paste0("rambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(rambda_component) <- paste0("rambda", indxM, t(indxM))

    colnames(Nc)  <- paste0("Nc", 1:dimens)
    colnames(N) <- paste0("N", 1:dimens)



    # Exact method
    current_lambda <- lambda0
    for (n in 2:size) {

      # Generate candidate arrivals
      # arrival due to mu
      candidate_arrival <- stats::rexp(dimens, rate = mu)
      #current_LAMBDA <- matrix(as.numeric(lambda_component[n-1, ]), nrow = dimens, byrow = TRUE)

      # arrival due to components

      matrixD <- 1 + beta * log(stats::runif(dimens^2)) / current_lambda
      candidate_arrival <- cbind(candidate_arrival, -1 / beta * log(pmax(matrixD, 0)))

      # The minimum is inter arrival time
      inter_arrival[n] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[n], arr.ind = TRUE) #row and col

      type[n] <- minIndex[1]  # row

      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival[n])
      decayed_lambda <- current_lambda * decayed

      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu + rowSums(decayed_lambda)

      # generate a mark for Hawkes
      # This quantity is added to the counting process.
      N[n, ] <- N[n-1, ]
      N[n, type[n]] <- N[n-1, type[n]] + 1

      if( !is.null(rmark) ){
        # mark may depends on other variables
        # mark[n] is a scalar
        mark[n] <- rmark(n = n, Nc = Nc, N = N,
                               lambda = lambda, lambda_component = lambda_component,
                               type = type)
      } else {
        mark[n] <- 1
      }

      Nc[n, ] <- Nc[n-1, ]
      Nc[n, type[n]] <- Nc[n-1, type[n]] + mark[n]

      # update lambda

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

        impact_mark[ , type[n]] <- impact_res[ , type[n]]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }


      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      rambda_component[n, ] <- t(new_lambda)
      rambda[n, ] <- mu + rowSums(new_lambda)

      current_lambda <- new_lambda
    }

    realization <- list(object, inter_arrival, cumsum(inter_arrival), type, mark, N, Nc, lambda, lambda_component, rambda, rambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "type", "mark", "N", "Nc", "lambda", "lambda_component", "rambda", "rambda_component")
    class(realization) <- c("hreal")

    return(realization)
  }
)
