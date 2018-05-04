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
#' @param n the number of observations.
#'
#'
#' @export
setMethod(
  f="hsim",
  signature(object = "hspec"),
  definition = function(object, lambda0 = NULL, n = 100){

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
    mark_hawkes <- object@mark_hawkes
    mark_lambda <- object@mark_lambda
    impact <- object@impact

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda0, c, numeric(length = n - 1)), ncol = dimens^2)
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    lambda   <- matrix(sapply(mu + rowSums_lambda0, c, numeric(length = n - 1)), ncol = dimens)

    Ng <- matrix(numeric(length = dimens * n), ncol = dimens)
    N  <- matrix(numeric(length = dimens * n), ncol = dimens)

    event_idx <- numeric(length = n)
    inter_arrival <- numeric(length = n)
    mark <- numeric(length = n)

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
    colnames(N)  <- paste0("N", 1:dimens)
    colnames(Ng) <- paste0("Ng", 1:dimens)

    # Exact method
    for (i in 2:n) {

      # Generate candidate arrivals
      # arrival due to mu
      candidate_arrival <- stats::rexp(dimens, rate = mu)
      current_LAMBDA <- matrix(as.numeric(lambda_component[i-1, ]), nrow = dimens, byrow = TRUE)

      # arrival due to components

      matrixD <- 1 + beta * log(stats::runif(dimens^2)) / current_LAMBDA
      candidate_arrival <- cbind(candidate_arrival, -1 / beta * log(pmax(matrixD, 0)))

      # The minimum is inter arrival time
      inter_arrival[i] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[i], arr.ind = TRUE) #row and col

      event_idx[i] <- minIndex[1]  # row

      # lambda decayed due to time, impact due to mark is not added yet
      decayled_lambda <- current_LAMBDA * exp(-beta * inter_arrival[i])
      lambda_component[i, ] <- t(decayled_lambda)
      lambda[i, ] <- mu + rowSums(decayled_lambda)

      # generate a mark for Hawkes
      # if !is.null

      mark[i] <- mark_hawkes(n = 1, i = i, N = N, Ng = Ng,
                             lambda = lambda, lambda_component = lambda_component,
                             event_idx = event_idx)

      # shoud be a matrix
      realized_mark_lambda <- mark_lambda(n = 1, i = i, N = N, Ng = Ng,
                                          lambda = lambda, lambda_component = lambda_component,
                                          event_idx = event_idx)

      Ng[i, ] <- Ng[i-1, ]
      Ng[i, event_idx[i]] <- Ng[i-1, event_idx[i]] + 1
      N[i, ] <- N[i-1, ]
      N[i, event_idx[i]] <- N[i-1, event_idx[i]] + mark[i]

      # update lambda

      new_lambda <- decayled_lambda + impact(realized_mark_lambda)
      lambda_component[i, ] <- t(new_lambda)
      lambda[i, ] <- mu + rowSums(new_lambda)

      # if (dimens == 1) {
      #   Impact <- alpha * (1 + (mark[i] - 1) * ETA )
      # } else {
      #   Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      #   Impact[ , event_idx[i]] <- alpha[ , event_idx[i]] * (1 + (mark[i] - 1) * ETA[ , event_idx[i]])
      # }

      # new_LAMBDA = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      #
      # Impact is added.
      new_lambda <- decayled_lambda + Impact
      lambda_component[i, ] <- t(new_lambda)
      lambda[i, ] <- mu + rowSums(new_lambda)
    }


    realization <- list(object, inter_arrival, cumsum(inter_arrival), event_idx, mark, N, Ng, lambda, lambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "event_idx", "mark", "N", "Ng", "lambda", "lambda_component")
    class(realization) <- c("hreal")

    return(realization)
  }
)
