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
    mark_hawkes <- object@mark_hawkes
    impact <- object@impact

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = dimens^2)
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    lambda   <- matrix(sapply(mu + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)

    N <- matrix(numeric(length = dimens * size), ncol = dimens)
    Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)

    event_idx <- numeric(length = size)
    inter_arrival <- numeric(length = size)
    mark <- numeric(length = size)

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
    colnames(Nc)  <- paste0("Nc", 1:dimens)
    colnames(N) <- paste0("N", 1:dimens)

    # Exact method
    for (n in 2:size) {

      # Generate candidate arrivals
      # arrival due to mu
      candidate_arrival <- stats::rexp(dimens, rate = mu)
      current_LAMBDA <- matrix(as.numeric(lambda_component[n-1, ]), nrow = dimens, byrow = TRUE)

      # arrival due to components

      matrixD <- 1 + beta * log(stats::runif(dimens^2)) / current_LAMBDA
      candidate_arrival <- cbind(candidate_arrival, -1 / beta * log(pmax(matrixD, 0)))

      # The minimum is inter arrival time
      inter_arrival[n] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[n], arr.ind = TRUE) #row and col

      event_idx[n] <- minIndex[1]  # row

      # lambda decayed due to time, impact due to mark is not added yet
      decayed_lambda <- current_LAMBDA * exp(-beta * inter_arrival[n])
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu + rowSums(decayed_lambda)

      # generate a mark for Hawkes
      # This quantity is added to the counting process.

      if( !is.null(mark_hawkes) ){
        # mark may depends on other variables
        # mark[n] is a scalar
        mark[n] <- mark_hawkes(n = n, Nc = Nc, N = N,
                               lambda = lambda, lambda_component = lambda_component,
                               event_idx = event_idx)
      } else {
        mark[n] <- 1
      }
#
#       if( !is.null(mark_lambda) ){
#
#         realized_mark_lambda <- mark_lambda(n = n, Nc = Nc, N = N, mark,
#                                             lambda = lambda, lambda_component = lambda_component,
#                                             event_idx = event_idx)
#
#         if( !is.matrix(realized_mark_lambda) ){
#           stop("mark_lambda should return a matrix.")
#         }
#         if( dim(realized_mark_lambda)[1] != dimens | dim(realized_mark_lambda)[2] != dimens) {
#           stop("The dimension of mark_lambda should be checked.")
#         }
#
#
#       } else {
#         realized_mark_lambda <- matrix(rep(0, dimens^2), nrow=dimens)
#       }
      # shoud be a matrix

      N[n, ] <- N[n-1, ]
      N[n, event_idx[n]] <- N[n-1, event_idx[n]] + 1
      Nc[n, ] <- Nc[n-1, ]
      Nc[n, event_idx[n]] <- Nc[n-1, event_idx[n]] + mark[n]

      # update lambda

      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens^2), nrow = dimens)
      impact_alpha[ , event_idx[n]] <- alpha[ , event_idx[n]]

      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens^2), nrow = dimens)

        impact_res <- impact(n = n, Nc = Nc, N = N, mark = mark,
                             lambda = lambda, lambda_component = lambda_component,
                             event_idx = event_idx)

        impact_mark[ , event_idx[n]] <- impact_res[ , event_idx[n]]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }

      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}

      lambda_component[n, ] <- t(new_lambda)
      lambda[n, ] <- mu + rowSums(new_lambda)

      # if (dimens == 1) {
      #   Impact <- alpha * (1 + (mark[n] - 1) * ETA )
      # } else {
      #   Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      #   Impact[ , event_idx[n]] <- alpha[ , event_idx[n]] * (1 + (mark[n] - 1) * ETA[ , event_idx[n]])
      # }


      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      #
      # Impact is added.

    }

    realization <- list(object, inter_arrival, cumsum(inter_arrival), event_idx, mark, Nc, N, lambda, lambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "event_idx", "mark", "Nc", "N", "lambda", "lambda_component")
    class(realization) <- c("hreal")

    return(realization)
  }
)
