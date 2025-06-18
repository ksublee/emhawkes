#' Realization of Hawkes Process
#'
#' @description
#' `hreal` is a list containing the following components:
#'
#' * `hspec`: An S4 object of class \code{\link{hspec-class}} that specifies the parameter values.
#' * `inter_arrival`: The time intervals between consecutive events.
#' * `arrival`: The cumulative sum of `inter_arrival` times.
#' * `type`: An integer representing the type of event.
#' * `mark`: The size of the mark, providing additional information associated with the event.
#' * `N`: A counting process that tracks the number of events.
#' * `Nc`: A counting process that tracks the number of events, weighted by mark.
#' * `lambda`: The left-continuous intensity process.
#' * `lambda_component`: The component of the intensity process, \eqn{\lambda_{ij}}, that excludes `mu`.
#' * `rambda`: The right-continuous intensity process.
#' * `rambda_component`: The right-continuous version of `lambda_component`.
#'
#' @name hreal
#'
#' @description
#' Functions for printing `hreal` objects are provided.
#'
#' @param x An S3 object of class `hreal`.
#' @param n The number of rows to display.
#' @param object An S3 object of class `hreal`.
#' @param ... Additional arguments passed to or from other methods.
#'
#'
#' @rdname hreal
#'
#' @export
print.hreal <- function(x, n=20, ...){
  options(digits=4)
  cat("-------------------------------------------------------\n")
  cat("Simulation result of exponential (marked) Hawkes model.\n")
  print(x$hspec)

  cat("Realized path :\n")
  mtrx <- as.matrix(x)
  dimens <- x$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)
  name_lambda_component <- colnames(x$lambda_component)

  len <- min(n, length(mtrx[,"arrival"]))

  if(is.null(x$mark)) {
    print(mtrx[1:len, c("arrival", name_N, name_lambda, name_lambda_component)])
  } else {
    print(mtrx[1:len, c("arrival", name_N, "mark", name_lambda, name_lambda_component)])
  }

  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("-------------------------------------------------------\n")
  options(digits=7)
}


#'
#'
#' @rdname hreal
#' @export
summary.hreal <- function(object, n=20, ...){

  options(digits=5)
  cat("-------------------------------------------------------\n")
  cat("Simulation result of exponential (marked) Hawkes model.\n")
  cat("Realized path :\n")
  mtrx <- as.matrix(object)
  dimens <- object$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)

  len <- min(n, length(mtrx[,"arrival"]))

  if(is.null(object$mark)) {
    print(mtrx[1:len, c("arrival", name_N, name_lambda)])
  } else {
    print(mtrx[1:len, c("arrival", name_N, "mark", name_lambda)])
  }

  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("-------------------------------------------------------\n")
  options(digits=7)
}

#'
#'
#' @rdname hreal
#' @export
as.matrix.hreal <- function(x, ...){

  mtrx <- numeric()
  for (i in 2:length(x)){
    mtrx <- cbind(mtrx, x[[i]])
    if(is.vector(x[[i]])){
      colnames(mtrx)[i-1] <- names(x)[i]
    }
  }
  mtrx
}


