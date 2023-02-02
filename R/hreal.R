#' Realization of Hawkes process
#'
#' @description
#' `hreal` is the list of the following:
#'
#' * `hspec` : S4 object \code{\link{hspec-class}} that specifies the parameter values.
#' * `inter_arrival` : the time between two consecutive events.
#' * `arrival` : cumulative sum of `inter_arrival`.
#' * `type` : integer, the type of event.
#' * `mark` : the size of mark, an additional information associated with event.
#' * `N` : counting process that counts the number of events.
#' * `Nc` : counting process that counts the number of events weighted by mark.
#' * `lambda` : intensity process, left-continuous version.
#' * `lambda_component` : the component of intensity process with `mu` not included.
#' * `rambda` : intensity process, right-continuous version.
#' * `rambda_component` : the right-continuous version of `lambda_component`.
#'
#'
#'
#' @name hreal

#' @description
#' Print functions for `hreal` are provided.
#'
#'
#' @param x S3-object of `hreal`.
#' @param n Number of rows to display.
#' @param object S3-object of `hreal`.
#' @param ... Further arguments passed to or from other methods.
#'
#'
#' @rdname hreal
#'
#' @export
print.hreal <- function(x, n=20, ...){
  options(digits=4)
  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
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

  cat("------------------------------------------\n")
  options(digits=7)
}


#'
#'
#' @rdname hreal
#' @export
summary.hreal <- function(object, n=20, ...){

  options(digits=5)
  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
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

  cat("------------------------------------------\n")
  options(digits=7)
}


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


