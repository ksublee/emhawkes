#' Generics for hreal
#'
#' Generic functions list for hreal:
#'
#' @name hreal
NULL

#'
#' Print function for hreal
#'
#' Print the realization of the Haweks model.
#'
#' @param x S3-object of hreal.
#' @param object S3-object of hreal.
#' @param n number of rows to diplay.
#' @param ... further arguments passed to or from other methods.
#' @rdname hreal
#' @export
print.hreal <- function(x, n=20, ...){
  options(digits=4)
  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
  print(x$hspec)

  cat("Realized path (with right continuous representation):\n")
  mtrx <- as.matrix(x)
  dimens <- x$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)
  name_lambda_component <- colnames(x$lambda_component)

  len <- min(n, length(mtrx[,"arrival"]))

  print(mtrx[1:len, c("arrival", name_N, name_lambda, name_lambda_component)])
  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("------------------------------------------\n")
  options(digits=7)
}

#' Summary function for hreal
#'
#' Print the summary of the Hawkes process realization.
#'
#' @rdname hreal
#' @export
summary.hreal <- function(object, n=20, ...){

  options(digits=5)
  cat("------------------------------------`------\n")
  cat("Simulation result of marked Hawkes model.\n")
  cat("Realized path (with right continuous representation):\n")
  mtrx <- as.matrix(object)
  dimens <- object$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)

  len <- min(n, length(mtrx[,"arrival"]))

  print(mtrx[1:len, c("arrival", name_N, name_lambda)])
  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("------------------------------------------\n")
  options(digits=7)
}


#' Matrix represetation of hreal
#'
#' Matrix like ouput of the realization of Hawkes model.
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


