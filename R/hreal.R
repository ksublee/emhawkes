#' Generic for hreal
#'
#' Generic functions list for hreal
#'
#' @name hreal
NULL

#'
#' Print function for hreal
#'
#' Print the summary of the realization of the Haweks model.
#'
#' @param res S3-object of hreal
#' @param n number of rows to diplay
#' @rdname hreal
#' @export
print.hreal <- function(res, n=20, ...){
  options(digits=4)
  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
  print(res$hspec)

  cat("Realized path (with right continuous representation):\n")
  mtrx <- as.matrix(res)
  dimens <- res$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)
  name_lambda_component <- colnames(res$lambda_component)

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

#' Summary function fo mhreal
#'
#' This function presents the summary of the Hawkes realization.
#'
#' @rdname hreal
#' @export
summary.hreal <- function(res, n=20){

  options(digits=5)
  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
  cat("Realized path (with right continuous representation):\n")
  mtrx <- as.matrix(res)
  dimens <- res$hspec@dimens
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
#' The realization of Hawkes model is represented by matrix like ouput.
#'
#' @rdname hreal
#' @export
as.matrix.hreal <- function(res){

  mtrx <- numeric()
  for (i in 2:length(res)){
    mtrx <- cbind(mtrx, res[[i]])
    if(is.vector(res[[i]])){
      colnames(mtrx)[i-1] <- names(res)[i]
    }
  }
  mtrx
}


