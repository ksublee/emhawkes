#' Realization of Hawkes Process
#'
#' @description
#' `hreal` is an S3 list returned by [hsim()], containing the following components:
#'
#' * `hspec`: An S4 object of class \code{\link{hspec-class}} that specifies the parameter values.
#' * `inter_arrival`: The time intervals between consecutive events.
#' * `arrival`: The cumulative sum of `inter_arrival` times.
#' * `type`: An event-type vector. Values 1, ..., `dimens` identify intensity rows.
#' * `mark`: A vector of event marks.
#' * `N`: A cumulative event-count matrix, with one column per event type.
#' * `Nc`: A cumulative mark-total matrix, with one column per event type.
#' * `lambda`: The left-continuous intensity process.
#' * `lambda_component`: The component of the intensity process, \eqn{\lambda_{ij}}, that excludes `mu`.
#' * `rambda`: The right-continuous intensity process.
#' * `rambda_component`: The right-continuous version of `lambda_component`.
#'
#' @details
#' The first row is the initial state, with zero `inter_arrival`, `type`, and `mark`.
#' Initial counts and mark totals are stored in the first rows of `N` and `Nc`.
#' Subsequent counting-process rows contain values after each event.
#' Kernel-component columns are ordered by intensity row, then kernel column.
#' [logLik()] with `infer = TRUE` returns an inference list, with likelihood and
#' residual fields, rather than an `hreal` simulation result.
#' @seealso [hsim()], [logLik()]
#' @name hreal
#'
#' @description
#' Functions for printing `hreal` objects are provided.
#'
#' @param x An S3 object of class `hreal`.
#' @param n The number of rows to display: a nonnegative integer or `Inf`.
#' @param object An S3 object of class `hreal`.
#' @param ... Additional arguments passed to or from other methods.
#'
#'
#' @rdname hreal
#'
#' @export
print.hreal <- function(x, n=20, ...){
  rows <- hreal_display_rows(n, length(x$arrival))
  old_options <- options(digits=4)
  on.exit(options(old_options), add = TRUE)
  cat("-------------------------------------------------------\n")
  cat("Simulation result of exponential (marked) Hawkes model.\n")
  print(x$hspec)

  cat("Realized path :\n")
  mtrx <- as.matrix(x)
  dimens <- x$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)
  name_lambda_component <- colnames(x$lambda_component)

  len <- length(rows)

  if(is.null(x$mark)) {
    print(mtrx[rows, c("arrival", name_N, name_lambda, name_lambda_component), drop = FALSE])
  } else {
    print(mtrx[rows, c("arrival", name_N, "mark", name_lambda, name_lambda_component), drop = FALSE])
  }

  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("-------------------------------------------------------\n")
  invisible(x)
}


#'
#'
#' @rdname hreal
#' @export
summary.hreal <- function(object, n=20, ...){

  rows <- hreal_display_rows(n, length(object$arrival))
  old_options <- options(digits=5)
  on.exit(options(old_options), add = TRUE)
  cat("-------------------------------------------------------\n")
  cat("Simulation result of exponential (marked) Hawkes model.\n")
  cat("Realized path :\n")
  mtrx <- as.matrix(object)
  dimens <- object$hspec@dimens
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)

  len <- length(rows)

  if(is.null(object$mark)) {
    print(mtrx[rows, c("arrival", name_N, name_lambda), drop = FALSE])
  } else {
    print(mtrx[rows, c("arrival", name_N, "mark", name_lambda), drop = FALSE])
  }

  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("-------------------------------------------------------\n")
  invisible(object)
}

hreal_display_rows <- function(n, size) {
  if (!is.numeric(n) || is.complex(n) || length(n) != 1L ||
      is.na(n) || n < 0 || n != floor(n)) {
    stop("n must be a nonnegative integer or Inf.")
  }
  seq_len(min(n, size))
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


