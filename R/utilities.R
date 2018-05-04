#' Unique naming coefficients from matrix
#'
#' With given matrix, look up the elements of the matrix
#' and assign the same name to the elements with the same value.
#' The name is based on \code{notation} and location of the element.
#'
#' @param M a square matrix
#' @param notation a string for name
#'
#' @return to covert to matrix, use byrow=TRUE
#'
#' @examples
#' A <- matrix(c(1,2,3,1), nrow = 2, byrow=TRUE)
#' name_unique_coef_mtrx(A, "alpha")
#'
#' @export
name_unique_coef_mtrx <- function(M, notation){
  reference <- character(length(M))

  if (ncol(M) == 1){
    k <- 1
    for (i in 1:nrow(M)){
      if (reference[k] == "")
        reference[which(M == M[i])] <- paste0(notation, toString(i))
      k <- k + 1
    }
  } else {
    k <- 1
    for  (i in 1:nrow(M)){
      for (j in 1:ncol(M)) {
        if (reference[k] == "")
          reference[which(t(M) == M[i,j])] <- paste0(notation, toString(i), toString(j))
        k <- k + 1
      }
    }
  }
  reference
}
