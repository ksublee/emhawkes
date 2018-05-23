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

unique_param_names <- function(M, notation){

  reference <- character(length(M))

  if (ncol(M) == 1){
    k <- 1
    for (i in 1:nrow(M)){
      if (reference[k] == ""){
        if(M[i] != 0)
          reference[which(M == M[i])] <- paste0(notation, toString(i))
      }
      k <- k + 1
    }
  } else {
    k <- 1
    for  (i in 1:nrow(M)){
      for (j in 1:ncol(M)) {
        if (reference[k] == ""){
          if(M[i,j] != 0)
            reference[which(t(M) == M[i,j])] <- paste0(notation, toString(i), toString(j))
        }
        k <- k + 1
      }
    }
  }
  reference

}

as.param <- function(M, prefix){

  unique_M <- M[full_names(M, prefix) == attr(M, "param.names") &
                  !is.na(attr(M, "param.names"))]
  names(unique_M) <- attr(M, "param.names")[full_names(M, prefix) == attr(M, "param.names") &
                                              !is.na(attr(M, "param.names"))]

  unique_M
}

param.names2 <- function(M, notation, sep=","){

  my_paste <- function(...){
    paste(..., sep=sep)
  }

  m <- as.vector(M)

  if (ncol(M) == 1){
    names(m) <- as.vector(outer(notation, 1:length(M), FUN = paste0))
  } else {
    names(m) <- as.vector(outer(as.vector(outer(notation, 1:4, FUN = paste0)),
                                1:4, FUN=my_paste))
  }


  m[!duplicated(m)]



}

full_names <- function(M, notation, sep=","){

  my_paste <- function(...){
    paste(..., sep=sep)
  }


  if (ncol(M) == 1){
    names_M <- outer(notation, 1:length(M), FUN = paste0)
  } else {
    names_M <- outer(as.vector(outer(notation, 1:nrow(M), FUN = paste0)),
                     1:nrow(M), FUN=my_paste)
  }

  names_M
}

param.names <- function(M, name_M, prefix="", sep=","){

  my_paste <- function(...){
    paste(..., sep=sep)
  }

  full_name <- as.vector(t(full_names(M, prefix, sep=sep)))

  names(full_name) <- as.vector(t(name_M))
  matrix(full_name[names(full_name)], nrow=nrow(M), byrow=T)

}


#Thanks to https://www.r-bloggers.com/hijacking-r-functions-changing-default-arguments/
hijack <- function (FUN, ...) {
  .FUN <- FUN
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  .FUN
}
