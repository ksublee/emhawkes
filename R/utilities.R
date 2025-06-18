
#convert type to N and Nc
type_to_N <- function(type, dimens, N0 = NULL){

  N <- matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(N) <- paste0("N", 1:dimens)

  if (!is.null(N0)) N[1,] <- N0

  N[cbind(2:length(type), type[2:length(type)])] <- 1
  N <- apply(N, 2, cumsum)

  N

}

type_mark_to_Nc <- function(type, mark, dimens, Nc0 = NULL){

  Nc <-  matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(Nc) <- paste0("Nc", 1:dimens)

  if (!is.null(Nc0)) Nc[1,] <- Nc0

  Nc[cbind(2:length(type), type[2:length(type)])] <- mark[2:length(mark)]
  Nc <- apply(Nc, 2, cumsum)

  Nc

}

# Unique naming coefficients from matrix
#
# With given matrix, look up the elements of the matrix
# and assign the same name to the elements with the same value.
# The name is based on \code{notation} and location of the element.
#
# @param M a matrix
# @param notation a string for name
#
# @return to covert to matrix, use byrow=TRUE
look_up_mtrx <- function(M, notation){
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
          reference[which(t(M) == M[i,j])] <- paste(paste0(notation, toString(i)), toString(j), sep=".")
        k <- k + 1
      }
    }
  }
  matrix(reference, nrow=nrow(M), byrow=TRUE)
}


full_names <- function(M, notation, sep="."){

  my_paste <- function(...){
    paste(..., sep=sep)
  }


  if (ncol(M) == 1){
    names_M <- outer(notation, 1:length(M), FUN = paste0)
  } else {
    names_M <- outer(as.vector(outer(notation, 1:nrow(M), FUN = paste0)),
                     1:ncol(M), FUN=my_paste)
  }

  names_M
}

# find unique vector
as.unique.vector <- function(M, notation, sep="."){

  my_paste <- function(...){
    paste(..., sep=sep)
  }

  m <- as.vector(t(M))

  if (ncol(M) == 1){
    names(m) <- as.vector(outer(notation, 1:length(M), FUN = paste0))
  } else {
    names(m) <- as.vector(t(outer(as.vector(outer(notation, 1:nrow(M), FUN = paste0)),
                                  1:ncol(M), FUN=my_paste)))
  }

  m[!duplicated(m)]


}


as.param <- function(M, prefix = "", reduced){

  if (is.null(M)) return(NULL)

  if (is.function(M)){
    # when the argument M is a function, extract the parameter from 'param' argument of M
    prs <- eval(formals(M)$param, envir = environment(M))

  } else{
    if(reduced){
      prs <- as.unique.vector(M, prefix)
    } else {
      prs <- as.vector(M)
      names(prs) <- as.vector(full_names(M, prefix))
    }
  }
  prs
}



#automatic parameter setting for mu, alpha, beta, eta with object
setting <- function(object) {

  list(
    mu    = eval_param(object@mu),
    alpha = eval_param(object@alpha),
    beta  = eval_param(object@beta),
    eta   = eval_param(object@eta)
  )

}


# Thanks to https://www.r-bloggers.com/hijacking-r-functions-changing-default-arguments/
# Copy function with changing default argument
hijack <- function (FUN, ...) {

  if(is.null(FUN)) return(NULL)

  .FUN <- FUN
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  .FUN
}


eval_param <- function(p) {
  if (is.function(p) && length(formals(p)) == 1)  evalf(p)
  else p
}


#function evaluation with default parameters
evalf <- function(FUN) {
  # Get the environment in which the function was defined
  fun_env <- environment(FUN)

  # Evaluate default parameters in that environment
  default_args <- lapply(formals(FUN), function(x) eval(x, envir = fun_env))

  # Execute the function with the default parameters
  result <- do.call(FUN, default_args)
  result

}



initialize_slot <- function(object_slot, pr_param, param_name = "", dimens = 0, reduced = TRUE) {

  if(is.null(object_slot)) return(NULL)

  if (is.function(object_slot)){

    slot0 <- hijack(object_slot, param = pr_param)

  } else{

    if(reduced){

      slot0 <- matrix(pr_param[look_up_mtrx(object_slot, param_name)], nrow=dimens)

    } else {

      slot0 <- matrix(pr_param, nrow=dimens)

    }
  }
  slot0
}
