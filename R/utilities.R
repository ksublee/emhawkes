
#convert type to N and Nc
type_to_N <- function(type, dimens, N0 = NULL){

  N <- matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(N) <- paste0("N", 1:dimens)

  if (!is.null(N0)) N[1,] <- N0

  events <- seq_len(max(0L, length(type) - 1L)) + 1L
  N[cbind(events, type[events])] <- 1
  for (j in seq_len(dimens)) N[, j] <- cumsum(N[, j])

  N

}

type_mark_to_Nc <- function(type, mark, dimens, Nc0 = NULL){

  Nc <-  matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(Nc) <- paste0("Nc", 1:dimens)

  if (!is.null(Nc0)) Nc[1,] <- Nc0

  events <- seq_len(max(0L, length(type) - 1L)) + 1L
  Nc[cbind(events, type[events])] <- mark[events]
  for (j in seq_len(dimens)) Nc[, j] <- cumsum(Nc[, j])

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
# Will be removed later
hijack <- function (FUN, ...) {

  if(is.null(FUN)) return(NULL)

  .FUN <- FUN
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  .FUN
}


# Bind parameter values to a model component function by modifying its default arguments.
#
# This utility takes a component function (e.g., impact, kernel, residual) and
# injects parameter values into its formal arguments as new defaults.
# The returned function preserves the original signature, so that during
# likelihood evaluation only runtime data arguments need to be supplied.
#
# This design enables a clean separation between:
#   - model parameters (bound once per optimizer iteration), and
#   - data/state inputs (supplied repeatedly inside the likelihood loop).
#
# Args:
#   FUN      : A function representing a model component.
#   ...      : Named values to be assigned as new default arguments of FUN.
#   .strict  : Logical; if TRUE, all names in ... must already exist in FUN's
#              formal arguments. This prevents silent bugs due to misspelled
#              parameter names and enforces the component interface contract.
#
# Returns:
#   A function with updated default arguments.
#
# Notes:
#   - The original function object is not modified outside this scope.
#   - The function environment is preserved.
#   - This operation is lightweight and intended to be called once per
#     parameter update (e.g., inside the optimizer), not inside event loops.
#
# Typical usage in the hspec reconstruction step:
#   impact0 <- bind_component(impact, param = current_param)
#   impact0(alpha, n, mark)  # param is already bound
#
bind_component <- function(FUN, ..., .strict = TRUE) {

  # Return NULL unchanged (used for optional components)
  if (is.null(FUN)) return(NULL)

  # Ensure a valid function input
  if (!is.function(FUN)) {
    stop("`FUN` must be a function, got: ", class(FUN)[1])
  }

  # Collect new default values
  defaults <- list(...)
  if (!length(defaults)) return(FUN)

  arg_names <- names(defaults)

  # All arguments must be explicitly named
  if (is.null(arg_names) || any(arg_names == "")) {
    stop("All arguments in `...` must be named")
  }

  # Extract current formal arguments
  fmls <- formals(FUN)
  fml_names <- names(fmls)

  # In strict mode, prevent introducing unknown arguments
  if (.strict) {
    unknown <- setdiff(arg_names, fml_names)
    if (length(unknown)) {
      stop(
        "Unknown argument(s): ",
        paste(unknown, collapse = ", "),
        "\nAvailable: ",
        paste(fml_names, collapse = ", ")
      )
    }
  }

  # Update defaults
  fmls[arg_names] <- defaults
  formals(FUN) <- fmls

  # Return the modified function
  FUN
}




eval_param <- function(p) {
  if (is.function(p) && identical(names(formals(p)), "param")) evalf(p)
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

    slot0 <- bind_component(object_slot, param = pr_param)

  } else{

    if(reduced){

      slot0 <- matrix(pr_param[look_up_mtrx(object_slot, param_name)], nrow=dimens)

    } else {

      slot0 <- matrix(pr_param, nrow=dimens)

    }
  }
  slot0
}

# Integrated exponential components, including the small-decay limit.
integrate_exp_kernel <- function(component, beta, dt) {
  scaled_time <- beta * dt
  weight <- -expm1(-scaled_time) / beta
  zero <- scaled_time == 0
  weight[zero] <- rep_len(dt, length(scaled_time))[zero]
  result <- component * weight
  result[component == 0] <- 0
  result
}

log_complement <- function(log_probability) {
  result <- log_probability
  small <- !is.na(log_probability) & log_probability < -log(2)
  result[small] <- log1p(-exp(log_probability[small]))
  result[!small] <- log(-expm1(log_probability[!small]))
  result
}

distribution_parameters <- function(FUN, param) {
  as.list(param[intersect(names(param), names(formals(FUN)))])
}

residual_density <- function(FUN, x, log = FALSE, params = list()) {
  args <- c(list(x), params)
  native_log <- "log" %in% names(formals(FUN))
  if (native_log) args$log <- log
  value <- do.call(FUN, args)
  if (log && !native_log) log(value) else value
}

residual_cdf <- function(FUN, q, lower.tail = TRUE, log.p = FALSE, params = list()) {
  args <- c(list(q), params)
  options <- names(formals(FUN))
  native_tail <- "lower.tail" %in% options
  native_log <- "log.p" %in% options
  if (native_tail) args$lower.tail <- lower.tail
  if (native_log) args$log.p <- TRUE
  value <- do.call(FUN, args)
  if (!native_tail && !lower.tail) {
    value <- if (native_log) log_complement(value) else log1p(-value)
  } else if (!native_log) {
    value <- log(value)
  }
  if (log.p) value else exp(value)
}

residual_quantile <- function(FUN, p, lower.tail = TRUE, log.p = FALSE, params = list()) {
  options <- names(formals(FUN))
  native_tail <- "lower.tail" %in% options
  native_log <- "log.p" %in% options
  value <- if (log.p) p else log(p)
  if (!native_tail && !lower.tail) value <- log_complement(value)
  args <- c(list(if (native_log) value else exp(value)), params)
  if (native_tail) args$lower.tail <- lower.tail
  if (native_log) args$log.p <- TRUE
  do.call(FUN, args)
}
