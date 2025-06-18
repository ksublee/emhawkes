#' @include hspec.R hmoment.R hllf.R utilities.R
NULL

#' Perform Maximum Likelihood Estimation
#'
#' This is a generic function named `hfit` designed for estimating the parameters
#' of the exponential Hawkes model. It is implemented as an S4 method for two main reasons:
#'
#' Model Representation: To represent the structure of the model as an `hspec` object.
#' The multivariate marked Hawkes model has numerous variations, and using an S4 class
#' allows for a flexible and structured approach.
#'
#' Optimization Initialization: To provide a starting point for numerical optimization.
#' The parameter values assigned to the `hspec` slots serve as initial values for the optimization process.
#'
#' This function utilizes the \code{\link[maxLik]{maxLik}} package for optimization.
#'
#' @param object An \code{\link{hspec-class}} object containing the parameter values.
#' @param inter_arrival A vector of inter-arrival times for events across all dimensions, starting with zero.
#' @param type A vector indicating the dimensions, represented by numbers like 1, 2, 3, etc., starting with zero.
#' @param mark A vector of mark (jump) sizes, starting with zero.
#' @param N A matrix representing counting processes.
#' @param Nc A matrix of counting processes weighted by mark sizes.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' Can be a numeric value or a matrix.
#' Must have the same number of rows and columns as \code{alpha} or \code{beta} in \code{object}.
#' @param N0 Initial values for the counting processes matrix \code{N}.
#' @param mylogLik A user-defined log-likelihood function, which must accept an `object` argument consistent with \code{object}.
#' @param reduced Logical; if `TRUE`, performs reduced estimation.
#' @param constraint Constraint matrices. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param method The optimization method to be used. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param grad A gradient matrix for the likelihood function. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param hess A Hessian matrix for the likelihood function. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param verbose Logical; if `TRUE`, prints the progress of the estimation process.
#' @param ... Additional parameters for optimization. Refer to \code{\link[maxLik]{maxLik}} for more details.
#'
#' @return \code{\link{maxLik}} object
#'
#' @docType methods
#' @rdname hfit
#' @export
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hsim,hspec-method}}
#'
#' @examples
#'
#' # example 1
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' summary(hfit(h, inter_arrival=res$inter_arrival, type=res$type))
#'
#'
#' # example 2
#' \donttest{
#' mu <- matrix(c(0.08, 0.08, 0.05, 0.05), nrow = 4)
#' alpha <- function(param = c(alpha11 = 0, alpha12 = 0.4, alpha33 = 0.5, alpha34 = 0.3)){
#'   matrix(c(param["alpha11"], param["alpha12"], 0, 0,
#'            param["alpha12"], param["alpha11"], 0, 0,
#'            0, 0, param["alpha33"], param["alpha34"],
#'            0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = TRUE)
#' }
#' beta <- matrix(c(rep(0.6, 8), rep(1.2, 8)), nrow = 4, byrow = TRUE)
#'
#' impact <- function(param = c(alpha1n=0, alpha1w=0.2, alpha2n=0.001, alpha2w=0.1),
#'                    n=n, N=N, ...){
#'
#'   Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
#'                   0, 0, param['alpha1n'], param['alpha1w'],
#'                   param['alpha2w'], param['alpha2n'], 0, 0,
#'                   param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=TRUE)
#'
#'   ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5
#'
#'   km <- matrix(c(!ind, !ind, !ind, !ind,
#'                  ind, ind, ind, ind,
#'                  ind, ind, ind, ind,
#'                  !ind, !ind, !ind, !ind), nrow = 4, byrow = TRUE)
#'
#'   km * Psi
#' }
#' h <- new("hspec",
#'          mu = mu, alpha = alpha, beta = beta, impact = impact)
#' hr <- hsim(h, size=100)
#' plot(hr$arrival, hr$N[,'N1'] - hr$N[,'N2'], type='s')
#' lines(hr$N[,'N3'] - hr$N[,'N4'], type='s', col='red')
#' fit <- hfit(h, hr$inter_arrival, hr$type)
#' summary(fit)
#' }
#'
#' # example 3
#' \donttest{
#' mu <- c(0.15, 0.15)
#' alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=TRUE)
#' beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2, byrow=TRUE)
#' rmark <- function(param = c(p=0.65), ...){
#'   rgeom(1, p=param[1]) + 1
#' }
#' impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
#'   ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#'   alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#' }
#' h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#'           rmark = rmark,
#'           impact=impact)
#' res <- hsim(h1, size=100, lambda_component0 = matrix(rep(0.1,4), nrow=2))
#'
#' fit <- hfit(h1,
#'             inter_arrival = res$inter_arrival,
#'             type = res$type,
#'             mark = res$mark,
#'             lambda_component0 = matrix(rep(0.1,4), nrow=2))
#' summary(fit)
#' }
#'# For more information, please see vignettes.
setGeneric("hfit", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL,
                            N = NULL, Nc = NULL,
                            lambda_component0 = NULL, N0 = NULL,
                            mylogLik = NULL,
                            reduced = TRUE,
                            grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
                            verbose = FALSE, ...) standardGeneric("hfit"))

#' @rdname hfit
setMethod(
  f="hfit",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL,
           N = NULL, Nc = NULL,
           lambda_component0 = NULL, N0 = NULL,
           mylogLik = NULL,
           reduced = TRUE,
           grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
           verbose = FALSE, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0")

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    if(verbose == TRUE && is.null(lambda_component0)){
      message("Initial values for intensity processes were not provided.
              Using internally determined default values for estimation.\n")
    }


    # Define mapping from slot name to prefix
    param_prefix <- list(
      mu = "mu",
      alpha = "alpha",
      beta = "beta",
      eta = "eta",
      impact = "",
      dmark = "",
      dresidual = "",
      presidual = "",
      qresidual = ""
    )


    # Extract each parameter group from the hspec object
    param_list <- lapply(names(param_prefix), function(slot_name) {
      prefix <- param_prefix[[slot_name]]
      slot_value <- slot(object, slot_name)
      as.param(slot_value, prefix = prefix, reduced = reduced)
    })


    # Extract parameter names for each group
    param_names <- lapply(param_list, names)

    # Create a list to store names for lookup during initialization
    names_list <- setNames(vector("list", length(param_prefix)), names(param_prefix))
    for (slot in names(param_prefix)) {
      idx <- match(slot, names(param_prefix))
      names_list[[slot]] <- param_names[[idx]]
    }


    # Flatten all parameters into a named numeric vector
    input_parameters <- unlist(param_list, recursive = FALSE, use.names = TRUE)

    # Rule: Parameters with the same name are treated as identical
    starting_point <- input_parameters[unique(names(input_parameters))]



    # Log-likelihood function passed to maxLik
    llh_function <- function(param){

      # This function may be run by maxLik repeatedly.
      # param may be automatically generated by optimizer

      # Reconstruct hspec0 for likelihood function with given param
      # object is defined in the bounding environment, i.e., hfit

      init_values <- lapply(names(param_prefix), function(slot_name) {
        initialize_slot(slot(object, slot_name), param[names_list[[slot_name]]],
                        param_name = slot_name, object@dimens, reduced)
      })

      names(init_values) <- paste0(names(param_prefix), "0")
      list2env(init_values, envir = environment())


      # hspec for logLik
      hspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0, eta = eta0,
                             impact = impact0, dmark = dmark0,
                             dresidual = dresidual0, presidual = presidual0, qresidual = qresidual0,
                             type_col_map = object@type_col_map)


      if (is.null(mylogLik)){

        logl <- logLik(hspec0, inter_arrival = inter_arrival, type = type,
                       mark = mark, N = N, Nc = Nc, N0 = N0,
                       lambda_component0 = lambda_component0,
                       showWarning = FALSE)

      } else {

        # arguments names for mylogLik
        args_needed <- names(formals(mylogLik))
        if ( !("object" %in% args_needed)) {
          stop('mylogLik needs \'object\' arguments with hspec')
        }
        args_logLik <- vector("list", length(args_needed))
        for (i in seq_along(args_needed)){
          # find necessary object for mylogLik arguments
          if(args_needed[i] == "object") args_logLik[[i]] <- hspec0
          else args_logLik[[i]] <- get(args_needed[i])
        }
        names(args_logLik) <- args_needed
        logl <- do.call(mylogLik, args = args_logLik)
      }

      if(verbose){
        cat("Parameters : ", param, "\n")
        cat("Log likelihood : ", logl, "\n")
      }

      logl

    }


    #llh_function(starting_point)
    maxLik::maxLik(logLik=llh_function,
                   start=starting_point, grad, hess, method = method, constraint = constraint, ...)

  }
)



