#' @include hspec.R hmoment.R hllf.R utilities.R hinterface.R
NULL

#' Perform Maximum Likelihood Estimation
#'
#' Estimate exponential Hawkes model parameters by maximum likelihood.
#' The supplied \code{\link{hspec-class}} object defines the model and the starting
#' values for optimization with \code{\link[maxLik]{maxLik}}.
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
#' @param Nc0 Initial values for the mark-weighted counting processes. This argument must be named.
#' @param mylogLik A user-defined log-likelihood function, which must accept an `object` argument consistent with \code{object}.
#' @param reduced Logical; if `TRUE`, equal entries within each parameter matrix
#' share one parameter. If `FALSE`, each matrix entry is estimated separately.
#' Function components share parameters by name. Conflicting starting values for
#' a shared name produce a warning; the first value is used.
#' @param constraints constraints matrices. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param method The optimization method to be used. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param grad A gradient matrix for the likelihood function. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param hess A Hessian matrix for the likelihood function. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param verbose Logical; if `TRUE`, prints the progress of the estimation process.
#' @param ... Additional parameters for optimization. Refer to \code{\link[maxLik]{maxLik}} for more details.
#'
#' @return \code{\link[maxLik]{maxLik}} object
#'
#' @details
#' Observations are checked once before optimization. The first row is the initial
#' state, and supplied counting processes contain values after each event.
#' Functions with a `param` argument must provide a default named numeric vector;
#' `numeric(0)` and functions without `param` are allowed for fixed components.
#' A custom `mylogLik` receives its declared arguments from `object`, the observation
#' arguments, `lambda_component0`, `param`, and `reduced`. Optional defaults are
#' preserved; `...` receives the remaining context. Additional inputs should be
#' captured in the function's environment or given defaults.
#'
#' When equal matrix entries are grouped, a message lists the affected components
#' once before optimization. Use `suppressMessages(hfit(...))` to hide it.
#' Setting `reduced = FALSE` separates matrix entries; function parameters still
#' share values by name. A numeric zero is a starting value, not a fixed zero.
#' To fix entries at zero, use a component function with constants at those positions.
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
                            grad = NULL, hess = NULL, constraints = NULL, method = "BFGS",
                            verbose = FALSE, ..., Nc0 = NULL) standardGeneric("hfit"))

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
           grad = NULL, hess = NULL, constraints = NULL, method = "BFGS",
           verbose = FALSE, ..., Nc0 = NULL){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0")

      lambda_component0 <- additional_argument[["lambda0"]]
      additional_argument[["lambda0"]] <- NULL

    }

    # legacy argument 처리
    if ("constraint" %in% names(additional_argument)) {

      if (!is.null(constraints)) {
        stop("Use only one of 'constraints' or 'constraint'.")
      }

      constraints <- additional_argument[["constraint"]]
      additional_argument[["constraint"]] <- NULL

      warning("'constraint' is deprecated; use 'constraints' instead.")
    }


    if(verbose == TRUE && is.null(lambda_component0)){
      message("Initial values for intensity processes were not provided.
              Using internally determined default values for estimation.\n")
    }


    parameter_map <- build_parameter_map(object, reduced)
    starting_point <- parameter_map$start
    callbacks <- if (is.null(mylogLik)) model_callback_plans(object) else NULL
    if (is.null(inter_arrival) && !is.null(mylogLik)) {
      # Custom objectives need not use event observations.
      observation <- list(inter_arrival = inter_arrival, type = type, mark = mark,
                          N = N, Nc = Nc, N0 = N0, Nc0 = Nc0)
    } else {
      observation <- normalize_observation(object, inter_arrival, type, mark, N, Nc, N0, Nc0,
                                           needs = required_count_paths(object, mylogLik))
    }
    kernel_cols <- ncol(eval_param(object@beta))
    lambda_component0 <- normalize_lambda_initial(lambda_component0, object@dimens, kernel_cols)
    if (is.null(mylogLik) && is.null(lambda_component0) && length(object@type_col_map)) {
      stop("In this model, please provide lambda_component0.")
    }
    custom_plan <- callback_plan(mylogLik,
      c("object", names(observation), "lambda_component0", "param", "reduced"), "mylogLik")
    if (!is.null(mylogLik) && (!is.function(mylogLik) || !("object" %in% names(formals(mylogLik))))) {
      stop("mylogLik must be a function with an 'object' argument.")
    }



    # Log-likelihood function passed to maxLik
    llh_function <- function(param){

      # This function may be run by maxLik repeatedly.
      # param may be automatically generated by optimizer

      # Reconstruct hspec0 for likelihood function with given param
      # object is defined in the bounding environment, i.e., hfit

      init_values <- restore_parameters(parameter_map, param)
      names(init_values) <- paste0(names(init_values), "0")
      list2env(init_values, envir = environment())

      # Invalid optimizer proposals must not become constructor errors.
      beta_values <- eval_param(beta0)
      if (!is.function(beta_values) &&
          (any(!is.finite(beta_values)) || any(beta_values <= 0))) return(-Inf)
      for (value in list(eval_param(mu0), eval_param(alpha0), eval_param(eta0))) {
        if (!is.function(value) && any(!is.finite(value))) return(-Inf)
      }
      mu_values <- eval_param(mu0)
      if (!is.function(mu_values) && any(mu_values < 0)) return(-Inf)


      # hspec for logLik
      hspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0, eta = eta0,
                             impact = impact0, dmark = dmark0,
                             dresidual = dresidual0, presidual = presidual0, qresidual = qresidual0,
                             type_col_map = object@type_col_map)


      if (is.null(mylogLik)){

        logl <- loglik_hspec(hspec0, observation, lambda_component0,
                            callbacks = callbacks, show_initial_message = FALSE)

      } else {

        context <- c(list(object = hspec0), observation,
                     list(lambda_component0 = lambda_component0, param = param, reduced = reduced))
        logl <- call_component(mylogLik, custom_plan, context)
      }

      if(verbose){
        cat("Parameters : ", param, "\n")
        cat("Log likelihood : ", logl, "\n")
      }

      logl

    }


    #llh_function(starting_point)
    #maxLik::maxLik(logLik=llh_function,
    #               start=starting_point, grad, hess, method = method, constraints = constraints, ...)

    ml_args <- c(
      list(
        logLik = llh_function,
        start  = starting_point,
        grad   = grad,
        hess   = hess,
        method = method,
        constraints = constraints
      ),
      additional_argument
    )

    message_parameter_sharing(parameter_map)
    do.call(maxLik::maxLik, ml_args)


  }
)



