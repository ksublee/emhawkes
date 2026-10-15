#' @include hspec.R hmoment.R hllf.R utilities.R hinterface.R huncertainty.R
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
#' @param N Optional cumulative event-count matrix, with one row per input row
#' and one column per event type. Each row contains counts after that event.
#' @param Nc Optional cumulative mark-total matrix with the same shape as `N`.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' A finite scalar is expanded to every kernel component. A matrix must match
#' the shape of `alpha` and `beta`. For the built-in likelihood, supply this
#' argument when `type_col_map` is set.
#' @param N0 Initial event counts: a finite scalar or one value per event type.
#' If omitted, use the first row of supplied `N`, or zero when `N` is absent.
#' @param Nc0 Initial mark totals: a finite scalar or one value per event type.
#' If omitted, use the first row of supplied `Nc`, or zero when `Nc` is absent.
#' This argument must be named.
#' @param mylogLik A user-defined log-likelihood function, which must accept an `object` argument consistent with \code{object}.
#' @param reduced Logical; if `TRUE`, equal entries within each parameter matrix
#' share one parameter. If `FALSE`, each matrix entry is estimated separately.
#' Function components share parameters by name. Conflicting starting values for
#' a shared name produce a warning; the first value is used.
#' @param constraints `NULL` or a constraint list passed to [maxLik::maxLik()].
#' @param method The optimization method to be used. Refer to \code{\link[maxLik]{maxLik}} for more details.
#' @param grad Optional gradient function passed to [maxLik::maxLik()]. Its first
#' argument is the optimization parameter vector. Takes precedence over `gradient`.
#' @param gradient Named-only choice, `"numeric"` (default) or `"analytic"`.
#' Request `"analytic"` for the Rcpp score of a Hawkes specification.
#' Parameter-only functions require component derivatives through `jacobian`.
#' Supports shared and fixed parameters, numeric mark effects, multi-kernel mappings,
#' and explicit or differentiable automatic initial states. Requires the `"auto"`
#' or `"cpp"` engine. Unsupported specifications produce an error.
#' @param jacobian Named-only list of derivatives for parameter-only `mu`, `alpha`,
#' `beta`, and `eta` functions when `gradient = "analytic"`. Each entry is a numeric
#' matrix or a function with the single argument `param` returning that matrix.
#' Rows follow `as.vector()` of the component matrix (column-major order); columns
#' are named with the component's local `param` names. Column order is arbitrary.
#' Required for each function with nonempty `param`. For an observed-state `impact`,
#' also supply `jacobian$impact`; see Details for its full-parameter convention.
#' Unused with a supplied `grad` or numerical differentiation.
#' @param hess Optional Hessian function passed to [maxLik::maxLik()]. Its first
#' argument is the optimization parameter vector. If `NULL`, use numerical differentiation.
#' @param final_hessian Named-only choice, `"maxLik"` (default) or `"stable"`.
#' The default retains the final Hessian and covariance calculation from `maxLik`.
#' Request `"stable"` to recalculate numerical final curvature for the built-in
#' Hawkes likelihood using scaled central differences and Richardson extrapolation.
#' User-supplied `hess`, `finalHessian = FALSE`, BHHH curvature, custom likelihoods,
#' and residual-distribution likelihoods retain their calculation choices.
#' @param verbose Logical; if `TRUE`, prints the progress of the estimation process.
#' @param ... Additional parameters for optimization. Refer to \code{\link[maxLik]{maxLik}} for more details.
#'
#' @return A \code{\link[maxLik]{maxLik}} object with an `emhawkes` attribute
#' containing the original specification, parameter map, normalized initial-state
#' arguments, constraints, `reduced` setting, supplied analytic Jacobians, and final
#' curvature diagnostics and recorded calculation methods.
#' Use [fit_diagnostics()] to inspect the methods and checks and [fitted_hspec()] to restore the estimated
#' specification and [fit_initial_state()] to retrieve its initial-state arguments.
#'
#' @details
#' Observations are checked once before optimization. The first row is the initial
#' state, and supplied counting processes contain values after each event.
#' Omit `type` for a univariate model and `mark` for unit marks. The observation
#' vectors have the same length. The likelihood ends at the last supplied event.
#'
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
#' To hold optimization parameters at their starting values, pass `fixed` through
#' `...` as described by [maxLik::maxLik()]. Fixing a shared parameter fixes all
#' entries that use it.
#'
#' The built-in likelihood uses the calculation engine selected by
#' `options(emhawkes.engine = "auto")` by default. It selects C++ for supported
#' numeric and callback models described below, and R otherwise. Setting `"r"`
#' forces the R likelihood for comparison or troubleshooting. Setting `"cpp"`
#' requires C++ support and errors on unsupported models; it is not needed for
#' automatic acceleration. Inspect `fit_diagnostics(fit)$computation$engine`
#' for actual engines and fallback reasons. A custom `mylogLik` bypasses this
#' selection. Standalone [logLik()] supports fewer C++ models than `hfit()`.
#' With `gradient = "analytic"`, Rcpp computes the likelihood and its gradient
#' together and reuses the result at the same parameter vector. Matrix scores are
#' summed according to the saved parameter map. Automatic initial states include
#' their parameter derivatives; a singular automatic initialization requires an
#' explicit `lambda_component0`. State-dependent baselines, impacts that depend on
#' latent intensities, mark densities, residual distributions, and custom likelihoods are not supported
#' by this option. Equality-constrained optimization is also excluded.
#' BHHH requires event scores and cannot use this total score.
#' A supplied `grad` takes precedence. A supplied `hess` is retained; otherwise
#' `maxLik` can calculate curvature by differentiating the analytic gradient.
#' Optimization paths and last digits can differ from numerical differentiation.
#'
#' Analytic function components must have the single argument `param` and return
#' plain numeric matrices with constant dimensions. Their values and Jacobians
#' must be deterministic, differentiable functions of `param`; captured constants
#' must remain fixed. Supply exact derivatives in `jacobian`; the package does not
#' automatically differentiate function bodies or verify the supplied formulas.
#' Fixed entries have zero derivative. Components with `param = numeric(0)` need
#' no Jacobian. Contributions for shared parameter names are added across components.
#' Missing or malformed Jacobians produce an error. For these functions, the
#' analytic option evaluates their matrices in R and the likelihood and component
#' scores in Rcpp. The default `emhawkes.engine = "auto"` also selects C++ for
#' supported function-valued models with numerical differentiation. No Jacobian
#' is required. Unsupported specifications use R; `"cpp"` requires C++ support,
#' and `"r"` selects the R likelihood. User callback errors propagate.
#' This option supports parameter-only `mu`, `alpha`, `beta`, and `eta` functions
#' returning plain numeric matrices with constant dimensions. They must be
#' deterministic functions of `param`, with captured values held fixed.
#' Observed-state `impact` functions are also supported without Jacobians.
#' Observed baselines and sequential state-dependent impacts are also supported
#' as described below. Mark densities and residual distributions are excluded.
#' Supplied `grad` and `hess` are retained. Standalone `logLik()` also uses C++
#' for parameter-only component functions, observed baselines, and observed or
#' sequential state-dependent impacts, including restored fitted models. Automatic
#' fallback through `...` assumes deterministic callbacks with fixed captures;
#' it may evaluate a callback again when an unsupported state read is detected.
#'
#' Observed-state `impact` functions may use `n`, `N`, `Nc`, `mark`, `type`,
#' `inter_arrival`, and the resolved `mu`, `alpha`, and `beta` matrices.
#' They return a finite, plain numeric matrix with the kernel dimensions. As in the R
#' likelihood, `N[n, ]` and `Nc[n, ]` include the current event; use row `n - 1`
#' for counts before it. The impact is added after evaluating that event's intensity,
#' only in columns selected by its type. The last event ends the likelihood window.
#' Automatic initialization uses the existing `mu`/`alpha`/`beta` rule; the
#' additional impact does not enter that rule.
#'
#' With `gradient = "analytic"`, `jacobian$impact` is a matrix or a function with
#' the same observed context.
#' Its `param` argument receives the full named optimization vector, whereas
#' the impact function itself receives its local `param` vector as before.
#' Rows follow `as.vector(impact(...))`; columns must include every fitted parameter
#' name, with zeros for no effect. Return the total derivative of the impact matrix,
#' including any dependence through `mu`, `alpha`, or `beta`. Contributions for
#' shared names belong in the same column. This convention differs from the local
#' Jacobians for parameter-only components. A Jacobian is required even when the
#' impact has no local parameters; a parameter-independent impact has a zero matrix.
#'
#' Declare needed `N` and `Nc` arguments explicitly in the impact function, or
#' supply these paths. Captured inputs must be fixed, and callbacks must be
#' deterministic. Ignored `...` is allowed; reading `lambda`, `lambda_component`,
#' or `lambda_component_n`, including through `...`, is unsupported by the
#' analytic option and observed-only precomputation. A callback
#' that forces all of `...` also forces these unavailable states and produces an
#' error. R evaluates impacts and Jacobians; Rcpp computes intensity recursions,
#' component scores, and reverse sensitivities to event impacts. Temporary storage
#' grows with observations times kernel components, without an additional
#' observations-by-components-by-parameters array. Standalone `logLik()` also
#' supports these observed-state callbacks without Jacobians. Default fitting uses
#' C++ for supported callbacks with
#' numerical differentiation and requires no Jacobian.
#'
#' With numerical differentiation and the C++ engine, `mu` may depend on observed
#' `n`, `N`, `Nc`, `mark`, `type`, and `inter_arrival`. Return a plain numeric
#' dimensions-by-one matrix, constant over each inter-event interval. The call
#' at index `n` supplies the baseline for the interval ending at event `n` and
#' that event's intensity; use counts at `n - 1` for its predictable baseline.
#' Initial values use `mu(n = 1)`. Latent-intensity dependencies and within-interval
#' changes are unsupported for these baselines. Ignored `...` is allowed.
#'
#' Numerical C++ fitting also supports an `impact` that explicitly declares
#' `lambda_component_n`, `lambda`, or `lambda_component`. Rcpp calls it after
#' computing the current event's intensity and before applying its jumps.
#' `lambda_component_n` is the dimensions-by-kernels matrix just before the jump.
#' `lambda` and `lambda_component` are allocated only when explicitly declared;
#' they contain the initial and current pre-jump rows, with future rows zero.
#' Component histories use the R engine's row-major column order and names.
#' Retained callback arguments are snapshots. Prefer `lambda_component_n` when
#' only the current state is needed. A full component history has storage
#' proportional to events times dimensions times kernels. If callbacks retain
#' history arguments, snapshot-preserving copies can make cumulative allocations
#' grow quadratically with event count;
#' this is separate from peak live memory. Current-state callbacks avoid this
#' history cost. Numerical differentiation repeats it at each likelihood evaluation.
#' Undeclared latent access through `...` does not select this sequential path.
#' When `mu` is an observed-baseline function, the impact's `mu` argument is the
#' parameter-bound function, as in the R engine. Callbacks must be deterministic
#' and keep captured inputs fixed. Analytic gradients do not support this feedback.
#'
#' The default uses numerical differentiation and does not run
#' additional curvature checks. With `final_hessian = "stable"`, final curvature is
#' checked after optimization without changing coefficients, log-likelihood,
#' the stored optimizer gradient, return code, or iteration count.
#' The stabilization selects difference steps and checks their consistency.
#' If default numerical Hawkes standard errors contain `NaN`, a short message
#' suggests this option.
#'
#' The requested stabilization compares information, inverse information, and
#' standard errors across difference steps. Failed checks, singular or indefinite
#' information, and equality or nearby inequality constraints yield `NA` covariance
#' for free parameters. Fixed parameters retain zero covariance. Supplied Hessians
#' are retained and checked for invertible positive information. No ridge or
#' eigenvalue clipping is applied. Use [fit_diagnostics()] to inspect the reason.
#' A numerically stable Hessian does not establish identifiability, stationarity,
#' a global maximum, or correct confidence-interval coverage. Weak excitation,
#' nearby kernel decay rates, short samples, boundary estimates, and nonsmooth
#' callbacks can make Wald inference unreliable. Consider profile likelihood and
#' a bootstrap that preserves the observation design and initial state.
#'
#' [fitted_hspec()] uses the saved parameter map to restore the estimated model.
#' [fit_initial_state()] retrieves the observation window's initial-state arguments.
#' These are suitable for re-evaluating the likelihood and starting a new simulation.
#' Keep observations separately and supply them again for inference. The stored
#' optimization objective can retain data in its function environment.
#' Save the fit with `saveRDS()` to retain its model, parameter map, and initial-state information.
#'
#' @docType methods
#' @rdname hfit
#' @importFrom maxLik maxLik
#' @export
#'
#' @seealso [fit_diagnostics()], [fitted_hspec()], [fit_initial_state()], [logLik()],
#' \code{\link{hspec-class}}, \code{\link{hsim,hspec-method}}
#'
#' @examples
#'
#' # example 1
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' fit <- hfit(h, inter_arrival = res$inter_arrival, type = res$type)
#' summary(fit)
#' # Optional Rcpp gradient for this numeric specification:
#' # hfit(h, res$inter_arrival, res$type, gradient = "analytic")
#' fit_diagnostics(fit)$status
#' estimated <- fitted_hspec(fit)
#' initial <- fit_initial_state(fit)
#' inferred <- do.call(logLik, c(list(object = estimated,
#'   inter_arrival = res$inter_arrival, type = res$type, infer = TRUE), initial))
#' c(fitted = as.numeric(logLik(fit)), restored = inferred$loglikelihood)
#'
#'
#' # A parameter-only function can keep selected entries fixed at zero.
#' alpha_fun <- function(param = c(a = 0.2, b = 0.1)) {
#'   matrix(c(param["a"], 0, param["b"], param["a"]), 2)
#' }
#' alpha_jac <- cbind(a = c(1, 0, 0, 1), b = c(0, 0, 1, 0))
#' h_fun <- new("hspec", mu = mu, alpha = alpha_fun, beta = beta)
#' # hfit(h_fun, res$inter_arrival, res$type, gradient = "analytic",
#' #      jacobian = list(alpha = alpha_jac))
#'
#' # Default auto selects C++ with numerical differentiation; no Jacobian:
#' # fit_fun <- hfit(h_fun, res$inter_arrival, res$type)
#'
#' # Observed counts can change the excitation after an event.
#' count_impact <- function(param = c(a = 0.2), n, N, ...) {
#'   matrix(param["a"] * (1 + 0.2 * (N[n, 1] %% 2)))
#' }
#' count_jac <- function(param, n, N, ...) {
#'   out <- matrix(0, 1, length(param), dimnames = list(NULL, names(param)))
#'   out[, "a"] <- 1 + 0.2 * (N[n, 1] %% 2)
#'   out
#' }
#' h_count <- new("hspec", mu = 0.3,
#'   alpha = function(param = numeric(0)) matrix(0), beta = 1.2,
#'   impact = count_impact)
#' # hfit(h_count, inter_arrival, gradient = "analytic",
#' #      jacobian = list(impact = count_jac))
#' # Or keep numerical differentiation without a Jacobian:
#' # hfit(h_count, inter_arrival)
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
                            verbose = FALSE, ..., Nc0 = NULL,
                            final_hessian = c("maxLik", "stable"),
                            gradient = c("numeric", "analytic"), jacobian = NULL) standardGeneric("hfit"))

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
           verbose = FALSE, ..., Nc0 = NULL,
           final_hessian = c("maxLik", "stable"),
           gradient = c("numeric", "analytic"), jacobian = NULL){

    final_hessian <- match.arg(final_hessian)
    gradient <- match.arg(gradient)
    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0")

      lambda_component0 <- additional_argument[["lambda0"]]
      additional_argument[["lambda0"]] <- NULL

    }

    # legacy argument 泥섎━
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



    execution <- new.env(parent = emptyenv())
    execution$used <- character(0)
    execution$fallback_reason <- NULL
    execution$gradient_attribute <- execution$hessian_attribute <- FALSE

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
                            callbacks = callbacks, show_initial_message = FALSE, execution = execution,
                            replay = FALSE)

      } else {

        context <- c(list(object = hspec0), observation,
                     list(lambda_component0 = lambda_component0, param = param, reduced = reduced))
        record_fit_engine(execution, "custom")
        logl <- call_component(mylogLik, custom_plan, context)
        if (!is.null(attr(logl, "gradient", exact = TRUE))) execution$gradient_attribute <- TRUE
        if (!is.null(attr(logl, "hessian", exact = TRUE))) execution$hessian_attribute <- TRUE
      }

      if(verbose){
        cat("Parameters : ", param, "\n")
        cat("Log likelihood : ", logl, "\n")
      }

      logl

    }


    gradient_used <- if (!is.null(grad)) "user" else gradient
    engine <- if (is.null(mylogLik)) {
      match.arg(getOption("emhawkes.engine", "auto"), c("auto", "r", "cpp"))
    } else NULL
    execution$requested <- engine
    if (identical(gradient_used, "analytic")) {
      if (!is.null(mylogLik)) stop("Analytic gradient does not support a custom likelihood.")
      if (!is.null(constraints$eqA) || !is.null(constraints$eqB)) {
        stop('Analytic gradient does not support equality constraints; use gradient = "numeric".')
      }
      if (engine == "r") stop('Analytic gradient requires the C++ engine; use emhawkes.engine = "auto" or "cpp".')
      if (identical(toupper(method), "BHHH") || isTRUE(additional_argument$bhhhHessian) ||
          identical(toupper(as.character(additional_argument$finalHessian)), "BHHH")) {
        stop("Analytic gradient supplies a total score, not the event scores required by BHHH.")
      }
      objective <- make_hawkes_objective(object, observation, lambda_component0, parameter_map, verbose, jacobian, execution = execution)
      llh_function <- objective$loglik
      grad <- objective$gradient
    } else if (!is.null(engine) && engine %in% c("auto", "cpp") && (!is.null(object@impact) || any(vapply(
      list(object@mu, object@alpha, object@beta, object@eta), is.function, logical(1))))) {
      # Reuse native recursions without scores or component Jacobians.
      # Numerical differentiation and supplied grad/hess remain with maxLik.
      if (engine == "cpp") {
        objective <- make_hawkes_objective(object, observation, lambda_component0,
          parameter_map, verbose, derivatives = FALSE, execution = execution)
        llh_function <- objective$loglik
      } else {
        execution$fallback_reason <- auto_hawkes_fit_reason(object, observation)
        objective <- if (is.null(execution$fallback_reason)) tryCatch(make_hawkes_objective(object, observation, lambda_component0,
          parameter_map, verbose, derivatives = FALSE, execution = execution),
          emhawkes_cpp_unsupported = function(e) {
            execution$fallback_reason <- conditionMessage(e)
            NULL
          }) else NULL
        if (!is.null(objective)) {
          r_likelihood <- llh_function
          native_available <- TRUE
          llh_function <- function(param) {
            if (!native_available) return(r_likelihood(param))
            tryCatch(objective$loglik(param),
              emhawkes_cpp_unsupported = function(e) {
                # Some callback capabilities can only be detected on use.
                # Stay on R for all subsequent proposals in this fit.
                native_available <<- FALSE
                execution$fallback_reason <- conditionMessage(e)
                r_likelihood(param)
              })
          }
        }
      }
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
    fit <- do.call(maxLik::maxLik, ml_args)
    if (!is.null(fit)) {
      final <- if (is.null(mylogLik)) {
        finalize_hawkes_fit(fit, object, constraints, hess, additional_argument,
                           method, final_hessian)
      } else {
        list(fit = fit, diagnostics = list(status = "not_checked", reason = "custom_likelihood"))
      }
      fit <- final$fit
      final$diagnostics$computation <- fit_computation(fit, execution, gradient,
        gradient_used, final_hessian, hess, additional_argument, final$diagnostics)
      attr(fit, "emhawkes") <- list(
        version = 1L,
        specification = object,
        parameter_map = parameter_map,
        initial_state = list(lambda_component0 = lambda_component0,
                             N0 = observation$N0, Nc0 = observation$Nc0),
        reduced = reduced,
        gradient = gradient_used,
        jacobian = if (gradient_used == "analytic") jacobian else NULL,
        custom_loglik = !is.null(mylogLik),
        constraints = constraints,
        diagnostics = final$diagnostics
      )
    }
    fit


  }
)



