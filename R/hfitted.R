#' @include hfit.R hinterface.R
NULL

fit_metadata <- function(fit) {
  if (inherits(fit, "maxLik") && !requireNamespace("maxLik", quietly = TRUE))
    stop("Reading this legacy maxLik fit requires the optional maxLik package. New fits can use optimizer = 'stats'.", call. = FALSE)
  metadata <- attr(fit, "emhawkes", exact = TRUE)
  if (!(inherits(fit, "maxLik") || inherits(fit, "emhawkes_stats_fit")) || !is.list(metadata) ||
      !identical(metadata$version, 1L) ||
      !methods::is(metadata$specification, "hspec") ||
      !is.list(metadata$parameter_map) || !is.list(metadata$initial_state)) {
    stop("Fit metadata is unavailable. Use a result from the current hfit().", call. = FALSE)
  }
  metadata
}

# Keep simulation-only functions and update their explicitly shared parameters.
bind_fitted_generator <- function(fun, coefficients) {
  if (is.null(fun) || !("param" %in% names(formals(fun)))) return(fun)
  param <- tryCatch(eval(formals(fun)$param, envir = environment(fun)),
                    error = function(e) NULL)
  if (!is.numeric(param) || !is.null(dim(param)) || is.null(names(param)) ||
      anyNA(names(param)) || any(names(param) == "") || anyDuplicated(names(param))) {
    return(fun)
  }
  shared <- intersect(names(param), names(coefficients))
  if (!length(shared)) return(fun)
  param[shared] <- coefficients[shared]
  bind_component(fun, param = param)
}

#' Restore an Estimated Hawkes Specification
#'
#' Reconstruct an `hspec` object from the coefficients and parameter map saved by
#' [hfit()]. The original specification and fit are not modified.
#'
#' @param fit A result from the current [hfit()].
#' @return An `hspec` object with the estimated parameter values.
#' @details
#' The saved map preserves matrix sharing under both `reduced` settings and named
#' parameters shared between functions. Function bodies and environments are
#' retained; their `param` defaults are updated. Fixed constants remain fixed.
#' The original kernel mapping, `rmark`, and `rresidual` are retained. Named default
#' parameters in these generators are updated when they match coefficient names;
#' other generator parameters keep their original defaults.
#'
#' With `optimizer = "maxLik"`, the fit keeps its `maxLik` classes and methods;
#' reading these fits requires maxLik. The `"stats"` optimizer returns an
#' `emhawkes_stats_fit` with its own coefficient, likelihood, covariance and
#' summary methods, requiring no maxLik API. Its `emhawkes` attribute stores
#' a versioned list with `specification`, `parameter_map`, `initial_state`,
#' `reduced`, and `custom_loglik`. New fits also retain `constraints` and final
#' `diagnostics`, available through [fit_diagnostics()]. New calculation records
#' are in `fit_diagnostics(fit)$computation`. Fits without the original model
#' metadata cannot be
#' reconstructed by this function. Function environments remain references to
#' the original environments; external state is not frozen.
#' Use `saveRDS()` and `readRDS()` to retain the fit's metadata between sessions.
#' Supply observations again to [logLik()] and keep them separately. The stored
#' optimization objective can retain data in its function environment.
#'
#' Restoration does not assert convergence or stationarity. Check the fit's
#' convergence information and whether the model is suitable for simulation.
#' A custom `mylogLik` is not part of the restored `hspec`.
#' @seealso [fit_initial_state()], [hfit()], [logLik()], [hsim()]
#' @export
#' @examplesIf requireNamespace("maxLik", quietly = TRUE)
#' h <- new("hspec", mu = 0.3, alpha = 0, beta = 1)
#' fit <- hfit(h, c(0, 0.2, 0.3, 0.5, 1), lambda_component0 = 0,
#'             fixed = c(FALSE, TRUE, TRUE))
#' estimated <- fitted_hspec(fit)
#' initial <- fit_initial_state(fit)
#' do.call(logLik, c(list(object = estimated,
#'                       inter_arrival = c(0, 0.2, 0.3, 0.5, 1)), initial))
fitted_hspec <- function(fit) {
  metadata <- fit_metadata(fit)
  coefficients <- stats::coef(fit)
  expected <- names(metadata$parameter_map$start)
  if (!is.numeric(coefficients) || !is.null(dim(coefficients)) ||
      length(coefficients) != length(expected) || any(!is.finite(coefficients)) ||
      (!is.null(names(coefficients)) &&
       (anyDuplicated(names(coefficients)) || !setequal(names(coefficients), expected)))) {
    stop("Fit coefficients do not match the saved parameter map or are not finite.", call. = FALSE)
  }
  if (is.null(names(coefficients))) names(coefficients) <- expected
  values <- restore_parameters(metadata$parameter_map, coefficients)
  original <- metadata$specification
  # Use the same constructor as optimization, including function-defined shapes.
  do.call(methods::new, c(list(Class = "hspec"), values,
    list(type_col_map = original@type_col_map,
         rmark = bind_fitted_generator(original@rmark, coefficients),
         rresidual = bind_fitted_generator(original@rresidual, coefficients))))
}

#' Retrieve Initial-State Arguments from a Hawkes Fit
#'
#' Retrieve the normalized initial-state arguments saved by [hfit()].
#'
#' @param fit A result from the current [hfit()].
#' @return A list with `lambda_component0`, `N0`, and `Nc0`, suitable for passing
#' to [logLik()] or [hsim()] with `do.call()`.
#' @details
#' Supplied initial intensity components exclude the baseline `mu` and are stored
#' as a matrix. If the argument was
#' omitted, the selected mean rule is stored as `"mean"`, or the selected zero
#' state as a matrix. Mean components are recomputed with the restored model's
#' parameters; the automatic selection is not repeated. This is not a matrix
#' evaluated at the original starting values.
#' Older saved fits without an initialization-policy record return a matrix
#' evaluated under their original initialization rule at the fitted parameters.
#' Their retained objective supplies any observation context needed by callbacks.
#'
#' Initial counts are normalized per dimension and include values inferred from
#' supplied `N` and `Nc` paths. With a custom likelihood and no event observations,
#' the count arguments are kept as supplied. These arguments describe the initial
#' state of the fitted observation window, not the state after its last event.
#' @seealso [fitted_hspec()], [hfit()], [logLik()], [hsim()]
#' @export
fit_initial_state <- function(fit) {
  metadata <- fit_metadata(fit)
  initial <- metadata$initial_state
  if (!("initialization" %in% names(metadata)) && is.null(initial$lambda_component0) &&
      !isTRUE(metadata$custom_loglik)) {
    object <- fitted_hspec(fit)
    observation <- get0("observation", envir = environment(fit_objective(fit)), inherits = TRUE)
    if (!is.list(observation)) observation <- list()
    values <- setting(object)
    context <- observation[intersect(names(observation), c("inter_arrival", "type", "mark", "N", "Nc"))]
    args <- c(list(object = object), context, values[c("mu", "alpha", "beta")], list(.legacy = TRUE))
    initial$lambda_component0 <- normalize_lambda_initial(do.call(get_lambda0, args),
      object@dimens, ncol(values$beta))
  }
  initial
}
