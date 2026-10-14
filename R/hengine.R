#' @useDynLib emhawkes, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

# Inspect capabilities without evaluating a user function.
native_hawkes_reason <- function(object, observation,
    values = list(mu = object@mu, alpha = object@alpha, beta = object@beta, eta = object@eta),
    allow_impact = FALSE) {
  if (any(vapply(values, is.function, logical(1)))) return("function-valued parameters")
  if (!is.null(object@impact) && !allow_impact) return("a user-defined impact function")
  if (!is.null(object@dmark)) return("a mark density function")
  if (!is.null(object@dresidual)) return("a residual distribution")
  # Custom numeric classes can override R arithmetic or subsetting.
  if (is.object(object@dimens) || any(vapply(values, is.object, logical(1)))) {
    return("classed numeric parameters")
  }
  if (is.object(object@type_col_map) ||
      any(vapply(object@type_col_map, is.object, logical(1)))) return("a classed kernel mapping")
  if (any(vapply(observation[c("inter_arrival", "type", "mark")], is.object, logical(1)))) {
    return("classed observation vectors")
  }
  dimensions <- object@dimens
  shape <- dim(values$beta)
  if (length(dimensions) != 1L || !is.finite(dimensions) || dimensions < 1 ||
      dimensions != floor(dimensions) || length(shape) != 2L ||
      shape[1] != dimensions || shape[2] < 1L ||
      !((is.matrix(values$mu) && identical(dim(values$mu), c(as.integer(dimensions), 1L))) ||
        (dimensions == 1 && is.numeric(values$mu) && is.null(dim(values$mu)) && length(values$mu) == 1L)) ||
      !identical(dim(values$alpha), shape) ||
      (!is.null(values$eta) && !identical(dim(values$eta), shape))) {
    return("inconsistent numeric parameter dimensions")
  }
  if (any(vapply(values, function(value) any(!is.finite(value)), logical(1))) ||
      any(values$beta <= 0) || any(values$mu < 0)) return("invalid numeric parameters")
  if (is.null(object@type_col_map)) {
    if (shape[2] != dimensions) return("a missing kernel column mapping")
  } else {
    mapping <- object@type_col_map
    if (length(mapping) != dimensions || !all(vapply(mapping, function(columns) {
      is.numeric(columns) && length(columns) > 0L && all(is.finite(columns)) &&
        all(columns == floor(columns) & columns >= 1 & columns <= shape[2]) &&
        !anyDuplicated(columns)
    }, logical(1)))) return("an invalid kernel column mapping")
  }
  NULL
}

loglik_hspec_cpp <- function(object, observation, lambda_component0, infer,
                             show_initial_message) {
  if (is.null(lambda_component0)) {
    if (length(object@type_col_map)) stop("In this model, please provide lambda_component0.")
    if (show_initial_message) {
      message("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
    }
    lambda_component0 <- get_lambda0(object,
      mark = observation$mark, type = observation$type,
      inter_arrival = observation$inter_arrival, N = observation$N, Nc = observation$Nc,
      mu = object@mu, alpha = object@alpha, beta = object@beta)
  }
  dimensions <- object@dimens
  kernels <- ncol(object@beta)
  initial <- normalize_lambda_initial(lambda_component0, dimensions, kernels)
  mapping <- if (is.null(object@type_col_map)) as.list(seq_len(dimensions)) else object@type_col_map
  mapping <- lapply(mapping, function(columns) as.integer(columns) - 1L)
  types <- observation$type
  if (types[1] != 0) types[1] <- 0L
  result <- hawkes_loglik_cpp(observation$inter_arrival, as.integer(types),
    observation$mark, as.numeric(object@mu), object@alpha, object@beta, object@eta,
    initial, mapping, infer)
  if (!infer) return(c(loglikelihood = result$loglikelihood))
  colnames(result$lambda) <- paste0("lambda", seq_len(dimensions))
  colnames(result$rambda) <- paste0("rambda", seq_len(dimensions))
  component_names <- paste0(rep(seq_len(dimensions), each = kernels), seq_len(kernels))
  colnames(result$lambda_component) <- paste0("lambda", component_names)
  colnames(result$rambda_component) <- paste0("rambda", component_names)
  residuals <- result$residuals
  result$residuals <- NULL
  names(residuals) <- paste0("res_process", seq_len(dimensions))
  c(result, residuals)
}

# Prepare the chain-rule adapter once. Component/Jacobian functions must be
# deterministic functions of their named local parameter vector.
hawkes_component_plan <- function(object, observation, parameter_map, jacobian, derivatives = TRUE) {
  component_names <- c("mu", "alpha", "beta", "eta")
  components <- parameter_map$components[component_names]
  functions <- names(components)[vapply(components, function(x) x$kind == "function", logical(1))]
  for (name in functions) {
    if (!identical(names(formals(components[[name]]$value)), "param")) {
      stop("Analytic gradient requires ", name, " to have only the 'param' argument.", call. = FALSE)
    }
  }
  if (is.null(jacobian)) jacobian <- list()
  allowed <- c(functions, if (!is.null(object@impact)) "impact")
  if (!is.list(jacobian) || (length(jacobian) &&
      (is.null(names(jacobian)) || anyNA(names(jacobian)) || any(names(jacobian) == "") ||
       anyDuplicated(names(jacobian)) || any(!names(jacobian) %in% allowed)))) {
    stop("jacobian must be a uniquely named list for function components mu, alpha, beta, eta, or impact.", call. = FALSE)
  }
  local_param <- function(name, param) {
    x <- components[[name]]
    stats::setNames(param[x$index], names(x$param))
  }
  values_at <- function(param) {
    values <- restore_parameters(parameter_map, param)
    for (name in functions) values[[name]] <- components[[name]]$value(param = local_param(name, param))
    values[component_names]
  }
  start_values <- values_at(parameter_map$start)
  shapes <- lapply(start_values, dim)
  validate_values <- function(values) {
    for (name in functions) {
      value <- values[[name]]
      if (!is.numeric(value) || !is.matrix(value) || is.object(value) ||
          !identical(dim(value), shapes[[name]])) {
        stop(name, " must return a plain numeric matrix of constant dimensions for the analytic gradient.", call. = FALSE)
      }
    }
    values
  }
  start_values <- validate_values(start_values)
  reason <- native_hawkes_reason(object, observation, start_values, allow_impact = TRUE)
  if (!is.null(reason)) stop("Analytic gradient does not support ", reason, ".", call. = FALSE)
  validate_jacobian <- function(value, name) {
    param_names <- names(components[[name]]$param)
    p <- length(param_names)
    if (!is.numeric(value) || !is.matrix(value) || is.object(value) ||
        !identical(dim(value), c(length(start_values[[name]]), p)) || any(!is.finite(value))) {
      stop("jacobian$", name, " must be a finite numeric matrix with one row per component entry and one column per local parameter.", call. = FALSE)
    }
    if (p && (is.null(colnames(value)) || anyNA(colnames(value)) ||
        anyDuplicated(colnames(value)) || !setequal(colnames(value), param_names))) {
      stop("jacobian$", name, " column names must match the component's param names.", call. = FALSE)
    }
    if (p) value <- value[, param_names, drop = FALSE]
    value
  }
  if (derivatives) for (name in functions) {
    if (is.null(jacobian[[name]])) {
      if (length(components[[name]]$param)) {
        stop("Analytic gradient requires jacobian$", name, " for this function component.", call. = FALSE)
      }
      jacobian[[name]] <- matrix(numeric(0), length(start_values[[name]]), 0L)
    }
    if (is.function(jacobian[[name]])) {
      if (!identical(names(formals(jacobian[[name]])), "param")) {
        stop("jacobian$", name, " must have only the 'param' argument.", call. = FALSE)
      }
      validate_jacobian(jacobian[[name]](param = local_param(name, parameter_map$start)), name)
    } else jacobian[[name]] <- validate_jacobian(jacobian[[name]], name)
  }
  list(start = start_values, functions = functions,
    values = function(param) validate_values(values_at(param)),
    jacobian = function(name, param) {
      value <- jacobian[[name]]
      if (is.function(value)) value <- validate_jacobian(value(param = local_param(name, param)), name)
      value
    })
}

# Hawkes objective: one native pass supplies both value and component scores.
# The one-point cache is local to a fit and is invalidated by any parameter change.
make_hawkes_objective <- function(object, observation, lambda_component0,
                                 parameter_map, verbose = FALSE, jacobian = NULL, derivatives = TRUE) {
  plan <- hawkes_component_plan(object, observation, parameter_map, jacobian, derivatives)
  impact_plan <- hawkes_impact_plan(object, observation, parameter_map, jacobian$impact,
                                    dim(plan$start$beta), derivatives)
  dimensions <- object@dimens
  kernels <- ncol(plan$start$beta)
  initial <- normalize_lambda_initial(lambda_component0, dimensions, kernels)
  if (is.null(initial) && length(object@type_col_map)) {
    stop("Analytic gradient requires lambda_component0 for a multi-kernel mapping.")
  }
  mapping <- if (is.null(object@type_col_map)) as.list(seq_len(dimensions)) else object@type_col_map
  mapping <- lapply(mapping, function(x) as.integer(x) - 1L)
  types <- observation$type
  types[1] <- 0L
  types <- as.integer(types)
  parameter_names <- names(parameter_map$start)
  cached_param <- cached_result <- NULL

  evaluate <- function(param) {
    if (!is.null(names(param))) param <- param[parameter_names]
    if (length(param) != length(parameter_names)) stop("Parameter vector does not match the model.")
    param <- unname(param)
    if (identical(param, cached_param)) return(cached_result)
    invalid <- list(loglikelihood = -Inf, gradient = stats::setNames(rep(NA_real_, length(param)), parameter_names))
    if (any(!is.finite(param))) return(invalid)
    values <- plan$values(param)
    if (any(!is.finite(unlist(values, use.names = FALSE))) ||
        any(values$beta <= 0) || any(values$mu < 0)) return(invalid)
    current_initial <- initial
    if (is.null(current_initial)) {
      if (dimensions == 1L) {
        gap <- values$beta[1] - values$alpha[1]
        current_initial <- matrix((values$mu[1] * values$beta[1] / gap - values$mu[1]) / 2)
        if (!all(is.finite(current_initial))) {
          stop("Analytic gradient cannot differentiate a singular automatic initial state; supply lambda_component0.", call. = FALSE)
        }
      } else {
        branching <- values$alpha / values$beta
        inverse <- tryCatch(solve(diag(dimensions) - branching), error = function(e) NULL)
        if (is.null(inverse)) {
          stop("Analytic gradient cannot differentiate a singular automatic initial state; supply lambda_component0.", call. = FALSE)
        }
        stationary <- inverse %*% values$mu
        current_initial <- matrix(rep(stationary, dimensions), nrow = dimensions, byrow = TRUE) * values$alpha / values$beta
      }
    }
    # Finite component values can still overflow the automatic initializer.
    if (any(!is.finite(current_initial))) return(invalid)
    impacts <- if (!is.null(impact_plan)) impact_plan$values(param, values) else NULL
    if (!is.null(impact_plan) && is.null(impacts)) return(invalid)
    result <- hawkes_loglik_cpp(observation$inter_arrival, types, observation$mark,
      as.numeric(values$mu), values$alpha, values$beta, values$eta,
      current_initial, mapping, infer = FALSE, gradient = derivatives, impact = impacts)
    if (!is.finite(result$loglikelihood)) {
      return(list(loglikelihood = result$loglikelihood, gradient = invalid$gradient))
    }
    if (!derivatives) {
      cached_param <<- param
      cached_result <<- list(loglikelihood = result$loglikelihood)
      return(cached_result)
    }
    score <- result$score
    if (is.null(initial) && is.finite(result$loglikelihood)) {
      if (dimensions == 1L) {
        q <- score$initial[1]
        score$mu <- score$mu + q * values$alpha[1] / (2 * gap)
        score$alpha <- score$alpha + q * values$mu[1] * values$beta[1] / (2 * gap^2)
        score$beta <- score$beta - q * values$mu[1] * values$alpha[1] / (2 * gap^2)
      } else {
        # C_ij = (alpha_ij / beta_ij) * stationary_j. This adjoint applies
        # the initial-state chain rule without a components-by-parameters array.
        w <- as.vector(crossprod(inverse, colSums(score$initial * branching)))
        g_branching <- sweep(score$initial + w, 2, as.vector(stationary), "*")
        score$mu <- score$mu + w
        score$alpha <- score$alpha + g_branching / values$beta
        score$beta <- score$beta - g_branching * branching / values$beta
      }
    }
    gradient <- stats::setNames(numeric(length(param)), parameter_names)
    for (name in c("mu", "alpha", "beta", "eta")) {
      indices <- parameter_map$components[[name]]$index
      if (length(indices) && name %in% plan$functions) {
        jac <- plan$jacobian(name, param)
        partial <- as.vector(score[[name]])
        contribution <- if (all(is.finite(partial))) as.vector(crossprod(jac, partial)) else {
          weighted <- jac * partial
          weighted[jac == 0] <- 0
          colSums(weighted)
        }
        gradient[indices] <- gradient[indices] + contribution
      } else if (length(indices)) {
        sums <- rowsum(as.vector(score[[name]]), indices, reorder = FALSE)
        positions <- as.integer(rownames(sums))
        gradient[positions] <- gradient[positions] + sums[, 1]
      }
    }
    cached_param <<- param
    cached_result <<- list(loglikelihood = result$loglikelihood, gradient = gradient)
    if (!is.null(impact_plan)) {
      # Line searches often ask only for a value. Retain native jump scores but
      # defer the expensive R Jacobian callbacks until a gradient is requested.
      cached_result$impact <<- list(values = values, score = result$impact_score)
    }
    cached_result
  }
  list(loglik = function(param) {
    value <- evaluate(param)$loglikelihood
    if (verbose) {
      cat("Parameters : ", param, "\n")
      cat("Log likelihood : ", value, "\n")
    }
    c(loglikelihood = value)
  }, gradient = if (derivatives) function(param) {
    result <- evaluate(param)
    if (!is.null(result$impact)) {
      result$gradient <- result$gradient + impact_plan$gradient(
        cached_param, result$impact$values, result$impact$score)
      result$impact <- NULL
      cached_result <<- result
    }
    result$gradient
  } else NULL)
}
