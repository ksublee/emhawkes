#' @useDynLib emhawkes, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

# Inspect capabilities without evaluating a user function.
native_hawkes_reason <- function(object, observation) {
  values <- list(mu = object@mu, alpha = object@alpha, beta = object@beta, eta = object@eta)
  if (any(vapply(values, is.function, logical(1)))) return("function-valued parameters")
  if (!is.null(object@impact)) return("a user-defined impact function")
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
  shape <- dim(object@beta)
  if (length(dimensions) != 1L || !is.finite(dimensions) || dimensions < 1 ||
      dimensions != floor(dimensions) || length(shape) != 2L ||
      shape[1] != dimensions || shape[2] < 1L ||
      !((is.matrix(object@mu) && identical(dim(object@mu), c(as.integer(dimensions), 1L))) ||
        (dimensions == 1 && is.numeric(object@mu) && is.null(dim(object@mu)) && length(object@mu) == 1L)) ||
      !identical(dim(object@alpha), shape) ||
      (!is.null(object@eta) && !identical(dim(object@eta), shape))) {
    return("inconsistent numeric parameter dimensions")
  }
  if (any(vapply(values, function(value) any(!is.finite(value)), logical(1))) ||
      any(object@beta <= 0) || any(object@mu < 0)) return("invalid numeric parameters")
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

# Numeric Hawkes objective: one native pass supplies both value and score.
# The one-point cache is local to a fit and is invalidated by any parameter change.
make_hawkes_objective <- function(object, observation, lambda_component0,
                                 parameter_map, verbose = FALSE) {
  reason <- native_hawkes_reason(object, observation)
  if (!is.null(reason)) stop("Analytic gradient does not support ", reason, ".", call. = FALSE)
  dimensions <- object@dimens
  kernels <- ncol(object@beta)
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
    values <- restore_parameters(parameter_map, param)
    invalid <- list(loglikelihood = -Inf, gradient = stats::setNames(rep(NA_real_, length(param)), parameter_names))
    if (any(!is.finite(param)) || any(values$beta <= 0) || any(values$mu < 0)) return(invalid)
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
    result <- hawkes_loglik_cpp(observation$inter_arrival, types, observation$mark,
      as.numeric(values$mu), values$alpha, values$beta, values$eta,
      current_initial, mapping, infer = FALSE, gradient = TRUE)
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
      if (length(indices)) {
        sums <- rowsum(as.vector(score[[name]]), indices, reorder = FALSE)
        positions <- as.integer(rownames(sums))
        gradient[positions] <- gradient[positions] + sums[, 1]
      }
    }
    cached_param <<- param
    cached_result <<- list(loglikelihood = result$loglikelihood, gradient = gradient)
    cached_result
  }
  list(loglik = function(param) {
    value <- evaluate(param)$loglikelihood
    if (verbose) {
      cat("Parameters : ", param, "\n")
      cat("Log likelihood : ", value, "\n")
    }
    c(loglikelihood = value)
  }, gradient = function(param) evaluate(param)$gradient)
}
