#' @include hengine.R
NULL

# Total Jacobians use the full named optimization vector. State Jacobians use
# as.vector(current components) for both rows and columns, without name mapping.
validate_full_jacobian <- function(value, rows, parameter_names, label, state = FALSE) {
  columns <- if (state) rows else length(parameter_names)
  if (!is.numeric(value) || !is.matrix(value) || is.object(value) ||
      !identical(dim(value), c(as.integer(rows), as.integer(columns))) || any(!is.finite(value)))
    stop(label, " must be a finite plain numeric matrix with ", rows, " rows and ",
      columns, " columns", if (!state) " for all fitted parameters" else " in component order", ".", call. = FALSE)
  if (!state && columns) {
    names <- colnames(value)
    if (is.null(names) || anyNA(names) || anyDuplicated(names) || !setequal(names, parameter_names))
      stop(label, " column names must match all fitted parameters.", call. = FALSE)
    value <- value[, parameter_names, drop = FALSE]
  }
  value
}

# Expand local matrix derivatives once per proposal, before the native event loop.
full_component_jacobian <- function(name, param, plan, parameter_map, rows) {
  component <- parameter_map$components[[name]]
  result <- matrix(0, rows, length(param), dimnames = list(NULL, names(parameter_map$start)))
  if (component$kind == "null") return(result)
  if (component$kind == "function") {
    if (length(component$index)) result[, component$index] <- plan$jacobian(name, param)
  } else if (length(component$index)) result[cbind(seq_len(rows), component$index)] <- 1
  result
}

# Explicit initial components are constants. Automatic components also depend
# on the observed baseline at n = 1 and on alpha/beta through stationarity.
state_initial_sensitivity <- function(values, jac, initial, dimensions, kernels) {
  parameters <- ncol(jac$mu)
  if (!is.null(initial)) return(list(value = initial, jacobian = matrix(0, dimensions * kernels, parameters)))
  mu <- values$mu; alpha <- values$alpha; beta <- values$beta
  if (dimensions == 1L) {
    gap <- beta[1] - alpha[1]
    if (!is.finite(gap) || gap == 0)
      stop("Analytic gradient cannot differentiate a singular automatic initial state; supply lambda_component0.", call. = FALSE)
    value <- matrix((mu[1] * beta[1] / gap - mu[1]) / 2)
    derivative <- alpha[1] / (2 * gap) * jac$mu +
      mu[1] * beta[1] / (2 * gap^2) * jac$alpha -
      mu[1] * alpha[1] / (2 * gap^2) * jac$beta
  } else {
    branching <- alpha / beta
    inverse <- tryCatch(solve(diag(dimensions) - branching), error = function(e) NULL)
    if (is.null(inverse)) stop("Analytic gradient cannot differentiate a singular automatic initial state; supply lambda_component0.", call. = FALSE)
    stationary <- inverse %*% mu
    value <- matrix(rep(stationary, dimensions), dimensions, byrow = TRUE) * branching
    derivative <- matrix(0, dimensions * kernels, parameters)
    for (p in seq_len(parameters)) {
      db <- (matrix(jac$alpha[, p], dimensions) - branching * matrix(jac$beta[, p], dimensions)) / beta
      ds <- inverse %*% (jac$mu[, p] + db %*% stationary)
      derivative[, p] <- as.vector(sweep(db, 2, as.vector(stationary), "*") +
        sweep(branching, 2, as.vector(ds), "*"))
    }
  }
  list(value = value, jacobian = derivative)
}

make_state_gradient_objective <- function(object, observation, lambda_component0,
    parameter_map, verbose = FALSE, jacobian = NULL, execution = NULL) {
  reason <- auto_hawkes_fit_reason(object, observation)
  if (!is.null(reason)) stop("Analytic gradient does not support ", reason, ".", call. = FALSE)
  latent <- c("lambda", "lambda_component", "lambda_component_n")
  formals_impact <- if (is.null(object@impact)) character(0) else names(formals(object@impact))
  if (any(formals_impact %in% c("lambda", "lambda_component")))
    stop("Analytic gradient does not support full intensity histories; use lambda_component_n.", call. = FALSE)
  stateful <- "lambda_component_n" %in% formals_impact
  if (stateful && (!is.list(jacobian$impact) || is.null(names(jacobian$impact)) ||
      anyDuplicated(names(jacobian$impact)) || !setequal(names(jacobian$impact), c("param", "state"))))
    stop("Current latent intensity feedback for analytic gradients requires jacobian$impact = list(param = ..., state = ...).", call. = FALSE)
  plan <- hawkes_component_plan(object, observation, parameter_map, jacobian, TRUE)
  dimensions <- object@dimens; kernels <- ncol(plan$start$beta); components <- dimensions * kernels
  initial <- normalize_lambda_initial(lambda_component0, dimensions, kernels)
  if (is.null(initial) && length(object@type_col_map))
    stop("Analytic gradient requires lambda_component0 for a multi-kernel mapping.", call. = FALSE)
  mapping <- if (is.null(object@type_col_map)) as.list(seq_len(dimensions)) else object@type_col_map
  mapping <- lapply(mapping, function(x) as.integer(x) - 1L)
  types <- as.integer(observation$type); types[1] <- 0L
  parameter_names <- names(parameter_map$start)
  available <- c("param", "n", "mark", "type", "inter_arrival", "N", "Nc", latent, "mu", "alpha", "beta")
  value_args <- if (is.null(object@impact)) NULL else callback_plan(object@impact, setdiff(available, "param"), "impact")
  jac_spec <- if (stateful) jacobian$impact else list(param = jacobian$impact)
  if (!is.null(object@impact) && is.null(jac_spec$param))
    stop("Analytic gradient requires jacobian$impact.", call. = FALSE)
  jac_args <- list()
  for (name in if (stateful) c("param", "state") else if (!is.null(object@impact)) "param" else character(0)) {
    value <- jac_spec[[name]]
    label <- paste0("jacobian$impact", if (stateful) paste0("$", name))
    if (is.function(value)) {
      if (any(names(formals(value)) %in% c("lambda", "lambda_component")) ||
          (!stateful && "lambda_component_n" %in% names(formals(value))))
        stop(label, " cannot depend on unsupported latent intensity histories.", call. = FALSE)
      jac_args[[name]] <- callback_plan(value, available, label)
    } else {
      jac_spec[[name]] <- validate_full_jacobian(value, components, parameter_names,
        label, state = name == "state")
      jac_args[[name]] <- character(0)
    }
  }
  context_at <- function(param, values) {
    context <- c(list(n = 1L), observation[c("mark", "type", "inter_arrival", "N", "Nc")],
      list(lambda = quote(stop_cpp_unsupported("Analytic gradients cannot read full intensity histories.")),
        lambda_component = quote(stop_cpp_unsupported("Analytic gradients cannot read full intensity histories.")),
        lambda_component_n = quote(stop_cpp_unsupported("Observed impacts cannot read latent intensity states."))),
      values[c("mu", "alpha", "beta")])
    if (plan$observed_mu) context$mu <- restore_parameters(parameter_map, param)$mu
    context
  }
  cached_param <- cached <- NULL
  prepare <- function(param, value_only = TRUE) {
    if (!is.null(names(param))) param <- param[parameter_names]
    if (length(param) != length(parameter_names)) stop("Parameter vector does not match the model.")
    param <- unname(param)
    if (identical(param, cached_param)) return(cached)
    invalid <- list(loglikelihood = -Inf, gradient = stats::setNames(rep(NA_real_, length(param)), parameter_names))
    if (any(!is.finite(param))) return(invalid)
    values <- plan$values(param)
    if (any(!is.finite(unlist(values, use.names = FALSE))) || any(values$beta <= 0) || any(values$mu < 0)) return(invalid)
    # Compute only values for line searches. Derivatives are deferred until asked.
    base_initial <- state_initial_sensitivity(values, list(mu = matrix(0, dimensions, 0),
      alpha = matrix(0, components, 0), beta = matrix(0, components, 0)), initial, dimensions, kernels)$value
    if (any(!is.finite(base_initial))) return(invalid)
    baseline <- if (plan$observed_mu) plan$baseline(param, values) else NULL
    if (plan$observed_mu && is.null(baseline)) return(invalid)
    context <- context_at(param, values)
    fun <- if (is.null(object@impact)) NULL else restore_parameters(parameter_map, param)$impact
    # Keep unavailable histories as error promises in value-only passes too.
    numeric_context <- context
    if (stateful) numeric_context["lambda_component_n"] <- list(NULL)
    numeric_args <- numeric_context[value_args]
    if (value_only) record_fit_engine(execution, "cpp")
    value <- if (value_only) hawkes_loglik_cpp(observation$inter_arrival, types, observation$mark,
      as.numeric(values$mu), values$alpha, values$beta, values$eta, base_initial,
      mapping, FALSE, impact_callback = fun, impact_arguments = if (!is.null(fun)) numeric_args else NULL,
      baseline = baseline, strict_callback = TRUE) else NULL
    if (!is.null(value$unsupported)) stop_cpp_unsupported(value$unsupported)
    result <- list(loglikelihood = value$loglikelihood,
      gradient = if (!is.null(value) && !is.finite(value$loglikelihood)) invalid$gradient else NULL,
      values = values, baseline = baseline, context = context, fun = fun,
      initial = base_initial)
    cached_param <<- param; cached <<- result
    result
  }
  evaluate_gradient <- function(param, value_only = TRUE) {
    result <- prepare(param, value_only)
    if (!is.null(result$gradient)) return(result$gradient)
    p <- cached_param; values <- result$values
    core <- lapply(c("mu", "alpha", "beta", "eta"), function(name) {
      rows <- if (name == "mu") dimensions else components
      if (name == "mu" && plan$observed_mu) matrix(0, dimensions, length(p)) else
        full_component_jacobian(name, p, plan, parameter_map, rows)
    })
    names(core) <- c("mu", "alpha", "beta", "eta")
    baseline_jac <- if (plan$observed_mu) plan$baseline_jacobian(p, values) else NULL
    if (plan$observed_mu) core$mu <- matrix(baseline_jac[, 1L], dimensions, length(p))
    init <- state_initial_sensitivity(values, core, initial, dimensions, kernels)
    if (any(!is.finite(init$jacobian))) stop("Nonfinite automatic initial-state derivative; supply lambda_component0.", call. = FALSE)
    context <- result$context
    full_context <- context; full_context$param <- stats::setNames(p, parameter_names)
    runtime_jac_args <- lapply(jac_args, function(args) full_context[args])
    sensitivity <- c(core, list(names = parameter_names, initial = init$jacobian,
      baseline = baseline_jac, stateful = stateful, impact_param = jac_spec$param,
      impact_state = jac_spec$state, impact_arguments = runtime_jac_args))
    record_fit_engine(execution, "cpp")
    scored <- hawkes_loglik_cpp(observation$inter_arrival, types, observation$mark,
      as.numeric(values$mu), values$alpha, values$beta, values$eta, result$initial,
      mapping, FALSE, gradient = TRUE, impact_callback = result$fun,
      impact_arguments = if (!is.null(result$fun)) context[value_args] else NULL,
      baseline = result$baseline, sensitivity = sensitivity)
    if (!is.null(scored$unsupported)) stop_cpp_unsupported(scored$unsupported)
    if (!is.null(result$loglikelihood) && !identical(scored$loglikelihood, result$loglikelihood))
      stop("Analytic likelihood changed between value and gradient passes; callbacks must be deterministic.", call. = FALSE)
    cached$loglikelihood <<- scored$loglikelihood
    cached$gradient <<- scored$gradient
    scored$gradient
  }
  list(loglik = function(param) {
    value <- prepare(param)$loglikelihood
    if (verbose) { cat("Parameters : ", param, "\n"); cat("Log likelihood : ", value, "\n") }
    c(loglikelihood = value)
  }, gradient = evaluate_gradient,
  # Curvature asks for scores at fresh vectors, without a preceding line search.
  # Obtain value and score together; retain the optimizer's value-first behavior.
  score = function(param) evaluate_gradient(param, value_only = FALSE))
}
