#' @include hinterface.R
NULL

# Observed-state impacts are independent of the latent intensity recursion.
# Their Jacobians differentiate the returned impact matrix with respect to the
# full fitted parameter vector, including dependencies through mu/alpha/beta.
hawkes_impact_plan <- function(object, observation, parameter_map, jacobian,
                               shape, derivatives = TRUE) {
  if (is.null(object@impact)) return(NULL)
  component <- parameter_map$components$impact
  parameter_names <- names(parameter_map$start)
  components <- prod(shape)
  size <- length(observation$inter_arrival)
  observed <- c("n", "mark", "type", "inter_arrival", "N", "Nc")
  latent <- c("lambda", "lambda_component", "lambda_component_n")
  # Preserve the existing positional order inside ... as well as argument names.
  available <- c("param", observed, latent, "mu", "alpha", "beta")
  check_callback <- function(fun, label) {
    if (any(names(formals(fun)) %in% latent)) {
      stop(label, " cannot depend on latent intensity states for the analytic gradient.", call. = FALSE)
    }
    callback_plan(fun, available, label)
  }
  value_args <- check_callback(component$value, "impact")
  if (!component$has_param) value_args <- setdiff(value_args, "param")
  validate_jacobian <- function(value) {
    if (!is.numeric(value) || !is.matrix(value) || is.object(value) ||
        !identical(dim(value), c(as.integer(components), length(parameter_names))) ||
        any(!is.finite(value))) {
      stop("jacobian$impact must be a finite numeric matrix with one row per component entry and one column per fitted parameter.", call. = FALSE)
    }
    if (length(parameter_names) && (is.null(colnames(value)) || anyNA(colnames(value)) ||
        anyDuplicated(colnames(value)) || !setequal(colnames(value), parameter_names))) {
      stop("jacobian$impact column names must match the full fitted parameter vector.", call. = FALSE)
    }
    if (length(parameter_names)) value <- value[, parameter_names, drop = FALSE]
    value
  }
  jac_args <- NULL
  if (derivatives) {
    if (is.null(jacobian)) stop("Analytic gradient requires jacobian$impact for an observed-state impact.", call. = FALSE)
    if (is.function(jacobian)) jac_args <- check_callback(jacobian, "jacobian$impact") else
      jacobian <- validate_jacobian(jacobian)
  }
  context_at <- function(param, values, full = FALSE) {
    context <- c(list(n = 1L), observation[c("mark", "type", "inter_arrival", "N", "Nc")],
                 values[c("mu", "alpha", "beta")])
    # do.call supplies these as promises. Ignored ... remains compatible; reading
    # an unavailable intensity through ... fails instead of silently using NULL.
    for (name in latent) context[[name]] <- quote(stop(
      "Analytic impact cannot read latent intensity states through '...'.", call. = FALSE))
    if (full) context$param <- stats::setNames(param, parameter_names) else if (component$has_param) {
      context$param <- stats::setNames(param[component$index], names(component$param))
    }
    context
  }
  list(values = function(param, values) {
    context <- context_at(param, values)
    result <- matrix(0, components, size)
    for (n in seq_len(size - 1L) + 1L) {
      context$n <- n
      value <- call_component(component$value, value_args, context)
      if (!is.numeric(value) || !is.matrix(value) || is.object(value) || !identical(dim(value), shape)) {
        stop("impact must return a plain numeric matrix matching the kernel dimensions.", call. = FALSE)
      }
      if (any(!is.finite(value))) return(NULL)
      result[, n] <- as.vector(value)
    }
    result
  }, gradient = function(param, values, score) {
    context <- context_at(param, values, full = TRUE)
    result <- stats::setNames(numeric(length(parameter_names)), parameter_names)
    for (n in seq_len(size - 1L) + 1L) {
      context$n <- n
      jac <- if (is.function(jacobian)) validate_jacobian(call_component(jacobian, jac_args, context)) else jacobian
      partial <- score[, n]
      if (all(is.finite(partial))) result <- result + as.vector(crossprod(jac, partial)) else {
        weighted <- jac * partial
        weighted[jac == 0] <- 0
        result <- result + colSums(weighted)
      }
    }
    result
  })
}
