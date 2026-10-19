#' @include utilities.R
NULL

# Resolve the existing sharing rules once, before optimization.
build_parameter_map <- function(object, reduced = TRUE) {
  if (!is.logical(reduced) || length(reduced) != 1L || is.na(reduced)) {
    stop("reduced must be TRUE or FALSE.")
  }
  prefixes <- c(mu = "mu", alpha = "alpha", beta = "beta", eta = "eta",
                impact = "", dmark = "", dresidual = "", presidual = "", qresidual = "")
  values <- lapply(names(prefixes), function(name) methods::slot(object, name))
  names(values) <- names(prefixes)
  build_component_parameter_map(values, prefixes, reduced)
}

# Model-independent registration; the hspec adapter keeps the saved map format.
build_component_parameter_map <- function(values, prefixes, reduced, named = character()) {
  components <- lapply(names(prefixes), function(name) {
    value <- values[[name]]
    if (is.null(value)) return(list(kind = "null"))
    if (name %in% named) return(list(kind = "parameters", param = value))
    if (is.function(value)) {
      has_param <- "param" %in% names(formals(value))
      param <- if (has_param) tryCatch(
        eval(formals(value)$param, envir = environment(value)),
        error = function(e) stop(name, " must supply a default numeric 'param' vector.", call. = FALSE)
      ) else numeric(0)
      if (!is.numeric(param) || !is.null(dim(param)) || any(!is.finite(param)) ||
          (length(param) && (is.null(names(param)) || anyNA(names(param)) ||
                             any(names(param) == "") || anyDuplicated(names(param))))) {
        stop(name, "$param must be a finite numeric vector with unique, nonempty names.")
      }
      return(list(kind = "function", value = value, param = param, has_param = has_param))
    }
    param <- as.param(value, prefix = prefixes[[name]], reduced = reduced)
    keys <- if (reduced) look_up_mtrx(value, prefixes[[name]]) else full_names(value, prefixes[[name]])
    list(kind = "matrix", param = param, keys = as.vector(keys), dim = dim(value))
  })
  names(components) <- names(prefixes)
  parameters <- unlist(lapply(components, function(x) x$param), use.names = FALSE)
  parameter_names <- unlist(lapply(components, function(x) names(x$param)), use.names = FALSE)
  names(parameters) <- parameter_names
  first <- !duplicated(parameter_names)
  start <- parameters[first]
  conflicting <- unique(parameter_names[parameters != start[match(parameter_names, names(start))]])
  if (length(conflicting)) {
    warning("Shared parameters have different starting values; using the first value for: ",
            paste(conflicting, collapse = ", "), call. = FALSE)
  }
  components <- lapply(components, function(x) {
    x$index <- match(if (x$kind == "matrix") x$keys else names(x$param), names(start))
    x
  })
  list(start = start, components = components)
}

restore_parameters <- function(map, param) {
  # maxLik uses this order; accept named vectors in any order for internal reuse.
  if (!is.null(names(param))) param <- param[names(map$start)]
  if (length(param) != length(map$start)) stop("Parameter vector does not match the model.")
  lapply(map$components, function(x) {
    if (x$kind == "null") return(NULL)
    if (x$kind == "matrix") return(matrix(param[x$index], nrow = x$dim[1], ncol = x$dim[2]))
    if (x$kind == "parameters") {
      value <- param[x$index]
      names(value) <- names(x$param)
      return(value)
    }
    if (!x$has_param) return(x$value)
    bound <- param[x$index]
    names(bound) <- names(x$param)
    bind_component(x$value, param = bound)
  })
}

#' Register and Restore Model Parameters
#'
#' @param components A named list of numeric matrices/scalars, functions or NULL.
#' @param prefixes Matrix name prefixes, a named character vector matching components.
#' @param named Names of components supplied as explicit named parameter vectors.
#'   Their names are retained and equal values do not imply sharing.
#' @param reduced Share equal starting values within each matrix?
#' @return A locked environment with `start`, the full starting vector, and
#'   `restore(param)`, returning a named list of matrices, functions and vectors.
#' @details Components are visited in list order. Matrix names and sharing match
#'   hfit. Function parameters come from their named numeric `param` defaults,
#'   evaluated in the original environment; function bodies are not called.
#'   Names shared across components identify one parameter. Conflicting starts
#'   warn and use the first value. Functions retain their original environments
#'   and binding semantics. No model object, optimizer, initial-state policy or
#'   random generator is supplied by this API.
#'
#'   Restore accepts the full vector in starting order, or exactly its names in
#'   any order. Nonfinite proposals are retained for the calling objective to
#'   handle. Fixed/free selectors and constraints use this full coordinate order
#'   and are handled by the optimizer layer. Create a new plan when definitions
#'   change. Captured external function environments remain references.
#' @examples
#' plan <- parameter_plan(list(mu = matrix(c(.2, .2), 2),
#'   residual = c(shape = 2, rate = 2)), named = "residual")
#' plan$start
#' plan$restore(plan$start)
#' @export
parameter_plan <- function(components, prefixes = names(components),
                            named = character(), reduced = TRUE) {
  valid_names <- function(x) is.character(x) && !anyNA(x) &&
    all(nzchar(x)) && !anyDuplicated(x)
  if (!is.list(components) || is.object(components) ||
      !length(components) || !valid_names(names(components)))
    stop("components must be a nonempty list with unique, nonempty names.", call. = FALSE)
  if (!is.logical(reduced) || length(reduced) != 1L || is.na(reduced))
    stop("reduced must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(prefixes) || length(prefixes) != length(components) ||
      anyNA(prefixes))
    stop("prefixes must contain one character prefix per component.", call. = FALSE)
  if (is.null(names(prefixes))) names(prefixes) <- names(components)
  if (!valid_names(names(prefixes)) || !setequal(names(prefixes), names(components)))
    stop("prefixes must match the component names.", call. = FALSE)
  prefixes <- prefixes[names(components)]
  if (!valid_names(named) || !all(named %in% names(components)))
    stop("named must identify distinct parameter-vector components.", call. = FALSE)
  for (name in names(components)) {
    value <- components[[name]]
    if (is.null(value)) next
    if (name %in% named) {
      if (!is.numeric(value) || is.complex(value) || is.object(value) ||
          !is.null(dim(value)) || any(!is.finite(value)) ||
          (length(value) && !valid_names(names(value))))
        stop(name, " must be a finite numeric vector with unique parameter names.", call. = FALSE)
    } else if (!is.function(value)) {
      if (!is.numeric(value) || is.complex(value) || is.object(value) ||
          !length(value) || any(!is.finite(value)) ||
          (!is.null(dim(value)) && length(dim(value)) != 2L) ||
          (is.null(dim(value)) && length(value) != 1L))
        stop(name, " must be a finite numeric matrix, scalar, function or NULL.", call. = FALSE)
      if (is.null(dim(value))) components[[name]] <- matrix(value, 1L, 1L)
    }
  }
  map <- build_component_parameter_map(components, prefixes, reduced, named)
  if (is.complex(map$start)) stop("Initial parameters must be real numeric values.", call. = FALSE)
  plan <- new.env(parent = emptyenv())
  plan$start <- if (is.null(map$start)) numeric() else map$start
  plan$restore <- function(param) {
    expected <- names(map$start)
    if (!is.numeric(param) || is.complex(param) || is.object(param) ||
        !is.null(dim(param)) || length(param) != length(map$start) ||
        (!is.null(names(param)) &&
          (!valid_names(names(param)) || !setequal(names(param), expected))))
      stop("param must contain the full parameter vector in starting order or with exactly its names.", call. = FALSE)
    restore_parameters(map, param)
  }
  class(plan) <- "parameter_plan"
  lockEnvironment(plan, bindings = TRUE)
  plan
}

message_parameter_sharing <- function(map) {
  shared <- vapply(map$components, function(x) {
    x$kind == "matrix" && anyDuplicated(x$index) > 0L
  }, logical(1))
  if (any(shared)) {
    message("hfit: equal starting values within each matrix share a parameter ",
            "(reduced = TRUE; ", paste(names(map$components)[shared], collapse = ", "), ").\n",
            "Use reduced = FALSE to estimate matrix entries separately.")
  }
  invisible(NULL)
}

# Compile argument selection without evaluating user defaults.
callback_plan <- function(FUN, available, label) {
  if (is.null(FUN) || !is.function(FUN)) return(NULL)
  fmls <- formals(FUN)
  required <- names(fmls)[vapply(fmls, function(x) identical(x, quote(expr = )), logical(1))]
  missing <- setdiff(required, c(available, "..."))
  if (length(missing)) {
    stop(label, " requires unsupported argument(s): ", paste(missing, collapse = ", "),
         ". Give extra arguments defaults or capture their values in the function's environment.", call. = FALSE)
  }
  if ("..." %in% names(fmls)) available else intersect(names(fmls), available)
}

model_callback_plans <- function(object, simulation = FALSE) {
  available <- c("n", "mark", "type", "inter_arrival", "N", "Nc", "lambda",
                 "lambda_component", "lambda_component_n", "mu", "alpha", "beta")
  slots <- c("mu", "impact", if (simulation) "rmark" else "dmark")
  plans <- lapply(slots, function(name) callback_plan(methods::slot(object, name), available, name))
  names(plans) <- slots
  plans
}

call_component <- function(FUN, plan, context) {
  do.call(FUN, context[plan])
}

# Prepare named argument selection once for repeated calls in the R likelihood.
# forceAndCall snapshots supplied values before invoking the callback; omitted
# optional arguments retain their own defaults and remain unevaluated.
prepare_component_call <- function(plan) {
  invoke <- function(FUN, n, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta) NULL
  arguments <- stats::setNames(lapply(plan, as.name), plan)
  body(invoke) <- as.call(c(list(quote(base::forceAndCall), length(plan), quote(FUN)), arguments))
  invoke
}

required_count_paths <- function(object, extra = NULL) {
  callbacks <- c(list(object@mu, object@impact, object@dmark), list(extra))
  unique(unlist(lapply(callbacks, function(fun) {
    if (!is.function(fun)) return(character(0))
    args <- names(formals(fun))
    if ("..." %in% args) c("N", "Nc") else intersect(args, c("N", "Nc"))
  }), use.names = FALSE))
}

normalize_count_initial <- function(value, dimens, name) {
  if (is.null(value)) return(rep(0, dimens))
  if (!is.numeric(value) || !(length(value) %in% c(1L, dimens)) || any(!is.finite(value))) {
    stop(name, " must be a finite scalar or contain one value per dimension.")
  }
  rep_len(as.numeric(value), dimens)
}

normalize_observation <- function(object, inter_arrival, type = NULL, mark = NULL,
                                  N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
                                  needs = required_count_paths(object)) {
  size <- length(inter_arrival)
  dimens <- object@dimens
  if (!is.numeric(inter_arrival) || !is.null(dim(inter_arrival)) || !size ||
      any(!is.finite(inter_arrival)) || any(inter_arrival < 0) || inter_arrival[1] != 0) {
    stop("inter_arrival must be a finite, nonnegative vector starting with an initial zero.")
  }
  if (is.null(type)) {
    if (dimens != 1L) stop("type is required for a multivariate model; supply one event type per row of inter_arrival, including the initial row.", call. = FALSE)
    type <- c(0, rep(1, size - 1L))
  }
  if (is.null(mark)) mark <- c(0, rep(1, size - 1L))
  if (!is.numeric(type) || !is.null(dim(type)) || length(type) != size || any(!is.finite(type)) ||
      any(type[-1] < 1 | type[-1] > dimens | type[-1] != floor(type[-1]))) {
    stop("type must match inter_arrival (", size, " rows) and contain integer event types from 1 to ",
         dimens, " after the initial row.", call. = FALSE)
  }
  if (!is.numeric(mark) || !is.null(dim(mark)) || length(mark) != size || any(!is.finite(mark))) {
    stop("mark must be finite and match inter_arrival.")
  }
  normalize_path <- function(path, initial, name) {
    if (!is.null(path)) {
      if (!is.matrix(path) || !is.numeric(path) ||
          !identical(dim(path), c(size, as.integer(dimens))) || any(!is.finite(path))) {
        stop(name, " must be a finite matrix with one row per observation and one column per dimension.")
      }
      if (is.null(initial)) initial <- path[1, ]
    }
    initial <- normalize_count_initial(initial, dimens, paste0(name, "0"))
    if (!is.null(path) || name %in% needs) {
      expected <- if (name == "N") type_to_N(type, dimens, initial) else
        type_mark_to_Nc(type, mark, dimens, initial)
      if (is.null(path)) path <- expected else {
        if (!isTRUE(all.equal(as.numeric(path), as.numeric(expected), tolerance = sqrt(.Machine$double.eps)))) {
          stop(name, " is inconsistent with type, mark, or its initial values.")
        }
        if (is.null(colnames(path))) colnames(path) <- colnames(expected)
      }
    }
    list(path = path, initial = initial)
  }
  counts <- normalize_path(N, N0, "N")
  marked_counts <- normalize_path(Nc, Nc0, "Nc")
  list(inter_arrival = inter_arrival, type = type, mark = mark,
       N = counts$path, Nc = marked_counts$path, N0 = counts$initial, Nc0 = marked_counts$initial)
}

normalize_lambda_initial <- function(value, dimens, kernel_cols, legacy = FALSE) {
  if (is.null(value)) return(NULL)
  size <- dimens * kernel_cols
  if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
    stop("lambda_component0 must contain finite numeric values.")
  }
  if (length(value) == 1L) return(matrix(value, dimens, kernel_cols))
  shape_ok <- is.null(dim(value)) || identical(dim(value), c(as.integer(dimens), as.integer(kernel_cols)))
  if (length(value) != size || !shape_ok) {
    expected <- paste0(dimens, " x ", kernel_cols)
    received <- if (is.null(dim(value))) paste0("length ", length(value)) else paste(dim(value), collapse = " x ")
    if (!legacy) stop("lambda_component0 must be a scalar or match the kernel matrix dimensions (",
                      expected, "; received ", received, ").", call. = FALSE)
    warning("lambda_component0 does not match the kernel dimensions (", expected, "; received ",
            received, "); recycling or truncating it is deprecated. Supply a scalar or a matrix with the correct dimensions.", call. = FALSE)
    value <- rep_len(as.numeric(value), size)
  }
  matrix(value, dimens, kernel_cols)
}
