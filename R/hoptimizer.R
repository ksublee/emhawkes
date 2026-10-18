#' @include hinterface.R
NULL

# Access legacy fits only through maxLik; the stats path never needs its API.
fit_active <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$active else maxLik::activePar(fit)
}
fit_objective <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$objectiveFn else maxLik::objectiveFn(fit)
}
fit_hessian <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$hessian else maxLik::hessian(fit)
}
fit_gradient <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$gradient else maxLik::gradient(fit)
}
fit_code <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$code else maxLik::returnCode(fit)
}
fit_message <- function(fit) {
  if (inherits(fit, "emhawkes_stats_fit")) fit$message else maxLik::returnMessage(fit)
}

# Saved legacy fits retain maxLik's classes. Install lazy bridges only when
# the provider has not registered its methods; loading it replaces the bridges.
# Merely loading emhawkes or using its stats optimizer never loads maxLik.
maxlik_saved_method <- function(generic) {
  if (!requireNamespace("maxLik", quietly = TRUE))
    stop("Reading this maxLik fit requires the optional maxLik package. New fits can use optimizer = 'stats'.", call. = FALSE)
  utils::getS3method(generic, "maxLik", envir = asNamespace("maxLik"))
}
maxlik_saved_coef <- function(object, ...) maxlik_saved_method("coef")(object, ...)
maxlik_saved_logLik <- function(object, ...) maxlik_saved_method("logLik")(object, ...)
maxlik_saved_vcov <- function(object, ...) maxlik_saved_method("vcov")(object, ...)
maxlik_saved_summary <- function(object, ...) maxlik_saved_method("summary")(object, ...)
maxlik_saved_print <- function(x, ...) maxlik_saved_method("print")(x, ...)

register_maxlik_saved_methods <- function() {
  if (isNamespaceLoaded("maxLik")) return(invisible(NULL))
  bridges <- list(coef = maxlik_saved_coef, logLik = maxlik_saved_logLik,
    vcov = maxlik_saved_vcov, summary = maxlik_saved_summary, print = maxlik_saved_print)
  for (generic in names(bridges)) {
    namespace <- asNamespace(if (generic %in% c("summary", "print")) "base" else "stats")
    if (is.null(utils::getS3method(generic, "maxLik", optional = TRUE, envir = namespace)))
      registerS3method(generic, "maxLik", bridges[[generic]], envir = namespace)
  }
  invisible(NULL)
}

stats_fit_active <- function(start, fixed) {
  active <- rep(TRUE, length(start))
  if (is.null(fixed)) return(active)
  valid <- is.null(dim(fixed)) && !anyNA(fixed)
  if (is.logical(fixed)) {
    valid <- valid && length(fixed) == length(start)
    if (valid) active <- !fixed
  } else if (is.numeric(fixed)) {
    valid <- valid && all(is.finite(fixed) & fixed == floor(fixed) &
                           fixed >= 1 & fixed <= length(start))
    if (valid) active[fixed] <- FALSE
  } else if (is.character(fixed)) {
    valid <- valid && all(fixed %in% names(start))
    if (valid) active[match(fixed, names(start))] <- FALSE
  } else valid <- FALSE
  if (!valid) stop("For optimizer = 'stats', fixed must be a full logical vector, parameter positions, or parameter names.", call. = FALSE)
  active
}

stats_fit_constraints <- function(constraints, start, active) {
  if (is.null(constraints)) return(NULL)
  if (!is.list(constraints) || length(constraints) != 2L ||
      is.null(names(constraints)) || anyDuplicated(names(constraints)) ||
      !setequal(names(constraints), c("ineqA", "ineqB"))) {
    stop("optimizer = 'stats' supports constraints = list(ineqA, ineqB); equality constraints require optimizer = 'maxLik'.", call. = FALSE)
  }
  A <- constraints$ineqA
  B <- constraints$ineqB
  if (!is.matrix(A) || !is.numeric(A) || ncol(A) != length(start) ||
      !is.numeric(B) || !is.null(dim(B)) || length(B) != nrow(A) ||
      any(!is.finite(A)) || any(!is.finite(B))) {
    stop("ineqA must be a finite numeric matrix with one column per optimization parameter; ineqB must be a finite vector with one entry per row.", call. = FALSE)
  }
  offset <- drop(A[, !active, drop = FALSE] %*% start[!active]) + B
  ui <- A[, active, drop = FALSE]
  movable <- rowSums(abs(ui)) > 0
  if (any(offset[!movable] < 0)) stop("Fixed parameters violate a linear inequality constraint.", call. = FALSE)
  ui <- ui[movable, , drop = FALSE]
  ci <- -offset[movable]
  if (length(ci) && any(drop(ui %*% start[active]) - ci <= 0)) {
    stop("constrOptim() requires a starting point strictly inside every movable linear inequality constraint.", call. = FALSE)
  }
  list(ui = ui, ci = ci)
}

stats_fit_derivative <- function(value, start, hessian = FALSE) {
  p <- length(start)
  if (!hessian && is.matrix(value) && ncol(value) == p) value <- colSums(value)
  if (hessian) {
    valid <- is.numeric(value) && is.matrix(value) && identical(dim(value), c(p, p))
    if (valid && (!is.null(rownames(value)) || !is.null(colnames(value)))) {
      valid <- !is.null(rownames(value)) && !is.null(colnames(value)) &&
        !anyDuplicated(rownames(value)) && !anyDuplicated(colnames(value)) &&
        setequal(rownames(value), names(start)) && setequal(colnames(value), names(start))
      if (valid) value <- value[names(start), names(start), drop = FALSE]
    }
  } else {
    valid <- is.numeric(value) && is.null(dim(value)) && length(value) == p
    if (valid && anyNA(names(value))) valid <- FALSE
    if (valid && !is.null(names(value)) && all(!is.na(names(value)) & nzchar(names(value)))) {
      valid <- !anyNA(names(value)) && !anyDuplicated(names(value)) &&
        setequal(names(value), names(start))
      if (valid) value <- value[names(start)]
    }
  }
  if (!valid || any(!is.finite(value))) {
    stop(if (hessian) "The stats optimizer requires a finite full-parameter Hessian matrix with matching dimensions and names." else
      "The stats optimizer requires a finite full-parameter score vector (or observation-by-parameter score matrix) with matching names.", call. = FALSE)
  }
  value
}

run_stats_fit <- function(fn, start, grad, hess, constraints, method, additional,
                          stable) {
  allowed <- c("fixed", "control", "finalHessian", "lower", "upper",
               "outer.iterations", "outer.eps")
  if (length(additional) && (is.null(names(additional)) || anyNA(names(additional)) ||
      any(names(additional) == "") || anyDuplicated(names(additional)) ||
      any(!names(additional) %in% allowed))) {
    stop("For optimizer = 'stats', ... accepts only fixed, control, finalHessian, lower, upper, outer.iterations and outer.eps; use control = list(maxit = ...) instead of iterlim.", call. = FALSE)
  }
  method <- match.arg(method, c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B"))
  control <- additional$control
  if (is.null(control)) control <- list()
  controls <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol",
                "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% controls))))
    stop("Use named stats::optim() control entries for optimizer = 'stats' (for example maxit, not iterlim).", call. = FALSE)
  if (!is.null(control$fnscale) && (!is.numeric(control$fnscale) ||
      length(control$fnscale) != 1L || !is.finite(control$fnscale) || control$fnscale <= 0))
    stop("control$fnscale must be positive: the stats backend minimizes negative log likelihood.", call. = FALSE)
  final <- additional$finalHessian
  if (is.null(final)) final <- TRUE
  if (length(final) != 1L || is.na(final) || !(is.logical(final) || is.numeric(final)) ||
      !final %in% c(0, 1)) stop("For optimizer = 'stats', finalHessian must be TRUE or FALSE.", call. = FALSE)
  for (name in c("outer.iterations", "outer.eps")) if (!is.null(additional[[name]])) {
    value <- additional[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0 ||
        (name == "outer.iterations" && value != floor(value)))
      stop(name, " must be a positive finite ", if (name == "outer.iterations") "integer." else "number.", call. = FALSE)
  }
  active <- stats_fit_active(start, additional$fixed)
  free <- start[active]
  linear <- stats_fit_constraints(constraints, start, active)
  constrained <- !is.null(linear) && nrow(linear$ui) > 0
  if (constrained && !method %in% c("BFGS", "Nelder-Mead"))
    stop("Linear inequalities with optimizer = 'stats' require method = 'BFGS' or 'Nelder-Mead'.", call. = FALSE)
  bounds <- !is.null(additional$lower) || !is.null(additional$upper)
  if (bounds && (method != "L-BFGS-B" || !is.null(constraints)))
    stop("lower/upper require method = 'L-BFGS-B' and cannot be combined with linear constraints.", call. = FALSE)
  bound <- function(x, default) {
    if (is.null(x)) return(rep(default, length(start)))
    if (!is.numeric(x) || !is.null(dim(x)) || anyNA(x) ||
        !length(x) %in% c(1L, length(start)))
      stop("lower/upper must be numeric scalars or vectors for all optimization parameters.", call. = FALSE)
    if (!is.null(names(x))) {
      if (length(x) != length(start) || anyDuplicated(names(x)) ||
          !setequal(names(x), names(start))) stop("Named bounds must match all parameter names.", call. = FALSE)
      x <- x[names(start)]
    }
    rep(x, length.out = length(start))
  }
  lower <- bound(additional$lower, -Inf)
  upper <- bound(additional$upper, Inf)
  if (any(lower > upper) || any(start < lower | start > upper) ||
      any(active & lower == upper))
    stop("The starting point must satisfy lower/upper; use fixed for parameters with equal bounds.", call. = FALSE)
  # Vector-valued optim controls refer to the free coordinates after reduction.
  for (name in c("parscale", "ndeps")) if (!is.null(control[[name]]) &&
      length(control[[name]]) != length(free))
    stop("control$", name, " must have one entry per free parameter.", call. = FALSE)
  expand <- function(x) { p <- start; p[active] <- x; p }
  initial <- fn(start)
  total <- function(value) {
    if (!is.numeric(value) || !length(value) || anyNA(value) || any(value == Inf))
      stop("The likelihood must return numeric log-likelihood contributions without NA or positive infinity.", call. = FALSE)
    sum(value)
  }
  if (!is.finite(total(initial))) stop("The starting log likelihood must be finite.", call. = FALSE)
  gradient_attribute <- is.null(grad) && !is.null(attr(initial, "gradient", exact = TRUE))
  score <- if (!is.null(grad)) function(p) stats_fit_derivative(grad(p), start) else if (
    gradient_attribute) function(p) stats_fit_derivative(attr(fn(p), "gradient", exact = TRUE), start) else NULL
  objective <- function(x) -total(fn(expand(x)))
  numeric_score <- function(x) {
    step <- .Machine$double.eps^(1/3) * pmax(abs(x), 1)
    result <- numeric(length(x))
    for (j in seq_along(x)) {
      above <- below <- x
      h <- step[j]
      if (constrained) {
        # constrOptim's barrier is defined strictly inside ui %*% x > ci.
        # Keep central differences there too, including near active boundaries.
        direction <- abs(linear$ui[, j])
        moving <- direction > 0
        if (any(moving)) {
          slack <- drop(linear$ui %*% x) - linear$ci
          h <- min(h, .5 * min(slack[moving] / direction[moving]))
        }
      }
      above[j] <- above[j] + h; below[j] <- below[j] - h
      # Stay inside box bounds when numerical derivatives are needed at a bound.
      above[j] <- min(above[j], upper[active][j])
      below[j] <- max(below[j], lower[active][j])
      if (above[j] == below[j])
        stop("Numerical differentiation has no representable feasible step; supply an analytic gradient or rescale the parameters.", call. = FALSE)
      result[j] <- (objective(above) - objective(below)) / (above[j] - below[j])
    }
    result
  }
  negative_score <- if (!is.null(score)) function(x) -score(expand(x))[active] else NULL
  args <- list(par = free, fn = objective, gr = negative_score, method = method, control = control)
  if (constrained) {
    args$theta <- args$par; args$par <- NULL
    args$f <- args$fn; args$fn <- NULL
    args$grad <- if (is.null(negative_score) && method == "BFGS") numeric_score else negative_score
    args$gr <- NULL
    args$ui <- linear$ui; args$ci <- linear$ci
    for (name in c("outer.iterations", "outer.eps")) if (!is.null(additional[[name]]))
      args[[name]] <- additional[[name]]
  } else {
    if (any(c("outer.iterations", "outer.eps") %in% names(additional)))
      stop("outer.iterations/outer.eps require movable linear inequality constraints.", call. = FALSE)
    if (method == "L-BFGS-B") { args$lower <- lower[active]; args$upper <- upper[active] }
  }
  result <- if (!length(free)) list(par = free, convergence = 0L,
      message = "All parameters fixed; likelihood evaluated without optimization.", counts = c(functions = 0L)) else
    do.call(if (constrained) stats::constrOptim else stats::optim, args)
  estimate <- expand(result$par)
  value <- fn(estimate)
  maximum <- total(value)
  gradient <- stats::setNames(numeric(length(start)), names(start))
  if (length(free)) {
    if (!is.null(score)) gradient <- score(estimate) else
      gradient[active] <- -numeric_score(result$par)
  }
  H <- NULL
  hessian_source <- "disabled"
  if (isTRUE(as.logical(final)) && length(free) && (!stable || !is.null(hess))) {
    H <- matrix(0, length(start), length(start), dimnames = list(names(start), names(start)))
    tagged <- attr(value, "hessian", exact = TRUE)
    if (!is.null(hess)) {
      H <- stats_fit_derivative(hess(estimate), start, TRUE)
      hessian_source <- "user"
    } else if (!is.null(tagged)) {
      H <- stats_fit_derivative(tagged, start, TRUE)
      hessian_source <- "attribute"
    } else {
      curvature <- tryCatch(-stats::optimHess(result$par, objective, negative_score,
        control = control[intersect(names(control), c("fnscale", "parscale", "ndeps"))]),
        error = function(e) {
          warning("Final stats Hessian is unavailable: ", conditionMessage(e), call. = FALSE)
          NULL
        })
      if (!is.null(curvature)) H[active, active] <- curvature else H <- NULL
      hessian_source <- if (is.null(H)) "unavailable" else "stats_numeric"
    }
  }
  # Box constraints use the same boundary diagnostic as linear inequalities.
  curvature_constraints <- constraints
  if (bounds) {
    indices <- which(is.finite(lower) | is.finite(upper))
    A <- matrix(0, 0L, length(start)); B <- numeric(0)
    for (i in indices) {
      row <- numeric(length(start)); row[i] <- 1
      if (is.finite(lower[i])) { A <- rbind(A, row); B <- c(B, -lower[i]) }
      if (is.finite(upper[i])) { A <- rbind(A, -row); B <- c(B, upper[i]) }
    }
    if (nrow(A)) curvature_constraints <- list(ineqA = A, ineqB = B)
  }
  covariance <- matrix(0, length(start), length(start), dimnames = list(names(start), names(start)))
  covariance[active, active] <- NA_real_
  boundary <- if (length(free)) fit_constraint_status(curvature_constraints, estimate, active,
                                                     pmax(abs(estimate[active]), .1)) else NULL
  if (!is.null(H) && is.null(boundary)) {
    information <- fit_information(H[active, active, drop = FALSE])
    if (information$positive) covariance[active, active] <- information$covariance
  }
  fit <- list(estimate = estimate, maximum = maximum, gradient = gradient, hessian = H,
    varcovar = covariance, code = result$convergence,
    message = if (is.null(result$message)) if (result$convergence == 0L) "Normal convergence." else
      "The stats optimizer did not converge; inspect convergence code and control settings." else result$message,
    counts = result$counts, active = active, objectiveFn = fn,
    method = method, optimizer = if (constrained) "constrOptim" else "optim",
    curvature_constraints = curvature_constraints, boundary = boundary,
    hessian_source = hessian_source, gradient_source = if (!is.null(grad)) "user" else
      if (gradient_attribute) "attribute" else "numeric",
    lower = lower, upper = upper, control = control)
  class(fit) <- "emhawkes_stats_fit"
  fit
}

#' @export
coef.emhawkes_stats_fit <- function(object, ...) object$estimate
#' @method logLik emhawkes_stats_fit
#' @export
logLik.emhawkes_stats_fit <- function(object, ...) {
  structure(object$maximum, class = "logLik", df = sum(object$active))
}
#' @export
vcov.emhawkes_stats_fit <- function(object, ...) object$varcovar
#' @export
summary.emhawkes_stats_fit <- function(object, ...) {
  estimates <- stats::coef(object)
  errors <- sqrt(diag(stats::vcov(object)))
  z <- rep(NA_real_, length(estimates)); names(z) <- names(estimates)
  usable <- object$active & is.finite(errors) & errors > 0
  z[usable] <- estimates[usable] / errors[usable]
  structure(list(coefficients = cbind(Estimate = estimates, `Std. Error` = errors,
    `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))),
    maximum = object$maximum, code = object$code, message = object$message,
    optimizer = object$optimizer, method = object$method, fixed = !object$active),
    class = "summary.emhawkes_stats_fit")
}
#' @export
print.emhawkes_stats_fit <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
#' @export
print.summary.emhawkes_stats_fit <- function(x, ...) {
  cat("Hawkes maximum likelihood fit (stats::", x$optimizer, ", ", x$method, ")\n", sep = "")
  cat("Log likelihood: ", format(x$maximum), "\n", sep = "")
  cat("Convergence code: ", x$code, " (", x$message, ")\n", sep = "")
  stats::printCoefmat(x$coefficients, ...)
  if (any(x$fixed)) cat("Fixed parameters: ", paste(rownames(x$coefficients)[x$fixed], collapse = ", "), "\n", sep = "")
  invisible(x)
}
