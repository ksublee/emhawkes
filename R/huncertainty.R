#' @include hinterface.R
NULL

# Differentiate in scaled coordinates, centering values before subtraction.
fit_scaled_derivatives <- function(fn, parameters, scale, step, value) {
  count <- length(parameters)
  gradient <- stats::setNames(numeric(count), names(parameters))
  hessian <- matrix(0, count, count,
                    dimnames = list(names(parameters), names(parameters)))
  evaluate <- function(offset) {
    answer <- tryCatch(sum(fn(parameters + scale * offset)),
                       error = function(e) NA_real_)
    if (length(answer) != 1L || !is.finite(answer)) return(NA_real_)
    answer - value
  }
  for (i in seq_len(count)) {
    delta_i <- numeric(count)
    delta_i[i] <- step
    above <- evaluate(delta_i)
    below <- evaluate(-delta_i)
    gradient[i] <- (above - below) / (2 * step)
    hessian[i, i] <- (above + below) / step^2
    if (i < count) for (j in (i + 1L):count) {
      delta_j <- numeric(count)
      delta_j[j] <- step
      hessian[i, j] <- hessian[j, i] <-
        ((evaluate(delta_i + delta_j) - evaluate(delta_i - delta_j)) -
         (evaluate(-delta_i + delta_j) - evaluate(-delta_i - delta_j))) /
        (4 * step^2)
    }
  }
  list(gradient = gradient, hessian = hessian)
}

# Differentiate the full analytic score, retaining only free coordinates.
fit_scaled_score <- function(score, parameters, active, scale, step, gradient) {
  free <- parameters[active]
  evaluate <- function(offset) {
    proposal <- parameters
    proposal[active] <- free + scale * offset
    answer <- score(proposal)
    if (!is.numeric(answer) || length(answer) != length(parameters) ||
        any(!is.finite(answer))) stop("Analytic curvature score is not finite or has the wrong length.")
    if (!is.null(names(answer))) {
      if (anyNA(names(answer)) || anyDuplicated(names(answer)) ||
          !setequal(names(answer), names(parameters)))
        stop("Analytic curvature score names do not match the fitted parameters.")
      answer <- answer[names(parameters)]
    }
    answer[active] * scale
  }
  if (is.null(gradient)) gradient <- evaluate(numeric(length(free)))
  H <- matrix(0, length(free), length(free), dimnames = list(names(free), names(free)))
  for (j in seq_along(free)) {
    offset <- numeric(length(free))
    offset[j] <- step
    H[, j] <- (evaluate(offset) - evaluate(-offset)) / (2 * step)
  }
  list(gradient = gradient, hessian = (H + t(H)) / 2)
}

fit_information <- function(hessian) {
  information <- -(hessian + t(hessian)) / 2
  if (any(!is.finite(information))) {
    return(list(positive = FALSE, reason = "nonfinite", minimum = NA_real_,
                condition = NA_real_, covariance = NULL))
  }
  values <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  positive <- min(values) > 1e-10 * max(abs(values))
  covariance <- if (positive) tryCatch(chol2inv(chol(information)),
                                       error = function(e) NULL) else NULL
  list(positive = positive && !is.null(covariance),
       reason = if (min(values) <= 0) "not_positive_definite" else "rank_deficient",
       minimum = min(values),
       condition = if (positive) max(values) / min(values) else Inf,
       covariance = covariance)
}

# Compare inverse curvature as well as Hessians, including weak directions.
fit_curvature_change <- function(first, second) {
  if (!first$information$positive || !second$information$positive) {
    return(c(hessian = Inf, covariance = Inf, standard_error = Inf))
  }
  C1 <- first$information$covariance
  C2 <- second$information$covariance
  se1 <- sqrt(diag(C1))
  se2 <- sqrt(diag(C2))
  c(hessian = norm(first$hessian - second$hessian, "F") /
      max(norm(second$hessian, "F"), .Machine$double.eps),
    covariance = norm(C1 - C2, "F") / max(norm(C2, "F"), .Machine$double.eps),
    standard_error = max(abs(se1 - se2) / pmax(se1, se2)))
}

stable_fit_curvature <- function(fit, score = NULL) {
  parameters <- stats::coef(fit)
  active <- fit_active(fit)
  free <- parameters[active]
  scale <- pmax(abs(free), .1)
  objective <- fit_objective(fit)
  fn <- function(candidate) {
    proposal <- parameters
    proposal[active] <- candidate
    sum(objective(proposal))
  }
  value <- tryCatch(fn(free), error = function(e) NA_real_)
  if (!is.finite(value) || !isTRUE(all.equal(value, as.numeric(stats::logLik(fit)),
                                           tolerance = 1e-8))) {
    return(list(status = "unavailable", reason = "objective_changed", scale = scale))
  }
  steps <- .004 / 2^(0:5)
  derivatives <- list()
  extrapolated <- list()
  checks <- list()
  chosen <- NULL
  change <- NULL
  score_gradient <- NULL
  for (index in seq_along(steps)) {
    derivative <- if (is.null(score)) {
      fit_scaled_derivatives(fn, free, scale, steps[index], value)
    } else {
      fit_scaled_score(score, parameters, active, scale, steps[index], score_gradient)
    }
    if (!is.null(score)) score_gradient <- derivative$gradient
    derivative$information <- fit_information(derivative$hessian)
    derivatives[[index]] <- derivative
    checks[[index]] <- data.frame(step = steps[index],
      minimum_scaled_information_eigenvalue = derivative$information$minimum,
      scaled_condition = derivative$information$condition,
      positive_definite = derivative$information$positive)
    if (index > 1L) {
      previous <- derivatives[[index - 1L]]
      extrapolation <- list(
        hessian = (4 * derivative$hessian - previous$hessian) / 3,
        gradient = (4 * derivative$gradient - previous$gradient) / 3)
      extrapolation$information <- fit_information(extrapolation$hessian)
      extrapolated[[index - 1L]] <- extrapolation
    }
    if (index > 2L) {
      change <- fit_curvature_change(extrapolated[[index - 2L]],
                                     extrapolated[[index - 1L]])
      central_change <- fit_curvature_change(derivatives[[index - 1L]], derivative)
      # Do not select a positive matrix from otherwise unstable steps.
      positive_steps <- all(vapply(derivatives[(index - 2L):index],
        function(item) item$information$positive, logical(1)))
      if (positive_steps && change["hessian"] <= .01 &&
          change["covariance"] <= .02 && change["standard_error"] <= .01 &&
          central_change["standard_error"] <= .02) {
        chosen <- extrapolated[[index - 1L]]
        break
      }
    }
  }
  common <- list(curvature_source = if (is.null(score)) "likelihood" else "analytic_score",
                 scale = scale, checks = do.call(rbind, checks),
                 change = change, original_hessian = fit_hessian(fit))
  if (is.null(chosen)) {
    reason <- if (any(vapply(derivatives, function(item) item$information$positive,
                             logical(1)))) "step_sensitive" else
      derivatives[[length(derivatives)]]$information$reason
    return(c(list(status = "unavailable", reason = reason), common))
  }
  information <- chosen$information
  covariance <- information$covariance * outer(scale, scale)
  hessian <- chosen$hessian / outer(scale, scale)
  gradient <- chosen$gradient / scale
  decrement <- sqrt(max(0, drop(crossprod(chosen$gradient,
                          information$covariance %*% chosen$gradient))))
  c(list(status = "ok", reason = NULL,
         method = if (is.null(score)) "scaled_central_Richardson" else "scaled_score_Richardson",
         hessian = hessian, covariance = covariance, gradient = gradient,
         scaled_condition = information$condition, newton_decrement = decrement,
         maximum_standardized_score = max(abs(gradient) * sqrt(diag(covariance)))), common)
}

fit_constraint_status <- function(constraints, parameters, active, scale) {
  if (is.null(constraints)) return(NULL)
  if (!is.list(constraints) ||
      !setequal(names(constraints), c("ineqA", "ineqB"))) return("unsupported_constraints")
  A <- as.matrix(constraints$ineqA)
  B <- constraints$ineqB
  if (!is.numeric(A) || ncol(A) != length(parameters) || length(B) != nrow(A) ||
      any(!is.finite(A)) || any(!is.finite(B))) return("unsupported_constraints")
  scaled_A <- sweep(A[, active, drop = FALSE], 2L, scale, "*")
  distance <- drop(A %*% parameters + B)
  # Mixed differences move at most two coordinates by the largest step.
  radius <- .004 * apply(abs(scaled_A), 1L, function(row) {
    sum(sort(row, decreasing = TRUE)[seq_len(min(2L, length(row)))])
  })
  movable <- rowSums(abs(scaled_A)) > 0
  if (any(distance[movable] <= radius[movable])) "constraint_near_boundary" else NULL
}

finalize_hawkes_fit <- function(fit, object, constraints, hess, additional, method,
                               final_hessian, score_factory = NULL) {
  unchecked <- function(reason) {
    list(fit = fit, diagnostics = list(status = "not_checked", reason = reason))
  }
  if (final_hessian == "optimizer" || (final_hessian == "maxLik" && inherits(fit, "emhawkes_stats_fit"))) {
    if (inherits(fit, "emhawkes_stats_fit")) {
      if (!is.null(fit$boundary)) {
        warning("Wald standard errors are unavailable: the estimate is near a constraint boundary. See fit_diagnostics().", call. = FALSE)
        return(list(fit = fit, diagnostics = list(status = "unavailable", reason = fit$boundary)))
      }
      return(unchecked("optimizer_requested"))
    }
    final_hessian <- "maxLik"
  }
  if (final_hessian == "maxLik") {
    # Inspect the stored curvature without re-evaluating the likelihood.
    final <- additional[["finalHessian"]]
    bhhh <- (is.character(final) && length(final) == 1L && tolower(final) == "bhhh") ||
      identical(additional[["bhhhHessian"]], TRUE) || identical(toupper(method), "BHHH")
    if (is.null(object@dresidual) && is.null(hess) && !bhhh) {
      standard_errors <- tryCatch(suppressWarnings(sqrt(diag(stats::vcov(fit)))),
                                  error = function(e) numeric(0))
      if (any(is.nan(standard_errors))) {
        message('hfit: some standard errors are NaN. Try final_hessian = "stable" ',
                'for Hawkes Hessian diagnostics, then inspect fit_diagnostics().')
      }
    }
    return(unchecked("maxLik_requested"))
  }
  if (!is.null(object@dresidual)) return(unchecked("residual_likelihood"))
  final <- additional[["finalHessian"]]
  if (!is.null(final) && length(final) == 1L && !is.na(final) &&
      (is.logical(final) || is.numeric(final)) && !as.logical(final)) {
    return(unchecked("finalHessian_false"))
  }
  if ((is.character(final) && length(final) == 1L && tolower(final) == "bhhh") ||
      identical(additional[["bhhhHessian"]], TRUE) ||
      identical(toupper(method), "BHHH")) return(unchecked("BHHH_requested"))
  parameters <- stats::coef(fit)
  active <- fit_active(fit)
  if (!is.numeric(parameters) || any(!is.finite(parameters)) || !any(active) ||
      !is.finite(as.numeric(stats::logLik(fit)))) return(unchecked("fit_unavailable"))
  scale <- pmax(abs(parameters[active]), .1)
  boundary <- fit_constraint_status(constraints, parameters, active, scale)
  if (!is.null(boundary)) {
    diagnostics <- list(status = "unavailable", reason = boundary, scale = scale)
  } else if (!is.null(hess)) {
    H <- fit_hessian(fit)[active, active, drop = FALSE]
    information <- fit_information(H * outer(scale, scale))
    diagnostics <- list(status = if (information$positive) "ok" else "unavailable",
      reason = if (information$positive) NULL else information$reason,
      method = "supplied", scale = scale, hessian = H,
      covariance = if (information$positive) information$covariance * outer(scale, scale) else NULL,
      scaled_condition = information$condition)
  } else {
    score_reason <- NULL
    diagnostics <- if (!is.null(score_factory)) tryCatch({
      objective <- score_factory()
      value <- sum(objective$loglik(parameters))
      if (!is.finite(value) || !isTRUE(all.equal(value, as.numeric(stats::logLik(fit)),
                                                tolerance = 1e-8)))
        stop("Analytic curvature likelihood does not match the fitted likelihood.")
      score <- if (is.function(objective$score)) objective$score else objective$gradient
      stable_fit_curvature(fit, score = score)
    }, error = function(e) {
      score_reason <<- conditionMessage(e)
      NULL
    }) else NULL
    if (is.null(diagnostics)) diagnostics <- stable_fit_curvature(fit)
    diagnostics$score_fallback_reason <- score_reason
  }
  if (identical(diagnostics$status, "ok") && identical(diagnostics$method, "supplied")) {
    gradient <- fit_gradient(fit)[active]
    diagnostics$gradient <- gradient
    usable_gradient <- is.numeric(gradient) && length(gradient) == sum(active) &&
      all(is.finite(gradient))
    diagnostics$newton_decrement <- if (usable_gradient) {
      sqrt(max(0, drop(crossprod(gradient, diagnostics$covariance %*% gradient))))
    } else NA_real_
    diagnostics$maximum_standardized_score <- if (usable_gradient)
      max(abs(gradient) * sqrt(diag(diagnostics$covariance))) else NA_real_
  }
  covariance <- matrix(0, length(parameters), length(parameters),
                        dimnames = list(names(parameters), names(parameters)))
  diagnostics$optimizer_code <- fit_code(fit)
  diagnostics$optimizer_message <- fit_message(fit)
  numerical_methods <- c("nr", "newton-raphson", "bfgsr", "bfgs-r", "bhhh")
  diagnostics$optimizer_success <- if (tolower(method) == "sann") NA else if (
    tolower(method) %in% numerical_methods) diagnostics$optimizer_code %in% c(1L, 2L, 8L) else
      diagnostics$optimizer_code == 0L
  diagnostics$free_parameters <- names(parameters)[active]
  diagnostics$warnings <- character(0)
  if (!is.na(diagnostics$optimizer_success) && !diagnostics$optimizer_success) {
    diagnostics$warnings <- c(diagnostics$warnings, "optimizer_not_converged")
    warning("The optimizer did not report normal convergence (code ",
            diagnostics$optimizer_code, "). Stable curvature does not establish a reliable fit; ",
            "see fit_diagnostics().", call. = FALSE)
  }
  if (identical(diagnostics$status, "ok")) {
    covariance[active, active] <- diagnostics$covariance
    if (diagnostics$method %in% c("scaled_central_Richardson", "scaled_score_Richardson")) {
      H <- matrix(0, length(parameters), length(parameters), dimnames = dimnames(covariance))
      H[active, active] <- diagnostics$hessian
      fit$hessian <- H
    }
    if (isTRUE(diagnostics$newton_decrement > .1)) {
      diagnostics$warnings <- c(diagnostics$warnings, "score_not_small")
      warning("The score is not close to zero. ",
              "Consider tighter optimization and multiple starting values; see fit_diagnostics().",
              call. = FALSE)
    }
    if (diagnostics$scaled_condition > 1e8) {
      diagnostics$warnings <- c(diagnostics$warnings, "weak_curvature")
      warning("The information matrix has weak curvature. Wald standard errors may be unreliable; ",
              "see fit_diagnostics().", call. = FALSE)
    }
  } else {
    covariance[active, active] <- NA_real_
    explanation <- switch(diagnostics$reason,
      nonfinite = "the curvature contains nonfinite values",
      not_positive_definite = "the information matrix is not positive definite",
      rank_deficient = "the information matrix has insufficient numerical rank",
      step_sensitive = "the curvature changes with the differentiation step",
      objective_changed = "the fitted likelihood could not be reproduced",
      unsupported_constraints = "the supplied constraints are not supported by this curvature check",
      constraint_near_boundary = "the estimate is near a constraint boundary",
      "the curvature check failed")
    warning("Wald standard errors are unavailable: ", explanation, " (", diagnostics$reason,
            "). See fit_diagnostics(); consider profile likelihood or bootstrap.", call. = FALSE)
  }
  fit$varcovar <- covariance
  list(fit = fit, diagnostics = diagnostics)
}

# Summarize recorded execution without evaluating the model again.
fit_computation <- function(fit, execution, gradient_requested, gradient_used,
                            final_hessian, hess, additional, diagnostics,
                            maxLik_hessian_skipped = FALSE) {
  gradient_reason <- NULL
  if (isTRUE(execution$gradient_attribute)) {
    if (gradient_used == "user") {
      gradient_used <- NA_character_
      gradient_reason <- "Both grad and likelihood gradient attributes were supplied; maxLik controls precedence."
    } else gradient_used <- "attribute"
  }
  H <- fit_hessian(fit)
  final <- additional[["finalHessian"]]
  disabled <- !is.null(final) && length(final) == 1L && !is.na(final) &&
    (is.logical(final) || is.numeric(final)) && !as.logical(final)
  reason <- if (final_hessian == "stable" && !identical(diagnostics$status, "ok"))
    diagnostics$reason else NULL
  used <- if (is.null(H)) {
    if (disabled) "disabled" else "unavailable"
  } else if (isTRUE(diagnostics$method %in% c("scaled_central_Richardson", "scaled_score_Richardson")) &&
             identical(diagnostics$status, "ok")) {
    diagnostics$method
  } else if (identical(attr(H, "type", exact = TRUE), "BHHH")) {
    "BHHH"
  } else if (isTRUE(execution$hessian_attribute)) {
    if (is.null(hess)) "attribute" else NA_character_
  } else if (!is.null(hess)) "user" else "maxLik_numeric"
  if (inherits(fit, "emhawkes_stats_fit")) {
    gradient_used <- if (identical(gradient_used, "analytic")) "analytic" else fit$gradient_source
    gradient_reason <- NULL
    if (!isTRUE(diagnostics$method %in% c("scaled_central_Richardson", "scaled_score_Richardson", "supplied"))) {
      used <- fit$hessian_source
      if (is.null(H) && !disabled && final_hessian == "stable") used <- "unavailable"
    }
  }
  if (is.na(used)) {
    reason <- "Both hess and likelihood Hessian attributes were supplied; maxLik controls precedence."
  } else if (used == "unavailable" && is.null(reason)) {
    reason <- if (inherits(fit, "emhawkes_stats_fit")) "The stats provider did not return a final Hessian." else
      "maxLik did not return a final Hessian."
  }
  list(recorded = TRUE,
    engine = list(requested = execution$requested, used = execution$used,
      fallback_reason = execution$fallback_reason),
    gradient = list(requested = gradient_requested, used = gradient_used,
      reason = gradient_reason),
    hessian = list(requested = final_hessian, finalHessian = final,
      used = used, reason = reason, maxLik_skipped = maxLik_hessian_skipped))
}

unrecorded_fit_computation <- function(metadata) {
  list(recorded = FALSE,
    engine = list(requested = NA_character_, used = NA_character_, fallback_reason = NULL),
    gradient = list(requested = NA_character_,
      used = if (is.null(metadata$gradient)) NA_character_ else metadata$gradient, reason = NULL),
    hessian = list(requested = NA_character_, finalHessian = NULL,
      used = NA_character_, reason = "not_recorded", maxLik_skipped = NA))
}

#' Inspect Fit Computation and Curvature Diagnostics
#'
#' Return the recorded calculation methods and final Hessian checks from [hfit()].
#' @param fit A result from [hfit()].
#' @return A list with recorded methods in `computation`, calculation status,
#' reason, free parameters, curvature source and analytic-score fallback reason,
#' scale, and, when available, difference-step checks, curvature changes,
#' recalculated score, scaled condition number, and Newton decrement.
#' @details
#' `computation` records the engine, optimizer gradient source and stored final
#' Hessian source without re-evaluating the model. Its `engine` list contains
#' `requested`, `used` and `fallback_reason`. `used` lists engines actually invoked
#' in first-use order, including any partial C++ pass before fallback. It can
#' contain both `"cpp"` and `"r"`. A custom likelihood uses `"custom"` and has no
#' requested built-in engine. Recording includes final curvature evaluations.
#'
#' `computation$gradient$used` is `"numeric"`, `"analytic"`, `"user"` or
#' `"attribute"`. It describes the returned optimizer score, separately from a recalculated
#' score in stable curvature checks. `computation$hessian$used` is
#' `"maxLik_numeric"`, `"scaled_central_Richardson"`, `"scaled_score_Richardson"`,
#' `"stats_numeric"`, `"user"`, `"attribute"`,
#' `"BHHH"`, `"disabled"` or `"unavailable"`. The Hessian list also retains the
#' requested `final_hessian` choice and any `finalHessian` override. Its `reason`
#' explains a skipped stabilization or unavailable result. `maxLik_skipped` records
#' whether stable curvature omitted maxLik's separate final Hessian. In that case
#' `original_hessian` is `NULL`. Numerical stable checks otherwise retain the
#' matrix returned by maxLik.
#' Older computation records may omit `maxLik_skipped`. If both a derivative
#' function and likelihood derivative attributes were supplied, the source is
#' `NA` with a reason because precedence depends on maxLik.
#'
#' Stable curvature records `curvature_source = "analytic_score"` when it
#' differentiates the built-in analytic score, including with a numerical optimizer
#' gradient. On success its method is `"scaled_score_Richardson"`. Likelihood
#' differences record `curvature_source = "likelihood"` and, on success,
#' `"scaled_central_Richardson"`. A score matrix that fails the step checks remains
#' unavailable rather than switching methods to obtain standard errors.
#' `score_fallback_reason` records an analytic-adapter or score-evaluation failure;
#' an explicitly requested R engine uses likelihood differences without trying
#' the analytic score. Both paths retain the same checks across difference steps.
#'
#' The saved record is a snapshot; later objective evaluations or option changes
#' do not update it. Older fits with model metadata but no computation record
#' return `recorded = FALSE` and `NA` for unrecorded methods. Any previously saved
#' gradient selection is retained. Their engine is not inferred from current options.
#'
#' `status = "ok"` means the numerical curvature checks passed, not that the
#' optimizer found a global maximum or that Wald intervals have correct coverage.
#' The numerical calculation checks information and standard-error stability
#' across steps. The stored optimizer gradient and coefficients are not changed.
#' For checked fits, optimizer nonconvergence is warned even when curvature is
#' unavailable. Inspect `optimizer_success`, `optimizer_code` and the additional
#' warning codes in `warnings` separately from the curvature `status` and `reason`.
#' A Newton decrement above 0.1 triggers a score warning. Its square divided by
#' two estimates the local quadratic improvement in log-likelihood.
#' An information condition number above `1e8` triggers a weak-curvature warning;
#' eigenvalues at or below `1e-10` times the largest absolute eigenvalue do not
#' support numerical inversion. These thresholds are numerical diagnostics.
#' With the default `final_hessian = "maxLik"`, these additional checks are not
#' performed and the status is `"not_checked"`. Request `final_hessian = "stable"`
#' to run them. Earlier fits without these checks also have status `"not_checked"`.
#' New computation records include `optimizer` with the requested provider,
#' actual R optimizer and method. The stats provider defaults to stable curvature;
#' `final_hessian = "optimizer"` retains its numerical or supplied Hessian.
#' @seealso [hfit()], [fitted_hspec()]
#' @export
#' @examplesIf requireNamespace("maxLik", quietly = TRUE)
#' h <- new("hspec", mu = .3, alpha = 0, beta = 1)
#' fit <- hfit(h, c(0, .2, .3, .5, 1), lambda_component0 = 0,
#'             fixed = c(FALSE, TRUE, TRUE), final_hessian = "stable")
#' fit_diagnostics(fit)$status
#' fit_diagnostics(fit)$computation
fit_diagnostics <- function(fit) {
  metadata <- fit_metadata(fit)
  diagnostics <- metadata$diagnostics
  if (is.null(diagnostics)) diagnostics <- list(status = "not_checked", reason = "older_fit")
  if (is.null(diagnostics$computation)) diagnostics$computation <- unrecorded_fit_computation(metadata)
  diagnostics
}
