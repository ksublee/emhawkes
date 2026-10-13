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

stable_fit_curvature <- function(fit) {
  parameters <- stats::coef(fit)
  active <- maxLik::activePar(fit)
  free <- parameters[active]
  scale <- pmax(abs(free), .1)
  objective <- maxLik::objectiveFn(fit)
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
  for (index in seq_along(steps)) {
    derivative <- fit_scaled_derivatives(fn, free, scale, steps[index], value)
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
  common <- list(scale = scale, checks = do.call(rbind, checks),
                 change = change, original_hessian = maxLik::hessian(fit))
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
  c(list(status = "ok", reason = NULL, method = "scaled_central_Richardson",
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
                               final_hessian) {
  unchecked <- function(reason) {
    list(fit = fit, diagnostics = list(status = "not_checked", reason = reason))
  }
  if (final_hessian == "maxLik") return(unchecked("maxLik_requested"))
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
  active <- maxLik::activePar(fit)
  if (!is.numeric(parameters) || any(!is.finite(parameters)) || !any(active) ||
      !is.finite(as.numeric(stats::logLik(fit)))) return(unchecked("fit_unavailable"))
  scale <- pmax(abs(parameters[active]), .1)
  boundary <- fit_constraint_status(constraints, parameters, active, scale)
  if (!is.null(boundary)) {
    diagnostics <- list(status = "unavailable", reason = boundary, scale = scale)
  } else if (!is.null(hess)) {
    H <- maxLik::hessian(fit)[active, active, drop = FALSE]
    information <- fit_information(H * outer(scale, scale))
    diagnostics <- list(status = if (information$positive) "ok" else "unavailable",
      reason = if (information$positive) NULL else information$reason,
      method = "supplied", scale = scale, hessian = H,
      covariance = if (information$positive) information$covariance * outer(scale, scale) else NULL,
      scaled_condition = information$condition)
  } else {
    diagnostics <- stable_fit_curvature(fit)
  }
  if (identical(diagnostics$status, "ok") && identical(diagnostics$method, "supplied")) {
    gradient <- maxLik::gradient(fit)[active]
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
  diagnostics$optimizer_code <- maxLik::returnCode(fit)
  diagnostics$optimizer_message <- maxLik::returnMessage(fit)
  numerical_methods <- c("nr", "newton-raphson", "bfgsr", "bfgs-r", "bhhh")
  diagnostics$optimizer_success <- if (tolower(method) == "sann") NA else if (
    tolower(method) %in% numerical_methods) diagnostics$optimizer_code %in% c(1L, 2L, 8L) else
      diagnostics$optimizer_code == 0L
  diagnostics$free_parameters <- names(parameters)[active]
  if (identical(diagnostics$status, "ok")) {
    covariance[active, active] <- diagnostics$covariance
    if (identical(diagnostics$method, "scaled_central_Richardson")) {
      H <- matrix(0, length(parameters), length(parameters), dimnames = dimnames(covariance))
      H[active, active] <- diagnostics$hessian
      fit$hessian <- H
    }
    diagnostics$warnings <- character(0)
    if (!is.na(diagnostics$optimizer_success) && !diagnostics$optimizer_success) {
      diagnostics$warnings <- c(diagnostics$warnings, "optimizer_not_converged")
      warning("The optimizer did not report normal convergence (code ",
              diagnostics$optimizer_code, "). Stable curvature does not establish a reliable fit; ",
              "see fit_diagnostics().", call. = FALSE)
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
    warning("Wald standard errors are unavailable (", diagnostics$reason,
            "). See fit_diagnostics(); consider profile likelihood or bootstrap.", call. = FALSE)
  }
  fit$varcovar <- covariance
  list(fit = fit, diagnostics = diagnostics)
}

#' Inspect Final Curvature Diagnostics
#'
#' Return the final Hessian checks requested through [hfit()].
#' @param fit A result from [hfit()].
#' @return A list describing the calculation status, reason, free parameters,
#' scale, and, when available, difference-step checks, curvature changes,
#' recalculated score, scaled condition number, and Newton decrement.
#' @details
#' `status = "ok"` means the numerical curvature checks passed, not that the
#' optimizer found a global maximum or that Wald intervals have correct coverage.
#' The numerical calculation checks information and standard-error stability
#' across steps. The stored optimizer gradient and coefficients are not changed.
#' A Newton decrement above 0.1 triggers a score warning. Its square divided by
#' two estimates the local quadratic improvement in log-likelihood.
#' An information condition number above `1e8` triggers a weak-curvature warning;
#' eigenvalues at or below `1e-10` times the largest absolute eigenvalue do not
#' support numerical inversion. These thresholds are numerical diagnostics.
#' With the default `final_hessian = "maxLik"`, these additional checks are not
#' performed and the status is `"not_checked"`. Request `final_hessian = "auto"`
#' to run them. Earlier fits without these checks also have status `"not_checked"`.
#' @seealso [hfit()], [fitted_hspec()]
#' @export
#' @examples
#' h <- new("hspec", mu = .3, alpha = 0, beta = 1)
#' fit <- hfit(h, c(0, .2, .3, .5, 1), lambda_component0 = 0,
#'             fixed = c(FALSE, TRUE, TRUE), final_hessian = "auto")
#' fit_diagnostics(fit)$status
fit_diagnostics <- function(fit) {
  diagnostics <- fit_metadata(fit)$diagnostics
  if (is.null(diagnostics)) list(status = "not_checked", reason = "older_fit") else diagnostics
}
