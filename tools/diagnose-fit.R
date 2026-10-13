# Development diagnostics; no fit coefficients or Hessians are replaced.
# Rscript tools/diagnose-fit.R <fit.rds> <output-directory>

central_fit_derivatives <- function(fn, parameters, scale, step) {
  count <- length(parameters)
  value <- as.numeric(fn(parameters))
  gradient <- setNames(numeric(count), names(parameters))
  hessian <- matrix(0, count, count,
                    dimnames = list(names(parameters), names(parameters)))
  for (i in seq_len(count)) {
    delta_i <- numeric(count)
    delta_i[i] <- scale[i] * step
    above <- as.numeric(fn(parameters + delta_i))
    below <- as.numeric(fn(parameters - delta_i))
    gradient[i] <- (above - below) / (2 * delta_i[i])
    hessian[i, i] <- ((above - value) + (below - value)) / delta_i[i]^2
    if (i < count) for (j in (i + 1L):count) {
      delta_j <- numeric(count)
      delta_j[j] <- scale[j] * step
      hessian[i, j] <- hessian[j, i] <-
        (as.numeric(fn(parameters + delta_i + delta_j)) -
         as.numeric(fn(parameters + delta_i - delta_j)) -
         as.numeric(fn(parameters - delta_i + delta_j)) +
         as.numeric(fn(parameters - delta_i - delta_j))) /
        (4 * delta_i[i] * delta_j[j])
    }
  }
  list(gradient = gradient, hessian = hessian)
}

diagnose_fit <- function(fit, steps = c(.003, .001, .0003, .0001), scale = NULL) {
  # Loading maxLik registers its S3 methods for serialized fits.
  objective <- maxLik::objectiveFn(fit)
  parameters <- stats::coef(fit)
  active <- maxLik::activePar(fit)
  if (!inherits(fit, "maxLik") || !is.function(objective) ||
      !is.numeric(parameters) || any(!is.finite(parameters)) || !any(active)) {
    stop("A maxLik fit with a stored objective and finite free coefficients is required.")
  }
  free <- parameters[active]
  if (is.null(scale)) scale <- pmax(abs(free), .1)
  if (length(scale) != length(free) || any(!is.finite(scale)) || any(scale <= 0) ||
      !length(steps) || any(!is.finite(steps)) || any(steps <= 0)) {
    stop("Supply positive finite scales and relative difference steps.")
  }
  fn <- function(candidate) {
    proposal <- parameters
    proposal[active] <- candidate
    sum(as.numeric(objective(proposal)))
  }
  value <- fn(free)
  if (!is.finite(value) || !isTRUE(all.equal(value, as.numeric(stats::logLik(fit)),
                                           tolerance = 1e-8))) {
    stop("The stored objective does not reproduce the fitted likelihood.")
  }
  stored_hessian <- maxLik::hessian(fit)
  if (is.null(stored_hessian)) {
    stored_hessian <- matrix(NA_real_, length(free), length(free),
                             dimnames = list(names(free), names(free)))
  } else {
    stored_hessian <- stored_hessian[active, active, drop = FALSE]
  }
  stored <- list(gradient = maxLik::gradient(fit)[active],
                 hessian = stored_hessian)
  unavailable_standard_errors <- setNames(rep(NA_real_, length(parameters)),
                                         names(parameters))
  stored_standard_errors <- tryCatch({
    covariance <- stats::vcov(fit)
    if (is.null(covariance)) unavailable_standard_errors else
      suppressWarnings(sqrt(diag(covariance)))
  }, error = function(error) unavailable_standard_errors)
  derivatives <- c(list(stored = stored), setNames(lapply(steps, function(step) {
    central_fit_derivatives(fn, free, scale, step)
  }), paste0("central_", format(steps, scientific = TRUE))))
  rows <- list()
  standard_errors <- list()
  for (method in names(derivatives)) {
    derivative <- derivatives[[method]]
    information <- -(derivative$hessian + t(derivative$hessian)) / 2
    finite <- all(is.finite(information))
    eigenvalues <- if (finite) eigen(information, symmetric = TRUE)$values else NA_real_
    scaled <- information * outer(scale, scale)
    scaled_eigenvalues <- if (finite) eigen(scaled, symmetric = TRUE)$values else NA_real_
    positive <- finite && min(scaled_eigenvalues) > 0
    # A finite SE selected from an indefinite matrix is not endorsed here.
    se <- if (positive) scale * sqrt(diag(chol2inv(chol(scaled)))) else rep(NA_real_, length(free))
    rows[[method]] <- data.frame(
      method = method, minimum_information_eigenvalue = min(eigenvalues),
      maximum_information_eigenvalue = max(eigenvalues),
      minimum_scaled_eigenvalue = min(scaled_eigenvalues),
      scaled_condition = if (positive) max(scaled_eigenvalues) / min(scaled_eigenvalues) else NA_real_,
      positive_definite = positive,
      maximum_scaled_hessian_asymmetry = max(abs((derivative$hessian - t(derivative$hessian)) * outer(scale, scale))),
      maximum_absolute_score = max(abs(derivative$gradient)),
      maximum_scaled_score = max(abs(derivative$gradient * scale)))
    standard_errors[[method]] <- data.frame(method = method, parameter = names(free),
                                           estimate = unname(free), standard_error = se)
  }
  list(code = maxLik::returnCode(fit), loglik = value, scale = scale,
       stored_standard_errors = stored_standard_errors,
       curvature = do.call(rbind, rows), standard_errors = do.call(rbind, standard_errors),
       derivatives = derivatives)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: Rscript tools/diagnose-fit.R <fit.rds> <output-directory>")
  suppressPackageStartupMessages(library(emhawkes))
  fit <- readRDS(args[1L])
  result <- diagnose_fit(fit)
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  write.csv(result$curvature, file.path(args[2L], "curvature.csv"), row.names = FALSE)
  write.csv(result$standard_errors, file.path(args[2L], "standard-errors.csv"), row.names = FALSE)
  saveRDS(result, file.path(args[2L], "diagnostics.rds"))
  capture.output(sessionInfo(), file = file.path(args[2L], "session-info.txt"))
  cat("Return code:", result$code, "\nLog-likelihood:", result$loglik, "\n")
  print(result$curvature, row.names = FALSE)
  cat("Check multiple steps and starting values. Constraints and nonsmooth objectives require separate analysis.\n")
}
