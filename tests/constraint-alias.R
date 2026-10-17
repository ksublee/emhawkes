# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)
capture_fit <- function(...) {
  warnings <- character()
  value <- withCallingHandlers(suppressMessages(hfit(...)), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(fit = value, warnings = warnings)
}
same_fit <- function(first, second) {
  for (field in c("estimate", "maximum", "gradient", "hessian", "code", "iterations", "varcovar"))
    stopifnot(identical(first$fit[[field]], second$fit[[field]]))
  stopifnot(identical(first$warnings, second$warnings),
    identical(fit_diagnostics(first$fit), fit_diagnostics(second$fit)),
    identical(attr(first$fit, "emhawkes")$constraints, attr(second$fit, "emhawkes")$constraints),
    identical(fit_initial_state(first$fit), fit_initial_state(second$fit)))
}
h <- new("hspec", mu = .5, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .5)
fixed <- c(FALSE, TRUE, TRUE)
restriction <- list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -.1)
old <- options(emhawkes.engine = "auto")
for (engine in c("r", "cpp")) for (curvature in c("maxLik", "stable")) {
  options(emhawkes.engine = engine)
  arguments <- list(object = h, inter_arrival = dt, lambda_component0 = 0,
    fixed = fixed, final_hessian = curvature)
  plural <- do.call(capture_fit, c(arguments, list(constraints = restriction)))
  singular <- do.call(capture_fit, c(arguments, list(constraint = restriction)))
  same_fit(plural, singular)
  stopifnot(!any(grepl("deprecated", singular$warnings, fixed = TRUE)),
    identical(attr(singular$fit, "emhawkes")$constraints, restriction))
  # Explicit NULL for the unused spelling does not conflict.
  same_fit(plural, do.call(capture_fit, c(arguments,
    list(constraint = restriction, constraints = NULL))))
  same_fit(plural, do.call(capture_fit, c(arguments,
    list(constraints = restriction, constraint = NULL))))
}
# Existing positional constraints/method/verbose arguments retain their positions.
positional <- do.call(capture_fit, c(list(h, dt, NULL, NULL, NULL, NULL, 0, NULL,
  NULL, TRUE, NULL, NULL, restriction, "BFGS", FALSE), list(fixed = fixed)))
named <- capture_fit(h, dt, lambda_component0 = 0, fixed = fixed, constraints = restriction)
same_fit(positional, named)
same_fit(capture_fit(h, dt, lambda_component0 = 0, fixed = fixed),
  capture_fit(h, dt, lambda_component0 = 0, fixed = fixed, constraint = NULL))
restored <- fitted_hspec(named$fit)
stopifnot(isTRUE(all.equal(as.numeric(logLik(named$fit)),
  as.numeric(logLik(restored, dt, lambda_component0 = 0)), tolerance = 1e-12)))

# Ambiguous simultaneous values fail before a custom objective can execute.
calls <- 0L
error <- tryCatch(hfit(h, constraint = restriction, constraints = restriction,
  mylogLik = function(object) { calls <<- calls + 1L; 0 }), error = identity)
stopifnot(inherits(error, "error"), calls == 0L, is.null(conditionCall(error)),
  grepl("only one of 'constraint' or 'constraints'", conditionMessage(error), fixed = TRUE))
# The alias also reaches the existing analytic equality-constraint guard.
equality <- list(eqA = matrix(c(1, 0, 0), 1), eqB = -2)
errors <- lapply(c("constraint", "constraints"), function(name) tryCatch(
  do.call(hfit, c(list(object = h, inter_arrival = dt, lambda_component0 = 0,
    gradient = "analytic"), setNames(list(equality), name))), error = identity))
stopifnot(all(vapply(errors, inherits, logical(1), "error")),
  identical(conditionMessage(errors[[1]]), conditionMessage(errors[[2]])),
  grepl("equality constraints", conditionMessage(errors[[1]]), fixed = TRUE))
options(old)

}
