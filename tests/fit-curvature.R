library(emhawkes)

quiet_fit <- function(..., final_hessian = "stable") {
  suppressMessages(hfit(..., final_hessian = final_hessian))
}
close <- function(x, y, tolerance = 1e-7) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}

# Independent Poisson information, with two fixed Hawkes coefficients.
h <- new("hspec", mu = 1, alpha = 0, beta = 1)
dt <- c(0, .5, 1.5, .25)
fixed <- c(FALSE, TRUE, TRUE)
legacy <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed,
                    final_hessian = "maxLik")
# The default is exactly the maxLik path, without stabilization or its warnings.
warning_messages <- character(0)
default <- withCallingHandlers(suppressMessages(hfit(h, dt,
  lambda_component0 = 0, fixed = fixed)), warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
for (field in c("estimate", "maximum", "gradient", "hessian", "varcovar", "code", "iterations")) {
  stopifnot(identical(default[[field]], legacy[[field]]))
}
stopifnot(length(warning_messages) == 0L, is.null(default$varcovar),
          identical(fit_diagnostics(default)$status, "not_checked"),
          is.null(fit_diagnostics(default)$checks))
close(vcov(default), vcov(legacy), tolerance = 0)

# Explicitly request stabilization for the checks below.
fit <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed)
for (field in c("estimate", "maximum", "gradient", "code", "iterations")) {
  stopifnot(identical(fit[[field]], legacy[[field]]))
}
mu <- coef(fit)[1]
close(mu, 3 / sum(dt), tolerance = 1e-5)
close(maxLik::hessian(fit)[1, 1], -3 / mu^2)
close(vcov(fit)[1, 1], mu^2 / 3)
stopifnot(all(vcov(fit)[fixed, ] == 0), all(vcov(fit)[, fixed] == 0),
          identical(class(fit), class(legacy)),
          identical(fit_diagnostics(fit)$status, "ok"),
          nrow(fit_diagnostics(fit)$checks) >= 3L)
close(summary(fit)$estimate[1, "Std. error"], mu / sqrt(3))

# Saved diagnostics, model restoration, and old metadata remain usable.
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
saved <- readRDS(path)
unlink(path)
stopifnot(identical(fit_diagnostics(saved), fit_diagnostics(fit)))
close(as.numeric(logLik(saved)), as.numeric(do.call(logLik,
  c(list(object = fitted_hspec(saved), inter_arrival = dt), fit_initial_state(saved)))))
old <- fit
attr(old, "emhawkes")$diagnostics <- NULL
stopifnot(identical(fit_diagnostics(old)$reason, "older_fit"))

# Supplied analytic curvature is retained; disabling final curvature is respected.
hessian <- function(param) diag(c(-3 / param[1]^2, 0, 0))
analytic <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed, hess = hessian)
close(maxLik::hessian(analytic), hessian(coef(analytic)))
stopifnot(identical(fit_diagnostics(analytic)$method, "supplied"))
for (flag in list(FALSE, 0)) {
  disabled <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed, finalHessian = flag)
  original <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed,
                        finalHessian = flag, final_hessian = "maxLik")
  stopifnot(is.null(maxLik::hessian(disabled)), is.null(disabled$varcovar),
            identical(fit_diagnostics(disabled)$reason, "finalHessian_false"))
  # maxLik versions return either NULL or an NA matrix when curvature is disabled.
  close(vcov(disabled), vcov(original), tolerance = 0)
}
for (flag in c("BHHH", "bhhh", "Bhhh")) {
  bhhh <- suppressWarnings(quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed,
                                     finalHessian = flag))
  stopifnot(identical(fit_diagnostics(bhhh)$reason, "BHHH_requested"),
            is.null(bhhh$varcovar))
}

# An optimizer iteration limit remains distinct from a good local Hessian.
warnings <- character(0)
early <- withCallingHandlers(quiet_fit(new("hspec", mu = 1, alpha = 0, beta = 1),
  dt, lambda_component0 = 0, fixed = fixed, iterlim = 1), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(identical(fit_diagnostics(early)$status, "ok"),
          !fit_diagnostics(early)$optimizer_success,
          "optimizer_not_converged" %in% fit_diagnostics(early)$warnings,
          any(grepl("normal convergence", warnings, fixed = TRUE)))
nr <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed, hess = hessian,
  grad = function(param) c(3 / param[1] - sum(dt), 0, 0), method = "NR")
stopifnot(maxLik::returnCode(nr) %in% c(1L, 2L, 8L), fit_diagnostics(nr)$optimizer_success)
close(coef(nr)[1], 4 / 3)
# A decay with no excitation or initial transient is genuinely unidentified.
warnings <- character(0)
singular <- withCallingHandlers(quiet_fit(h, dt, lambda_component0 = 0,
  fixed = c(FALSE, TRUE, FALSE)), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(length(warnings) == 1L,
          grepl("Wald standard errors are unavailable", warnings, fixed = TRUE),
          identical(fit_diagnostics(singular)$status, "unavailable"),
          all(is.na(vcov(singular)[c(1, 3), c(1, 3)])),
          all(vcov(singular)[2, ] == 0),
          all(is.na(summary(singular)$estimate[c(1, 3), "Std. error"])))

# A constraint on the Poisson rate changes inference at its boundary.
boundary <- suppressWarnings(quiet_fit(new("hspec", mu = 2.5, alpha = 0, beta = 1),
  dt, lambda_component0 = 0, fixed = fixed,
  constraints = list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -2),
  control = list(reltol = 1e-10)))
stopifnot(abs(coef(boundary)[1] - 2) < 1e-3,
          identical(fit_diagnostics(boundary)$reason, "constraint_near_boundary"),
          is.na(vcov(boundary)[1, 1]))
constraint_status <- getFromNamespace("fit_constraint_status", "emhawkes")
for (factor in c(1e-9, 1, 1e9)) {
  stopifnot(identical(constraint_status(
    list(ineqA = factor * matrix(c(1, 0, 0), 1), ineqB = factor * -2),
    c(2, 0, 1), !fixed, 2), "constraint_near_boundary"))
}

# Non-Hawkes objectives are not silently assigned the Hawkes calculation.
custom <- quiet_fit(h, mylogLik = function(object) -(object@mu - 1)^2 -
  object@alpha^2 - (object@beta - 1)^2)
stopifnot(identical(fit_diagnostics(custom)$reason, "custom_likelihood"))
residual <- new("hspec", mu = 1, alpha = 0, beta = 1,
  dresidual = function(x) dexp(x), presidual = function(q) pexp(q),
  qresidual = function(p) qexp(p))
residual_fit <- quiet_fit(residual, dt, lambda_component0 = 0, fixed = fixed)
stopifnot(identical(fit_diagnostics(residual_fit)$reason, "residual_likelihood"),
          is.null(residual_fit$varcovar))

# A kink can have positive finite-difference curvature at every step, yet
# inverse curvature is step dependent and must not be accepted.
curvature <- getFromNamespace("stable_fit_curvature", "emhawkes")
kink <- maxLik::maxLik(function(p) -abs(p[1] - .5) - .5 * (p[2] - 1)^2,
  start = c(kink = .5, smooth = 1), method = "BFGS")
result <- curvature(kink)
stopifnot(identical(result$status, "unavailable"),
          identical(result$reason, "step_sensitive"),
          all(result$checks$positive_definite))

# A deterministic callback checks that the default and explicit maxLik paths
# make the same calls, while stable preserves the optimizer result.
evaluations <- 0L
count_impact <- function(param = numeric(0), n) {
  evaluations <<- evaluations + 1L
  matrix(0)
}
counted <- new("hspec", mu = 1, alpha = 0, beta = 1, impact = count_impact)
run_counted <- function(mode = NULL) {
  evaluations <<- 0L
  arguments <- list(object = counted, inter_arrival = dt,
                    lambda_component0 = 0, fixed = fixed)
  if (!is.null(mode)) arguments$final_hessian <- mode
  result <- suppressMessages(do.call(hfit, arguments))
  list(fit = result, evaluations = evaluations)
}
plain <- run_counted()
legacy_counted <- run_counted("maxLik")
requested <- run_counted("stable")
stopifnot(plain$evaluations == legacy_counted$evaluations,
          identical(plain$fit$estimate, requested$fit$estimate),
          identical(plain$fit$maximum, requested$fit$maximum),
          identical(fit_diagnostics(requested$fit)$status, "ok"))

# Stable score differences preserve the numerical optimizer and match known
# Poisson curvature. Explicit R execution retains likelihood differences.
score_fit <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed)
stopifnot(fit_diagnostics(score_fit)$method == "scaled_score_Richardson",
  fit_diagnostics(score_fit)$computation$gradient$used == "numeric",
  fit_diagnostics(score_fit)$curvature_source == "analytic_score")
old_engine <- options(emhawkes.engine = "r")
likelihood_fit <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed)
options(old_engine)
stopifnot(fit_diagnostics(likelihood_fit)$method == "scaled_central_Richardson")
for (field in c("estimate", "maximum", "gradient", "code", "iterations"))
  stopifnot(identical(score_fit[[field]], likelihood_fit[[field]]))
close(vcov(score_fit), vcov(likelihood_fit), 1e-7)
close(vcov(score_fit)[1, 1], coef(score_fit)[1]^2 / 3, 1e-9)
analytic_score_fit <- quiet_fit(h, dt, lambda_component0 = 0, fixed = fixed,
  gradient = "analytic")
stopifnot(fit_diagnostics(analytic_score_fit)$method == "scaled_score_Richardson",
  fit_diagnostics(analytic_score_fit)$computation$gradient$used == "analytic")

# Missing component Jacobians fall back without changing optimizer settings.
stopifnot(fit_diagnostics(requested$fit)$method == "scaled_central_Richardson",
  grepl("jacobian", fit_diagnostics(requested$fit)$score_fallback_reason))

# A correlated quadratic checks scaling, fixed coordinates, and score names.
Q <- matrix(c(4, 1, 1, 3), 2)
quadratic <- maxLik::maxLik(function(p) {
  z <- p[c(1, 3)] - c(.5, 20)
  -sum(z * (Q %*% z)) / 2
}, start = c(a = .5, fixed = 7, b = 20), fixed = c(FALSE, TRUE, FALSE), method = "BFGS")
score_calls <- 0L
exact_score <- function(p) {
  score_calls <<- score_calls + 1L
  stopifnot(p[2] == 7)
  g <- setNames(numeric(3), names(p))
  g[c(1, 3)] <- -drop(Q %*% (p[c(1, 3)] - c(.5, 20)))
  g[c(3, 2, 1)]
}
curvature <- getFromNamespace("stable_fit_curvature", "emhawkes")
result <- curvature(quadratic, score = exact_score)
stopifnot(result$status == "ok", result$method == "scaled_score_Richardson",
  score_calls == 1L + 2L * 2L * nrow(result$checks))
close(result$hessian, -Q, 1e-11)
close(result$covariance, solve(Q), 1e-11)

# Runtime score failure falls back to the existing likelihood route. An
# unidentified direction remains unavailable even when a score is supplied.
finalize <- getFromNamespace("finalize_hawkes_fit", "emhawkes")
for (bad in list(function(p) stop("score failed"), function(p) rep(NA_real_, length(p)),
                 function(p) 1, function(p) setNames(rep(0, length(p)), rep("wrong", length(p))))) {
  finalized <- finalize(legacy, h, NULL, NULL, list(), "BFGS", "stable",
    score_factory = function() list(loglik = maxLik::objectiveFn(legacy), gradient = bad))
  stopifnot(finalized$diagnostics$status == "ok",
    finalized$diagnostics$method == "scaled_central_Richardson",
    nzchar(finalized$diagnostics$score_fallback_reason))
}
failed <- curvature(kink, score = function(p) c(-sign(p[1] - .5), -(p[2] - 1)))
stopifnot(failed$status == "unavailable", failed$reason == "step_sensitive")
