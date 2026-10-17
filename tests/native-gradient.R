# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
build_map <- getFromNamespace("build_parameter_map", "emhawkes")
restore <- getFromNamespace("restore_parameters", "emhawkes")
normalize <- getFromNamespace("normalize_observation", "emhawkes")
objective <- getFromNamespace("make_hawkes_objective", "emhawkes")
close <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tolerance)))
}
fails <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), ignore.case = TRUE))
}

# Independent O(N^2) likelihood: sum all earlier events rather than recurse.
# Event order remains significant when two supplied times coincide.
history_loglik <- function(dt, type, mark, mu, alpha, beta, eta, initial, mapping, ...) {
  time <- cumsum(dt)
  horizon <- tail(time, 1)
  events <- seq_along(time)[-1]
  integral <- function(rate, duration) -expm1(-rate * duration) / rate
  value <- -sum(mu) * horizon - sum(initial * integral(beta, horizon))
  for (n in events) {
    state <- initial * exp(-beta * time[n])
    for (k in events[events < n]) {
      columns <- mapping[[type[k]]] + 1L
      jump <- alpha[, columns, drop = FALSE]
      if (!is.null(eta)) jump <- jump + eta[, columns, drop = FALSE] * (mark[k] - 1)
      state[, columns] <- state[, columns, drop = FALSE] +
        jump * exp(-beta[, columns, drop = FALSE] * (time[n] - time[k]))
    }
    value <- value + log(mu[type[n]] + sum(state[type[n], ]))
    columns <- mapping[[type[n]]] + 1L
    jump <- alpha[, columns, drop = FALSE]
    if (!is.null(eta)) jump <- jump + eta[, columns, drop = FALSE] * (mark[n] - 1)
    value <- value - sum(jump * integral(beta[, columns, drop = FALSE], horizon - time[n]))
  }
  value
}

central <- function(fn, x, step) {
  vapply(seq_along(x), function(k) {
    plus <- minus <- x
    h <- step * max(1, abs(x[k]))
    plus[k] <- plus[k] + h
    minus[k] <- minus[k] - h
    (fn(plus) - fn(minus)) / (2 * h)
  }, numeric(1))
}

args <- list(dt = c(0, .3, 0, .7, .2, 1.1, .05),
             type = c(0L, 1L, 2L, 2L, 1L, 2L, 1L),
             mark = c(0, 1, 3, .5, 2, 1, 2.5), mu = c(.8, .9),
             alpha = matrix(c(.08, -.02, .12, .06, .08, .06, .03, .09), 2),
             beta = matrix(c(.7, 1.1, 2.3, .9, .7, 1.6, 3.1, .9), 2),
             eta = matrix(c(-.01, .02, .03, -.005, -.01, .02, .015, .01), 2),
             initial = matrix(c(.11, -.03, .07, .13, -.02, .08, .02, .05), 2),
             mapping = list(c(2L, 0L, 3L), c(3L, 1L)), infer = FALSE)
snapshot <- serialize(args, NULL)
set.seed(178)
seed <- .Random.seed
scored <- do.call(native, c(args, list(gradient = TRUE)))
stopifnot(identical(seed, .Random.seed), identical(snapshot, serialize(args, NULL)))
close(scored$loglikelihood, do.call(history_loglik, args), 1e-12)
close(scored$loglikelihood, do.call(native, args)$loglikelihood, 0)
for (component in c("mu", "alpha", "beta", "eta", "initial")) {
  reference <- function(value) {
    candidate <- args
    candidate[[component]][] <- value
    do.call(history_loglik, candidate)
  }
  for (step in c(1e-4, 2e-5)) {
    close(scored$score[[component]], central(reference, args[[component]], step), 2e-7)
  }
}
inferred_args <- args
inferred_args$infer <- TRUE
inferred <- do.call(native, c(inferred_args, list(gradient = TRUE)))
for (component in names(scored$score)) close(inferred$score[[component]], scored$score[[component]])
close(inferred$lambda, do.call(native, inferred_args)$lambda, 0)

# A nearly zero decay rate must retain the integral's beta derivative.
# For zero excitation it has the explicit beta -> 0 limit below.
small <- list(dt = c(0, .2, .5, 1), type = c(0L, 1L, 1L, 1L),
              mark = c(0, 1, 1, 1), mu = .8, alpha = matrix(0),
              beta = matrix(1e-16), eta = NULL, initial = matrix(.3),
              mapping = list(0L), infer = FALSE, gradient = TRUE)
score <- do.call(native, small)$score
arrival <- cumsum(small$dt)[-1]
close(score$beta, -sum(arrival) * .3 / 1.1 + .3 * tail(arrival, 1)^2 / 2, 1e-12)
small$beta[] <- 1e200
small$initial[] <- 0
score <- do.call(native, small)$score
stopifnot(all(is.finite(unlist(score))))
close(score$beta, 0, 0)
close(score$mu, 3 / .8 - sum(small$dt), 1e-12)

# Products may be representable even when an intermediate beta^2 or dt^2
# overflows. Relative comparisons also protect very small nonzero derivatives.
extreme <- list(dt = c(0, 1), type = c(0L, 1L), mark = c(0, 1),
                mu = .8, alpha = matrix(0), beta = matrix(1e300), eta = NULL,
                initial = matrix(1e308), mapping = list(0L),
                infer = FALSE, gradient = TRUE)
score <- do.call(native, extreme)$score
stopifnot(abs(score$beta[1] / 1e-292 - 1) < 1e-12)
extreme$dt[2] <- 1e300
extreme$beta[] <- 1e-300
extreme$initial[] <- 1e-308
score <- do.call(native, extreme)$score
expected <- 1e292 * (1 - 2 / exp(1))
stopifnot(abs(score$beta[1] / expected - 1) < 1e-12)
# Zero components have exactly zero beta sensitivity despite the large weight.
extreme$initial[] <- 0
score <- do.call(native, extreme)$score
stopifnot(all(is.finite(unlist(score))))
close(score$beta, 0, 0)
# A subnormal positive event intensity has an overflowing baseline derivative.
# Unaffected component derivatives must remain zero, not zero times infinity.
extreme$dt[2] <- 1
extreme$beta[] <- 1
extreme$mu <- 1e-320
score <- do.call(native, extreme)$score
stopifnot(is.infinite(score$mu[1]), score$mu[1] > 0,
          !any(is.nan(unlist(score))))
close(score$alpha, 0, 0)
close(score$beta, 0, 0)

# The initial row alone has zero likelihood and score.
empty <- args
empty$dt <- 0
empty$type <- 0L
empty$mark <- 0
zero <- do.call(native, c(empty, list(gradient = TRUE)))
close(zero$loglikelihood, 0, 0)
stopifnot(all(unlist(zero$score) == 0))

# Mapped gradients must add every shared matrix contribution, preserve names,
# and differentiate parameter-dependent automatic initialization.
check_mapped <- function(object, dt, type = NULL, mark = NULL, initial = NULL) {
  observation <- normalize(object, dt, type, mark)
  mapping <- if (is.null(object@type_col_map)) as.list(seq_len(object@dimens)) else object@type_col_map
  mapping <- lapply(mapping, function(x) as.integer(x) - 1L)
  for (reduced in c(TRUE, FALSE)) {
    map <- build_map(object, reduced)
    fun <- objective(object, observation, initial, map)
    reference <- function(param) {
      values <- restore(map, param)
      init <- initial
      if (is.null(init)) {
        if (object@dimens == 1L) {
          init <- matrix(values$mu[1] * values$alpha[1] /
                           (2 * (values$beta[1] - values$alpha[1])))
        } else {
          branching <- values$alpha / values$beta
          stationary <- as.vector(solve(diag(object@dimens) - branching, values$mu))
          init <- sweep(branching, 2, stationary, "*")
        }
      } else init <- matrix(init, nrow(values$beta), ncol(values$beta))
      history_loglik(observation$inter_arrival, as.integer(observation$type),
                     observation$mark, as.vector(values$mu), values$alpha,
                     values$beta, values$eta, init, mapping)
    }
    param <- map$start
    score <- fun$gradient(param)
    stopifnot(identical(names(score), names(param)))
    close(fun$loglik(param), reference(param), 1e-12)
    for (step in c(1e-4, 2e-5)) close(score, central(reference, param, step), 3e-7)
    direction <- rep(c(.3, -.5, .2), length.out = length(param))
    step <- 2e-5
    close(sum(score * direction),
          (reference(param + step * direction) - reference(param - step * direction)) / (2 * step), 3e-7)
    # Interleaved value and gradient calls must not reuse a different point.
    moved <- param + 1e-3 * direction
    close(fun$loglik(moved), reference(moved), 1e-12)
    close(fun$gradient(param), score, 0)
    close(fun$gradient(moved), central(reference, moved, 2e-5), 3e-7)
    close(fun$loglik(param), reference(param), 1e-12)
    bad <- param
    bad[startsWith(names(bad), "beta")] <- 0
    stopifnot(identical(as.numeric(fun$loglik(bad)), -Inf))
    stopifnot(all(is.na(fun$gradient(bad))))
  }
}

multi <- new("hspec", mu = args$mu, alpha = args$alpha, beta = args$beta,
             eta = args$eta, type_col_map = lapply(args$mapping, function(x) x + 1L))
check_mapped(multi, args$dt, args$type, args$mark, args$initial)
univariate <- new("hspec", mu = .4, alpha = .2, beta = 1.3, eta = .03)
check_mapped(univariate, c(0, .1, .4, 0, .7, 1.2), mark = c(0, 2, 1, 3, 1, 2))
bivariate <- new("hspec", mu = c(.4, .4), alpha = matrix(c(.12, .07, .07, .1), 2),
                 beta = matrix(c(.9, 1.4, .9, 2), 2))
check_mapped(bivariate, args$dt, args$type)

# The Poisson optimum and information are known without either Hawkes engine.
# Numeric remains the default; analytic honors fixed parameters and restoration.
previous <- options(emhawkes.engine = "cpp")
poisson <- new("hspec", mu = 1, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .5, 1, 1)
fit_one <- function(...) suppressMessages(hfit(poisson, dt, lambda_component0 = 0,
  fixed = c(FALSE, TRUE, TRUE), ...))
numeric <- fit_one(gradient = "numeric")
default <- fit_one()
for (field in c("estimate", "maximum", "gradient", "hessian", "code", "iterations")) {
  stopifnot(identical(default[[field]], numeric[[field]]))
}
analytic <- fit_one(gradient = "analytic", method = "NR")
close(coef(analytic), c(5 / 3, 0, 1), 1e-7)
close(logLik(analytic), logLik(numeric), 1e-8)
close(vcov(analytic)[1, 1], (5 / 3)^2 / 5, 1e-5)
restored <- fitted_hspec(analytic)
close(logLik(analytic), do.call(logLik, c(list(object = restored, inter_arrival = dt),
                                         fit_initial_state(analytic))), 1e-10)

# Optional final curvature and a supplied Hessian retain their established roles.
stable <- fit_one(gradient = "analytic", method = "NR", final_hessian = "stable")
for (field in c("estimate", "maximum", "gradient", "code", "iterations")) {
  stopifnot(identical(stable[[field]], analytic[[field]]))
}
close(vcov(stable)[1, 1], (5 / 3)^2 / 5, 1e-7)
hessian_calls <- 0L
user_hessian <- function(param) {
  hessian_calls <<- hessian_calls + 1L
  diag(c(-5 / param[1]^2, 0, 0))
}
with_hessian <- fit_one(gradient = "analytic", hess = user_hessian,
                       method = "NR", final_hessian = "stable")
stopifnot(hessian_calls > 0L, identical(fit_diagnostics(with_hessian)$method, "supplied"))
close(maxLik::hessian(with_hessian)[1, 1], user_hessian(coef(with_hessian))[1, 1], 0)
close(coef(with_hessian), coef(analytic), 1e-7)

# User-supplied derivatives take precedence, including with the R engine.
calls <- 0L
user_gradient <- function(param) {
  calls <<- calls + 1L
  c(5 / param[1] - sum(dt), 0, 0)
}
options(emhawkes.engine = "r")
override <- fit_one(gradient = "analytic", grad = user_gradient, method = "NR")
stopifnot(calls > 0L)
close(coef(override), coef(analytic), 1e-7)
fails(fit_one(gradient = "analytic"), "analytic|gradient")
options(emhawkes.engine = "cpp")
fails(fit_one(gradient = "analytic", method = "BHHH"), "BHHH")
fails(fit_one(gradient = "analytic", finalHessian = "BHHH"), "BHHH")
fails(fit_one(gradient = "analytic", method = "NR", bhhhHessian = TRUE), "BHHH")
fails(fit_one(gradient = "analytic", constraints = list(eqA = matrix(c(1, 0, 0), 1), eqB = -1)), "equality")

# Capability checks cannot silently treat arbitrary callbacks as constants.
callback <- poisson
callback@impact <- function(...) stop("Do not evaluate this callback")
fails(suppressMessages(hfit(callback, dt, lambda_component0 = 0,
                            gradient = "analytic")), "impact|analytic|gradient")
callback <- new("hspec", mu = function(param = c(m = 1)) matrix(param[1]), alpha = 0, beta = 1)
fails(suppressMessages(hfit(callback, dt, lambda_component0 = 0,
                            gradient = "analytic")), "function|analytic|gradient")
fails(fit_one(gradient = "analytic", mylogLik = function(object) 0), "custom|mylogLik|analytic|gradient")
residual <- set_residual(poisson, c(rate = 2), dresidual = dexp, presidual = pexp, qresidual = qexp)
fails(suppressMessages(hfit(residual, dt, lambda_component0 = 0,
                            gradient = "analytic")), "residual|analytic|gradient")
singular <- new("hspec", mu = c(.4, .5), alpha = diag(2), beta = matrix(1, 2, 2))
fails(suppressMessages(hfit(singular, c(0, .3, .2), type = c(0, 1, 2),
                            gradient = "analytic")), "lambda_component0|singular")
options(previous)

}
