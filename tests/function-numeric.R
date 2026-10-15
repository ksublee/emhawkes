library(emhawkes)

fit_engine <- function(engine, ...) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(hfit(...))
}
close <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tolerance)))
}
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}
map_for <- getFromNamespace("build_parameter_map", "emhawkes")
restore <- getFromNamespace("restore_parameters", "emhawkes")
normalize <- getFromNamespace("normalize_observation", "emhawkes")
objective <- getFromNamespace("make_hawkes_objective", "emhawkes")

# The Poisson special case independently fixes the optimum and covariance.
poisson <- new("hspec", mu = function(param = c(log_rate = 0)) matrix(exp(param)),
  alpha = function(param = numeric(0)) matrix(0),
  beta = function(param = numeric(0)) matrix(1))
dt <- c(0, .2, .3, .5, 1, 1)
fit_poisson <- function(engine, ...) fit_engine(engine, poisson, dt,
  lambda_component0 = 0, method = "NR", ...)
check_route <- function() {
  native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
  calls <- 0L
  assignInNamespace("hawkes_loglik_cpp", function(...) {
    calls <<- calls + 1L
    native(...)
  }, "emhawkes")
  on.exit(assignInNamespace("hawkes_loglik_cpp", native, "emhawkes"))
  auto <- fit_poisson("auto")
  stopifnot(calls == 0L)
  compiled <- fit_poisson("cpp", jacobian = stop("numeric fitting must not force jacobian"))
  stopifnot(calls > 0L, attr(compiled, "emhawkes")$gradient == "numeric",
            is.null(attr(compiled, "emhawkes")$jacobian))
  close(coef(compiled), log(5 / 3), 1e-7)
  close(vcov(compiled), 1 / 5, 1e-6)
  close(coef(auto), coef(compiled), 1e-7)
  close(logLik(auto), logLik(compiled), 1e-10)
  compiled
}
snapshot <- serialize(poisson, NULL)
set.seed(324)
seed <- .Random.seed
fit <- check_route()
stopifnot(identical(snapshot, serialize(poisson, NULL)), identical(seed, .Random.seed))
restored <- fitted_hspec(fit)
stopifnot(is.function(restored@mu), is.function(restored@alpha))
previous <- options(emhawkes.engine = "r")
close(suppressMessages(do.call(logLik, c(list(object = restored, inter_arrival = dt),
  fit_initial_state(fit)))), logLik(fit), 1e-10)
options(previous)
stable <- fit_poisson("cpp", final_hessian = "stable")
close(coef(stable), coef(fit), 0)
close(vcov(stable), 1 / 5, 1e-6)

# User derivatives retain precedence; no component Jacobians are needed.
calls <- c(gradient = 0L, hessian = 0L)
grad <- function(param) {
  calls["gradient"] <<- calls["gradient"] + 1L
  5 - 3 * exp(param)
}
hess <- function(param) {
  calls["hessian"] <<- calls["hessian"] + 1L
  matrix(-3 * exp(param))
}
overridden <- fit_poisson("cpp", gradient = "analytic", grad = grad, hess = hess)
close(coef(overridden), log(5 / 3), 1e-7)
stopifnot(all(calls > 0L), attr(overridden, "emhawkes")$gradient == "user")
fixed <- new("hspec", mu = function(param = c(rate1 = 0, rate2 = 0)) matrix(exp(param), 2),
  alpha = function(param = numeric(0)) matrix(0, 2, 2),
  beta = function(param = numeric(0)) matrix(1, 2, 2))
fixed_fit <- fit_engine("cpp", fixed, dt, type = c(0, 1, 2, 1, 2, 1),
  lambda_component0 = 0, fixed = c(FALSE, TRUE), method = "NR")
close(coef(fixed_fit), c(0, 0), 1e-8)
fixed_reference <- fit_engine("r", fixed, dt, type = c(0, 1, 2, 1, 2, 1),
  lambda_component0 = 0, fixed = c(FALSE, TRUE), method = "NR")
close(vcov(fixed_fit), vcov(fixed_reference), 1e-10)
# maxLik's default finite-difference Hessian has cancellation near log_rate = 0.
close(diag(vcov(fixed_fit)), c(1 / 3, 0), 1e-3)
automatic <- fit_engine("cpp", fixed, dt, type = c(0, 1, 2, 1, 2, 1), method = "NR")
close(coef(automatic), c(0, log(2 / 3)), 1e-7)

# Nonlinear components, shared names, fixed zeros, mark effects and unequal rates.
mu <- function(param = c(scale = log(.3))) matrix(c(1, 1.25) * exp(param), 2)
alpha <- function(param = c(scale = log(.3), a = log(.08))) {
  matrix(c(.2 * exp(param["scale"]), 0, exp(param["a"]), .2 * exp(param["scale"])), 2)
}
beta <- function(param = c(scale = log(.3), decay = log(1.2))) {
  matrix(c(1, 1.5, 2, 3) * exp(param["decay"]) + exp(param["scale"]), 2)
}
eta <- function(param = c(mark_gain = .02)) matrix(c(param, 0, param, param), 2)
h <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = eta)
intervals <- c(0, .3, 0, .7, .2, 1.1, .05)
types <- c(0, 1, 2, 2, 1, 2, 1)
marks <- c(0, 1, 3, .5, 2, 1, 2.5)
reference <- function(model, map, point, initial) {
  candidate <- model
  values <- restore(map, point)
  for (name in c("mu", "alpha", "beta", "eta")) slot(candidate, name) <- values[[name]]
  previous <- options(emhawkes.engine = "r")
  on.exit(options(previous))
  as.numeric(suppressWarnings(suppressMessages(logLik(candidate, intervals,
    type = types, mark = marks, lambda_component0 = initial))))
}
check <- function(model, initial) {
  for (reduced in c(TRUE, FALSE)) {
    map <- map_for(model, reduced)
    fun <- objective(model, normalize(model, intervals, types, marks), initial,
      map, derivatives = FALSE)
    stopifnot(is.null(fun$gradient))
    for (point in list(map$start, map$start + rep(c(.01, -.02), length.out = length(map$start)))) {
      close(fun$loglik(point), reference(model, map, point, initial), 1e-12)
      close(fun$loglik(rev(point)), reference(model, map, point, initial), 1e-12)
    }
    invalid <- map$start
    invalid[1] <- Inf
    stopifnot(identical(as.numeric(fun$loglik(invalid)), -Inf))
  }
}
check(h, matrix(c(.1, -.02, .03, .12), 2))
check(h, NULL)
mapped <- h
mapped@type_col_map <- list(c(2, 1), 2)
check(mapped, .1)

# Singular multivariate initialization keeps the existing zero-state fallback.
singular <- new("hspec", mu = c(.4, .5),
  alpha = function(param = c(a = .5)) diag(rep(param, 2)),
  beta = function(param = numeric(0)) matrix(1, 2, 2))
singular_map <- map_for(singular, TRUE)
singular_fun <- objective(singular, normalize(singular, intervals, types, marks),
  NULL, singular_map, derivatives = FALSE)
point <- singular_map$start
point["a"] <- 1
close(suppressWarnings(singular_fun$loglik(point)), reference(singular, singular_map, point, NULL), 1e-12)

# A singular univariate automatic state is an invalid numeric proposal.
unary <- new("hspec", mu = function(param = c(m = .4)) matrix(param),
  alpha = function(param = c(a = .5)) matrix(param), beta = 1)
unary_map <- map_for(unary, TRUE)
unary_point <- unary_map$start
unary_point["a"] <- 1
unary_fun <- objective(unary, normalize(unary, dt), NULL, unary_map, derivatives = FALSE)
stopifnot(identical(as.numeric(unary_fun$loglik(unary_point)), -Inf))
unary_explicit <- objective(unary, normalize(unary, dt), .1, unary_map, derivatives = FALSE)
stopifnot(is.finite(unary_explicit$loglik(unary_point)))

# Unsupported state callbacks, latent impacts, distributions, classes and shape changes
# are rejected before they can silently enter the native loop.
bad <- poisson
bad@mu <- function(param = c(log_rate = 0), lambda) matrix(exp(param))
fails(fit_engine("cpp", bad, dt, lambda_component0 = 0), "latent intensity")
bad <- poisson
bad@impact <- function(n, ...) matrix(..6)
fails(fit_engine("cpp", bad, dt, lambda_component0 = 0), "latent intensity")
bad <- poisson
bad@dmark <- function(...) 1
fails(fit_engine("cpp", bad, dt, lambda_component0 = 0), "mark density")
bad <- poisson
bad@mu <- function(param = c(log_rate = 0)) {
  structure(matrix(exp(param)), class = c("custom", "matrix", "array"))
}
fails(fit_engine("cpp", bad, dt, lambda_component0 = 0), "plain numeric matrix")
bad <- poisson
bad@mu <- function(param = c(log_rate = -.5)) {
  if (param > 0) matrix(1, 2) else matrix(exp(param))
}
map <- map_for(bad, TRUE)
fun <- objective(bad, normalize(bad, dt), 0, map, derivatives = FALSE)
fails(fun$loglik(c(log_rate = .5)), "constant dimensions")
# A custom likelihood continues to bypass native model selection.
custom <- fit_engine("cpp", poisson, mylogLik = function(object, param) -(param[1] - .2)^2,
  method = "NR")
stopifnot(abs(coef(custom) - .2) < 1e-6, attr(custom, "emhawkes")$custom_loglik)
