# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

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
central <- function(fn, x, step) vapply(seq_along(x), function(k) {
  up <- down <- x
  delta <- step * max(1, abs(x[k]))
  up[k] <- up[k] + delta
  down[k] <- down[k] - delta
  (fn(up) - fn(down)) / (2 * delta)
}, numeric(1))
previous <- options(emhawkes.engine = "auto")

# Nonlinear transformations, fixed zeros, shared names across functions, and
# deliberately reordered Jacobian columns exercise the full parameter map.
mu <- function(param = c(scale = -.7, offset = -.9)) matrix(exp(param), 2)
alpha <- function(param = c(a = -.8, b = -1.4)) {
  matrix(c(exp(param["a"]) / 2, 0, exp(param["b"]), exp(param["a"]) / 3), 2)
}
beta <- function(param = c(scale = -.7, decay = .4)) {
  matrix(c(exp(param["decay"]), 2 + exp(param["scale"]),
           exp(param["decay"]) + exp(param["scale"]), 2 * exp(param["decay"])), 2)
}
eta <- function(param = c(a = -.8, mark_effect = .03)) {
  matrix(c(param["mark_effect"], 0, exp(param["a"]) / 20, param["mark_effect"]), 2)
}
jac <- list(
  mu = function(param) {
    value <- diag(exp(param))
    colnames(value) <- names(param)
    value[, c("offset", "scale"), drop = FALSE]
  },
  alpha = function(param) cbind(b = c(0, 0, exp(param["b"]), 0),
    a = c(exp(param["a"]) / 2, 0, 0, exp(param["a"]) / 3)),
  beta = function(param) cbind(decay = c(1, 0, 1, 2) * exp(param["decay"]),
    scale = c(0, 1, 1, 0) * exp(param["scale"])),
  eta = function(param) cbind(a = c(0, 0, exp(param["a"]) / 20, 0),
    mark_effect = c(1, 0, 0, 1)))
h <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = eta)
dt <- c(0, .3, 0, .7, .2, 1.1, .05)
type <- c(0, 1, 2, 2, 1, 2, 1)
mark <- c(0, 1, 3, .5, 2, 1, 2.5)

test_dt <- dt
test_type <- type
test_mark <- mark
check <- function(h, jac, initial, dt = test_dt, type = test_type, mark = test_mark) {
  snapshot <- serialize(h, NULL)
  set.seed(72)
  seed <- .Random.seed
  obs <- normalize(h, dt, type, mark)
  for (reduced in c(TRUE, FALSE)) {
    map <- build_map(h, reduced)
    fun <- objective(h, obs, initial, map, jacobian = jac)
    # Keep the original functions and use the R likelihood implementation.
    reference <- function(param) {
      candidate <- h
      values <- restore(map, param)
      for (name in c("mu", "alpha", "beta", "eta")) slot(candidate, name) <- values[[name]]
      old <- options(emhawkes.engine = "r")
      on.exit(options(old))
      as.numeric(suppressMessages(logLik(candidate, dt, type, mark, lambda_component0 = initial)))
    }
    theta <- map$start
    direction <- rep(c(.2, -.3, .4), length.out = length(theta))
    for (point in list(theta, theta + .07 * direction)) {
      close(fun$loglik(point), reference(point), 1e-12)
      score <- fun$gradient(point)
      stopifnot(identical(names(score), names(theta)))
      for (step in c(1e-4, 2e-5)) close(score, central(reference, point, step), 3e-7)
      close(sum(score * direction),
        (reference(point + 2e-5 * direction) - reference(point - 2e-5 * direction)) / 4e-5, 3e-7)
      close(fun$gradient(rev(point)), score, 0)
    }
    close(fun$loglik(theta), reference(theta), 1e-12)
    invalid <- theta
    invalid[1] <- Inf
    stopifnot(identical(as.numeric(fun$loglik(invalid)), -Inf), all(is.na(fun$gradient(invalid))))
  }
  stopifnot(identical(snapshot, serialize(h, NULL)), identical(seed, .Random.seed))
}
check(h, jac, matrix(c(.1, -.02, .03, .12), 2))
check(h, jac, NULL)

unary <- new("hspec", mu = function(param = c(m = -.8)) matrix(exp(param)),
  alpha = function(param = c(a = -1.1)) matrix(exp(param)),
  beta = function(param = c(b = .4)) matrix(exp(param)))
uj <- list(mu = function(param) matrix(exp(param), dimnames = list(NULL, "m")),
  alpha = function(param) matrix(exp(param), dimnames = list(NULL, "a")),
  beta = function(param) matrix(exp(param), dimnames = list(NULL, "b")))
check(unary, uj, NULL, type = rep(1, length(dt)))
# Constant Jacobian, numeric/function mixture, and no-parameter component.
mixed <- new("hspec", mu = c(.4, .5), alpha = alpha,
  beta = function(param = numeric(0)) matrix(c(1, 2, 3, 4), 2))
check(mixed, jac["alpha"], .1)
linear <- new("hspec", mu = .5,
  alpha = function(param = c(a = .1)) matrix(c(param["a"], 0, param["a"]), 1),
  beta = matrix(c(.6, 1.5, 3), 1), type_col_map = list(c(3, 1, 2)))
check(linear, list(alpha = matrix(c(1, 0, 1), dimnames = list(NULL, "a"))), .1,
      type = rep(1, length(dt)))
# Overlapping/reordered kernel mappings with function-valued mark effects.
multi <- h
multi@type_col_map <- list(c(2, 1), 2)
check(multi, jac, .1)

# Numeric matrix names can also be shared with local function parameters.
shared_numeric <- new("hspec", mu = .4,
  alpha = function(param = c(mu1 = .4)) matrix(param["mu1"] / 2), beta = 1.3)
check(shared_numeric, list(alpha = matrix(.5, dimnames = list(NULL, "mu1"))), .1,
      type = rep(1, length(dt)))
# Unknown function means use zero, avoiding an unrelated mean approximation.
overflow_model <- new("hspec",
  mu = function(param = c(s = -1)) matrix(exp(param), 2),
  alpha = function(param = c(a = -2)) diag(rep(exp(param), 2)),
  beta = function(param = c(s = -1)) matrix(1 + exp(param), 2, 2))
overflow_jac <- list(mu = function(param) matrix(rep(exp(param), 2), dimnames = list(NULL, "s")),
  alpha = function(param) matrix(c(1, 0, 0, 1) * exp(param), dimnames = list(NULL, "a")),
  beta = function(param) matrix(rep(exp(param), 4), dimnames = list(NULL, "s")))
overflow <- objective(overflow_model, normalize(overflow_model, dt, type), NULL,
  build_map(overflow_model, TRUE), jacobian = overflow_jac)
overflow_zero <- objective(overflow_model, normalize(overflow_model, dt, type), 0,
  build_map(overflow_model, TRUE), jacobian = overflow_jac)
close(overflow$loglik(c(s = 599, a = 158)), overflow_zero$loglik(c(s = 599, a = 158)), 0)
close(overflow$gradient(c(s = 599, a = 158)), overflow_zero$gradient(c(s = 599, a = 158)), 0)
# No-event data have zero horizon and therefore zero likelihood and score.
empty <- objective(unary, normalize(unary, 0), 0, build_map(unary, TRUE), jacobian = uj)
close(empty$loglik(build_map(unary, TRUE)$start), 0, 0)
close(empty$gradient(build_map(unary, TRUE)$start), c(0, 0, 0), 0)

obs <- normalize(unary, dt)
map <- build_map(unary, TRUE)
make <- function(j = uj, model = unary) objective(model, obs, .1, build_map(model, TRUE), jacobian = j)
fails(make(NULL), "requires jacobian")
fails(make(list(alpha = uj$alpha, beta = uj$beta)), "jacobian.mu")
fails(make(list(unknown = diag(2))), "named list")
fails(make(c(uj, uj[1])), "named list")
bad <- uj
bad$mu <- matrix(1, 2, 1, dimnames = list(NULL, "m"))
fails(make(bad), "one row")
bad$mu <- matrix(Inf, dimnames = list(NULL, "m"))
fails(make(bad), "finite")
bad$mu <- matrix(1, dimnames = list(NULL, "wrong"))
fails(make(bad), "column names")
bad$mu <- structure(matrix(1, dimnames = list(NULL, "m")), class = c("custom", "matrix", "array"))
fails(make(bad), "numeric matrix")
bad$mu <- function(param, ...) matrix(1)
fails(make(bad), "only.*param")
bad_jac <- jac
bad_jac$alpha <- matrix(0, 4, 2, dimnames = list(NULL, c("a", "a")))
fails(objective(h, normalize(h, dt, type, mark), .1, build_map(h, TRUE),
                jacobian = bad_jac), "column names")
bad <- unary
bad@mu <- function(param = c(m = -.8), n = 1) matrix(exp(param))
fails(make(model = bad)$gradient(build_map(bad, TRUE)$start), "incorrect dimensions")
bad <- unary
bad@mu <- function(param = c(m = -.8)) if (param[1] > 0) matrix(1, 2) else matrix(exp(param))
fun <- make(model = bad)
moved <- map$start
moved["m"] <- .1
fails(fun$loglik(moved), "constant dimensions")
bad <- uj
bad$mu <- function(param) if (param[1] > 0) matrix(1, 2) else uj$mu(param)
fun <- make(bad)
fails(fun$gradient(moved), "one row")
fun <- make()
moved["m"] <- 1000
stopifnot(identical(as.numeric(fun$loglik(moved)), -Inf), all(is.na(fun$gradient(moved))))

# The same-point cache must avoid repeated component/Jacobian evaluations.
counts <- c(value = 0L, jacobian = 0L)
counted <- unary
counted@mu <- function(param = c(m = -.8)) {
  counts["value"] <<- counts["value"] + 1L
  matrix(exp(param))
}
cj <- uj
cj$mu <- function(param) {
  counts["jacobian"] <<- counts["jacobian"] + 1L
  uj$mu(param)
}
fun <- make(cj, counted)
invisible(fun$loglik(map$start))
seen <- counts
invisible(fun$gradient(map$start))
invisible(fun$loglik(map$start))
stopifnot(identical(counts, seen))
invisible(fun$gradient(map$start + .01))
stopifnot(all(counts > seen))

# A transformed Poisson model has an independently known optimum and variance.
poisson <- new("hspec", mu = function(param = c(log_rate = 0)) matrix(exp(param)),
  alpha = function(param = numeric(0)) matrix(0),
  beta = function(param = numeric(0)) matrix(1))
pj <- list(mu = function(param) matrix(exp(param), dimnames = list(NULL, "log_rate")))
fit_one <- function(...) suppressMessages(hfit(poisson, c(0, .2, .3, .5, 1, 1),
  lambda_component0 = 0, method = "NR", ...))
numeric_fit <- fit_one()
analytic_fit <- fit_one(gradient = "analytic", jacobian = pj)
close(coef(analytic_fit), log(5 / 3), 1e-7)
close(logLik(analytic_fit), logLik(numeric_fit), 1e-9)
close(vcov(analytic_fit), 1 / 5, 1e-6)
stable <- fit_one(gradient = "analytic", jacobian = pj, final_hessian = "stable")
close(vcov(stable), 1 / 5, 1e-7)
stopifnot(identical(coef(stable), coef(analytic_fit)))
saved <- tempfile(fileext = ".rds")
saveRDS(analytic_fit, saved)
reloaded <- readRDS(saved)
unlink(saved)
stopifnot(identical(attr(reloaded, "emhawkes")$jacobian, pj))
restored <- fitted_hspec(reloaded)
stopifnot(is.function(restored@mu), is.function(restored@alpha))
close(logLik(reloaded), suppressMessages(do.call(logLik,
  c(list(object = restored, inter_arrival = c(0, .2, .3, .5, 1, 1)), fit_initial_state(reloaded)))), 1e-12)
# Fixed optimization parameters keep their values and zero covariance.
fixed_model <- new("hspec",
  mu = function(param = c(log_rate1 = 0, log_rate2 = 0)) matrix(exp(param), 2),
  alpha = function(param = numeric(0)) matrix(0, 2, 2),
  beta = function(param = numeric(0)) matrix(1, 2, 2))
fixed_jac <- list(mu = function(param) {
  out <- diag(exp(param)); colnames(out) <- names(param); out
})
fixed_fit <- suppressMessages(hfit(fixed_model, c(0, .2, .3, .5, 1, 1),
  type = c(0, 1, 2, 1, 2, 1), lambda_component0 = 0, fixed = c(FALSE, TRUE),
  method = "NR", gradient = "analytic", jacobian = fixed_jac))
close(coef(fixed_fit), c(0, 0), 1e-9)
close(diag(vcov(fixed_fit)), c(1 / 3, 0), 1e-6)

# A user grad bypasses the adapter, so it does not require a component Jacobian.
options(emhawkes.engine = "r")
overridden <- fit_one(gradient = "analytic", grad = function(param) 5 - 3 * exp(param))
close(coef(overridden), coef(analytic_fit), 1e-8)
options(previous)

}
