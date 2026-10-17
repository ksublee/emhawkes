# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)
old <- options(emhawkes.engine = "auto")
evaluate <- function(engine, object, dt, type = NULL, mark = NULL, initial = .1, infer = FALSE) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(logLik(object, dt, type = type, mark = mark,
    lambda_component0 = initial, infer = infer))
}
close <- function(a, b, tolerance = 1e-10) stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance)))
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), ignore.case = TRUE))
}
# Verify real dispatch, rather than just agreement after an automatic R fallback.
native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
native_calls <- 0L
assignInNamespace("hawkes_loglik_cpp", function(...) {
  native_calls <<- native_calls + 1L
  native(...)
}, "emhawkes")
poisson <- new("hspec", mu = function(param = c(log_rate = log(.7))) matrix(exp(param)),
  alpha = function(param = numeric(0)) matrix(0),
  beta = function(param = numeric(0)) matrix(1))
dt <- c(0, .2, .3, .5, 1, 1)
expected <- 5 * log(.7) - .7 * 3
for (engine in c("r", "auto", "cpp")) {
  native_calls <- 0L
  actual <- evaluate(engine, poisson, dt, initial = 0, infer = TRUE)
  close(as.numeric(actual$loglikelihood), expected)
  close(as.vector(actual$lambda), rep(.7, length(dt)))
  close(actual$res_process1, .7 * dt[-c(1, 2)])
  stopifnot(native_calls == if (engine == "r") 0L else 1L)
}
wrap <- function(value) {
  force(value)
  function(param = value) matrix(param, nrow = nrow(value))
}
set.seed(7132)
for (d in c(1L, 2L, 4L)) for (k in c(d, 2L * d)) {
  mu <- matrix(seq(.2, .4, length.out = d), d)
  alpha <- matrix(runif(d * k, .02, .07), d)
  beta <- matrix(runif(d * k, .8, 2.2), d)
  eta <- matrix(runif(d * k, -.005, .005), d)
  h <- new("hspec", mu = wrap(mu), alpha = wrap(alpha), beta = wrap(beta), eta = wrap(eta),
    type_col_map = if (k == d) NULL else lapply(seq_len(d), function(j) seq(j, k, by = d)))
  intervals <- c(0, 0, rexp(39))
  types <- c(0L, sample.int(d, 40L, replace = TRUE))
  marks <- c(0, sample(1:3, 40L, replace = TRUE))
  for (initial in list(.1, matrix(seq(.02, .08, length.out = d * k), d),
                       if (k == d) NULL else .2)) for (infer in c(FALSE, TRUE)) {
    reference <- evaluate("r", h, intervals, types, marks, initial, infer)
    for (engine in c("auto", "cpp")) {
      native_calls <- 0L
      close(evaluate(engine, h, intervals, types, marks, initial, infer), reference)
      stopifnot(native_calls == 1L)
    }
  }
  close(evaluate("cpp", h, 0, type = 0, initial = .1, infer = TRUE),
        evaluate("r", h, 0, type = 0, initial = .1, infer = TRUE))
}
# Overlapping/reordered kernel columns and event types absent from the data.
h@type_col_map <- list(c(4L, 1L), c(2L, 4L), c(3L, 1L), c(8L, 5L))
close(evaluate("cpp", h, c(0, .2, .4), c(0, 2, 2), initial = .1, infer = TRUE),
      evaluate("r", h, c(0, .2, .4), c(0, 2, 2), initial = .1, infer = TRUE))
# Reuse resolved values on fallback. Preserve default evaluation counts for
# explicit and automatic initialization, including the existing second evaluation.
calls <- c(mu = 0L, alpha = 0L, beta = 0L, eta = 0L)
component <- function(name, value) {
  force(name); force(value)
  function(param = value) { calls[name] <<- calls[name] + 1L; matrix(param) }
}
h <- new("hspec", mu = component("mu", .3), alpha = component("alpha", .1),
  beta = component("beta", 1.2), eta = component("eta", .01))
for (initial in list(.1, NULL)) {
  calls[] <- 0L
  reference <- evaluate("r", h, dt, initial = initial)
  reference_calls <- calls
  for (engine in c("auto", "cpp")) {
    calls[] <- 0L
    close(evaluate(engine, h, dt, initial = initial), reference)
    stopifnot(identical(calls, reference_calls))
  }
}
classed <- h
classed@alpha <- function(param = .1) {
  calls["alpha"] <<- calls["alpha"] + 1L
  structure(matrix(param), class = c("custom", "matrix", "array"))
}
for (initial in list(.1, NULL)) {
  calls[] <- 0L
  reference <- evaluate("r", classed, dt, initial = initial)
  reference_calls <- calls
  calls[] <- 0L; native_calls <- 0L
  close(evaluate("auto", classed, dt, initial = initial), reference)
  stopifnot(identical(calls, reference_calls), native_calls == 0L)
}
fails(evaluate("cpp", classed, dt), "plain numeric matrices")
vector_mu <- h
vector_mu@mu <- function(param = .3) param
close(evaluate("auto", vector_mu, dt), evaluate("r", vector_mu, dt))
fails(evaluate("cpp", vector_mu, dt), "plain numeric matrices")
bad <- h
bad@alpha <- function(param = .1) stop("replay callback sentinel")
fails(evaluate("auto", bad, dt), "replay callback sentinel")
fails(evaluate("cpp", bad, dt), "replay callback sentinel")
# Unsupported observation classes and mark/residual densities retain R semantics.
for (name in c("dmark", "dresidual")) {
  bad <- poisson
  slot(bad, name) <- if (name == "dmark") function(...) 1 else function(x) dexp(x)
  native_calls <- 0L
  close(evaluate("auto", bad, dt), evaluate("r", bad, dt))
  stopifnot(native_calls == 0L)
  fails(evaluate("cpp", bad, dt), "does not support")
}
# Classed observations stay on R, before resolving functions for native use.
classed_dt <- structure(dt, class = "emhawkes_replay_times")
native_calls <- 0L
close(evaluate("auto", poisson, classed_dt), evaluate("r", poisson, classed_dt))
stopifnot(native_calls == 0L)
fails(evaluate("cpp", poisson, classed_dt), "does not support")
# Function-valued models use zero components when no mean formula is supported.
singular <- new("hspec", mu = wrap(matrix(.2, 2)),
  alpha = wrap(matrix(.5, 2, 2)), beta = wrap(matrix(1, 2, 2)))
reference <- suppressWarnings(evaluate("r", singular, c(0, .2, .3), c(0, 1, 2), initial = NULL, infer = TRUE))
close(suppressWarnings(evaluate("cpp", singular, c(0, .2, .3), c(0, 1, 2), initial = NULL, infer = TRUE)), reference)
invalid_initial <- poisson
invalid_initial@alpha <- function(param = numeric(0)) matrix(1)
for (engine in c("r", "auto", "cpp")) {
  close(evaluate(engine, invalid_initial, dt, initial = NULL),
        evaluate(engine, invalid_initial, dt, initial = 0))
  fails(evaluate(engine, invalid_initial, dt, initial = "mean"), "Mean initialization is unavailable")
}
# Fitted defaults and saved initial-state arguments replay with the native engine.
fit <- suppressMessages(hfit(poisson, dt, lambda_component0 = 0, method = "NR"))
fitted <- fitted_hspec(fit)
native_calls <- 0L
replay <- suppressMessages(do.call(logLik, c(list(object = fitted, inter_arrival = dt,
  infer = TRUE), fit_initial_state(fit))))
close(as.numeric(replay$loglikelihood), as.numeric(logLik(fit)))
stopifnot(native_calls == 1L)
# Inputs and RNG state are preserved; native replay never creates a seed.
snapshot <- serialize(fitted, NULL)
set.seed(71); seed <- .Random.seed
invisible(evaluate("cpp", fitted, dt, initial = 0, infer = TRUE))
stopifnot(identical(seed, .Random.seed), identical(snapshot, serialize(fitted, NULL)))
rm(.Random.seed, envir = .GlobalEnv)
invisible(evaluate("cpp", fitted, dt, initial = 0))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
options(old)

}
