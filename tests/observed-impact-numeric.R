library(emhawkes)
previous <- options(emhawkes.engine = "r")
old_jit <- compiler::enableJIT(0)
build_map <- emhawkes:::build_parameter_map
normalize <- emhawkes:::normalize_observation
objective <- emhawkes:::make_hawkes_objective
close <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tolerance)))
}
fails <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), ignore.case = TRUE))
}
fit_engine <- function(engine, ...) {
  saved <- options(emhawkes.engine = engine)
  on.exit(options(saved))
  suppressWarnings(suppressMessages(hfit(...)))
}
# Numeric core slots must also dispatch to the adapter when an impact is present.
impact <- function(param = c(q = .06), n, N, Nc, mark, type, inter_arrival,
                   mu, alpha, beta, offset = .005, ...) {
  matrix(param["q"] * (1 + N[n, 1] %% 2) + .003 * Nc[n - 1, 1] +
    .002 * mark[n] + .001 * type[n] + .001 * inter_arrival[n] +
    offset * (mu[1] + alpha[1] + beta[1]))
}
h <- new("hspec", mu = .5, alpha = .08, beta = 1.4, impact = impact)
dt <- c(0, .2, .7, .1, 0, .6, .2, 1.2, .5, .3)
marks <- c(0, 1, 2, 1, 3, .5, 2, 1, 1, 2)
arguments <- list(object = h, inter_arrival = dt, mark = marks, N0 = 2, Nc0 = 3,
  lambda_component0 = .1, method = "BFGS", control = list(iterlim = 6))
check_route <- function() {
  native <- emhawkes:::hawkes_loglik_cpp
  calls <- 0L
  assignInNamespace("hawkes_loglik_cpp", function(...) {
    args <- list(...)
    stopifnot(!is.null(args$impact), identical(args$gradient, FALSE))
    calls <<- calls + 1L
    native(...)
  }, "emhawkes")
  on.exit(assignInNamespace("hawkes_loglik_cpp", native, "emhawkes"))
  auto <- do.call(fit_engine, c(list(engine = "auto"), arguments))
  stopifnot(calls == 0L)
  cpp <- do.call(fit_engine, c(list(engine = "cpp", jacobian = quote(stop("unused Jacobian"))), arguments))
  stopifnot(calls > 0L, attr(cpp, "emhawkes")$gradient == "numeric",
            is.null(attr(cpp, "emhawkes")$jacobian))
  close(coef(cpp), coef(auto), 1e-8)
  close(logLik(cpp), logLik(auto), 1e-11)
  close(cpp$gradient, auto$gradient, 1e-6)
  close(cpp$hessian, auto$hessian, 2e-3)
  stopifnot(identical(cpp$iterations, auto$iterations), identical(cpp$code, auto$code))
  cpp
}
snapshot <- serialize(h, NULL)
set.seed(625)
seed <- .Random.seed
fit <- check_route()
stopifnot(identical(snapshot, serialize(h, NULL)), identical(seed, .Random.seed))
restored <- fitted_hspec(fit)
close(suppressMessages(do.call(logLik, c(list(object = restored, inter_arrival = dt,
  mark = marks), fit_initial_state(fit)))), logLik(fit), 1e-11)

# Independent full-history calculation: signed jumps, unequal rates, zero waiting
# times and overlapping/reordered columns, with nonzero count initial states.
types <- c(0, 1, 2, 1, 2, 2, 1, 2, 1, 2)
model <- new("hspec", mu = c(.8, .9), alpha = matrix(c(.08, -.02, .04, .07, .01, .03), 2),
  beta = matrix(c(.7, 1.1, 2, .9, 3, 1.5), 2), eta = matrix(.01, 2, 3),
  type_col_map = list(c(3, 1), c(3, 2)),
  impact = function(param = c(q = .02), n, N, Nc, mark, ...) {
    matrix(param["q"] * (N[n, 1] - N[n - 1, 2]) + .001 * Nc[n, 2] +
      .002 * mark[n] + c(-.01, .02, 0, .01, -.02, .01), 2)
  })
initial <- matrix(c(.1, -.01, .03, .06, .02, .01), 2)
obs <- normalize(model, dt, types, marks, N0 = c(2, 3), Nc0 = c(4, 5))
full_history <- function(values, initial) {
  mu <- as.vector(values$mu)
  beta <- values$beta
  times <- cumsum(dt)
  horizon <- tail(times, 1)
  integral <- function(t) -expm1(-beta * t) / beta
  result <- -sum(mu) * horizon - sum(initial * integral(horizon))
  jumps <- lapply(2:length(dt), function(n) {
    extra <- values$impact(n = n, N = obs$N, Nc = obs$Nc, mark = marks)
    jump <- values$alpha + values$eta * (marks[n] - 1) + extra
    unused <- setdiff(seq_len(ncol(beta)), model@type_col_map[[types[n]]])
    jump[, unused] <- 0
    jump
  })
  for (n in 2:length(dt)) {
    state <- initial * exp(-beta * times[n])
    if (n > 2) for (k in 2:(n - 1)) state <- state + jumps[[k - 1]] * exp(-beta * (times[n] - times[k]))
    result <- result + log(mu[types[n]] + sum(state[types[n], ])) -
      sum(jumps[[n - 1]] * integral(horizon - times[n]))
  }
  result
}
fails(objective(model, obs, NULL, build_map(model, TRUE), derivatives = FALSE),
  "C[+][+] likelihood requires lambda_component0")
for (reduced in c(TRUE, FALSE)) {
  map <- build_map(model, reduced)
  fun <- objective(model, obs, initial, map, derivatives = FALSE)
  stopifnot(is.null(fun$gradient))
  for (point in list(map$start, map$start + .001)) {
    values <- emhawkes:::restore_parameters(map, point)
    close(fun$loglik(point), full_history(values, initial), 1e-12)
    close(fun$loglik(rev(point)), full_history(values, initial), 1e-12)
    restored_model <- model
    for (name in names(values)) slot(restored_model, name) <- values[[name]]
    close(fun$loglik(point), suppressMessages(as.numeric(logLik(restored_model, dt, types,
      marks, N0 = c(2, 3), Nc0 = c(4, 5), lambda_component0 = initial))), 1e-12)
  }
}

# Shared parameters across core callbacks and impact, with automatic initial state.
shared <- new("hspec", mu = function(param = c(q = -.5)) matrix(exp(param)),
  alpha = function(param = numeric(0)) matrix(.05), beta = 1.3,
  impact = function(param = c(q = -.5), n, N, mu, alpha, beta, ...) {
    matrix(.02 * exp(param["q"]) * (1 + N[n, 1] %% 2) + .01 * mu + .002 * alpha / beta)
  })
map <- build_map(shared, TRUE)
stopifnot(sum(names(map$start) == "q") == 1L)
fun <- objective(shared, normalize(shared, dt), NULL, map, derivatives = FALSE)
for (point in list(map$start, map$start + .01)) {
  candidate <- shared
  for (name in names(map$components)) slot(candidate, name) <- emhawkes:::restore_parameters(map, point)[[name]]
  close(fun$loglik(point), suppressMessages(as.numeric(logLik(candidate, dt))), 1e-12)
}
# A callback without local parameters also needs no Jacobian. Call it at every
# actual event, including the final event, and never at the initial row.
seen <- integer()
fixed <- h
fixed@impact <- function(n, mark, ...) { seen <<- c(seen, n); matrix(.002 * mark[n]) }
map <- build_map(fixed, TRUE)
fun <- objective(fixed, normalize(fixed, dt, mark = marks), .1, map, derivatives = FALSE)
value <- fun$loglik(map$start)
stopifnot(identical(seen, 2:length(dt)))
close(fun$loglik(map$start), value, 0)
stopifnot(identical(seen, 2:length(dt)))
empty <- objective(fixed, normalize(fixed, 0), 0, map, derivatives = FALSE)
close(empty$loglik(map$start), 0, 0)

# Analytic scores still reject feedback. Numerical observed-only callbacks
# must declare their state needs rather than implicitly reading them via ... .
for (name in c("lambda", "lambda_component", "lambda_component_n")) {
  bad <- h
  bad@impact <- function(...) matrix(0)
  formals(bad@impact) <- c(setNames(alist(x = ), name), alist(... = ))
  fails(fit_engine("cpp", bad, dt, lambda_component0 = .1, gradient = "analytic"),
    "latent intensity.*analytic")
}
for (callback in list(function(n, ...) matrix(list(...)$lambda_component_n),
                      function(n, ...) matrix(..6),
                      function(n, ...) { force(..8); matrix(0) })) {
  bad <- h
  bad@impact <- callback
  fails(fit_engine("cpp", bad, dt, lambda_component0 = .1), "latent intensity.*C[+][+]")
}
bad <- h
bad@impact <- function(n) 0
fails(fit_engine("cpp", bad, dt, lambda_component0 = .1), "plain numeric matrix")
bad@impact <- function(n) structure(matrix(0), class = c("custom", "matrix", "array"))
fails(fit_engine("cpp", bad, dt, lambda_component0 = .1), "plain numeric matrix")
bad@impact <- function(n) matrix(Inf)
map <- build_map(bad, TRUE)
fun <- objective(bad, normalize(bad, dt), .1, map, derivatives = FALSE)
stopifnot(identical(as.numeric(fun$loglik(map$start)), -Inf))
options(previous)
compiler::enableJIT(old_jit)
