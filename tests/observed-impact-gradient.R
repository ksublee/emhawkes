# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)
# JIT compilation mutates closure bytecode; exclude it from serialization checks.
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
central <- function(fn, x, step = 2e-5) vapply(seq_along(x), function(k) {
  up <- down <- x
  delta <- step * max(1, abs(x[k]))
  up[k] <- up[k] + delta
  down[k] <- down[k] - delta
  (fn(up) - fn(down)) / (2 * delta)
}, numeric(1))
previous <- options(emhawkes.engine = "auto")

# Independent full-history likelihood, including every selected historical jump.
history_ll <- function(dt, type, mark, mu, alpha, beta, eta, initial, mapping, impact, ...) {
  times <- cumsum(dt)
  horizon <- tail(times, 1)
  W <- function(rate, time) -expm1(-rate * time) / rate
  result <- -sum(mu) * horizon - sum(initial * W(beta, horizon))
  for (n in seq_len(length(dt) - 1L) + 1L) {
    state <- initial * exp(-beta * times[n])
    for (k in seq_len(n - 2L) + 1L) {
      columns <- mapping[[type[k]]] + 1L
      jump <- alpha + eta * (mark[k] - 1) + matrix(impact[, k], nrow(beta))
      state[, columns] <- state[, columns, drop = FALSE] +
        jump[, columns, drop = FALSE] * exp(-beta[, columns, drop = FALSE] * (times[n] - times[k]))
    }
    result <- result + log(mu[type[n]] + sum(state[type[n], ]))
    columns <- mapping[[type[n]]] + 1L
    jump <- alpha + eta * (mark[n] - 1) + matrix(impact[, n], nrow(beta))
    result <- result - sum(jump[, columns, drop = FALSE] *
      W(beta[, columns, drop = FALSE], horizon - times[n]))
  }
  result
}
args <- list(dt = c(0, .3, 0, .7, .2, 1.1, .05),
  type = c(0L, 1L, 2L, 2L, 1L, 2L, 1L), mark = c(0, 1, 3, .5, 2, 1, 2.5),
  mu = c(.8, .9), alpha = matrix(c(.08, -.02, .12, .06, .03, .09), 2),
  beta = matrix(c(.7, 1.1, 2.3, .9, 3.1, .9), 2),
  eta = matrix(c(-.01, .02, .03, -.005, .015, .01), 2),
  initial = matrix(c(.11, -.03, .07, .13, -.02, .08), 2),
  mapping = list(c(2L, 0L), c(2L, 1L)), infer = FALSE,
  impact = matrix(seq(-.03, .1, length.out = 42), 6))
snapshot <- serialize(args, NULL)
set.seed(98)
seed <- .Random.seed
native <- function(x, gradient = TRUE) do.call(emhawkes:::hawkes_loglik_cpp, c(x, list(gradient = gradient)))
result <- native(args)
close(result$loglikelihood, do.call(history_ll, args), 1e-12)
close(result$loglikelihood, native(args, FALSE)$loglikelihood, 0)
for (name in c("mu", "alpha", "beta", "eta", "initial", "impact")) {
  reference <- function(value) {
    candidate <- args
    candidate[[name]][] <- value
    do.call(history_ll, candidate)
  }
  score <- if (name == "impact") result$impact_score else result$score[[name]]
  for (step in c(1e-4, 2e-5)) close(score, central(reference, args[[name]], step), 3e-7)
}
# Initial-row and terminal-event jumps do not affect this observation window.
close(result$impact_score[, c(1, ncol(args$impact))], 0 * args$impact[, c(1, ncol(args$impact))], 0)
for (n in 2:length(args$dt)) {
  ignored <- setdiff(seq_len(nrow(args$impact)),
    as.vector(outer(seq_len(nrow(args$beta)), nrow(args$beta) * args$mapping[[args$type[n]]], "+")))
  close(result$impact_score[ignored, n], rep(0, length(ignored)), 0)
}
inferred <- args
inferred$infer <- TRUE
close(native(inferred)$impact_score, result$impact_score, 0)
stopifnot(identical(seed, .Random.seed), identical(snapshot, serialize(args, NULL)))
bad <- args
bad$impact <- matrix(0, 1, 1)
fails(native(bad), "impact dimensions")
bad$impact <- args$impact
bad$impact[1] <- Inf
fails(native(bad), "impact value")

# Shared parameters and all core dependencies contribute to the impact's total
# Jacobian. Its param argument uses the full fitted vector; the impact uses local param.
mu <- function(param = c(k = -.7)) matrix(c(exp(param["k"]), .8), 2)
alpha <- function(param = c(a = -2, b = -3)) matrix(c(exp(param["a"]), 0, exp(param["b"]), exp(param["a"])), 2)
beta <- function(param = c(k = -.7, d = .3)) matrix(c(exp(param["d"]), 1 + exp(param["k"]), exp(param["d"]), 2 * exp(param["d"])), 2)
weight <- function(n, N, Nc, mark, type, inter_arrival) {
  .2 + log1p(N[n, 1] + Nc[n, 2]) / (1 + N[n - 1, 2]) +
    .05 * mark[n] + .01 * type[n] + .02 * inter_arrival[n]
}
impact <- function(param = c(q = .15, k = -.7), n, N, Nc, mark, type,
                   inter_arrival, mu, alpha, beta, ...) {
  w <- weight(n, N, Nc, mark, type, inter_arrival)
  param["q"] * w * alpha + .01 * mu %*% matrix(1, 1, 2) +
    .002 * beta^2 + .003 * exp(param["k"])
}
impact_jac <- function(param, n, N, Nc, mark, type, inter_arrival, mu, alpha, beta, ...) {
  w <- weight(n, N, Nc, mark, type, inter_arrival)
  out <- matrix(0, 4, length(param), dimnames = list(NULL, names(param)))
  out[, "q"] <- w * as.vector(alpha)
  out[, "a"] <- param["q"] * w * c(1, 0, 0, 1) * exp(param["a"])
  out[, "b"] <- param["q"] * w * c(0, 0, 1, 0) * exp(param["b"])
  out[, "k"] <- .01 * c(1, 0, 1, 0) * exp(param["k"]) +
    .004 * as.vector(beta) * c(0, 1, 0, 0) * exp(param["k"]) + .003 * exp(param["k"])
  out[, "d"] <- .004 * as.vector(beta) * c(1, 0, 1, 2) * exp(param["d"])
  out[, rev(names(param)), drop = FALSE]
}
jac <- list(mu = function(param) matrix(c(exp(param), 0), dimnames = list(NULL, "k")),
  alpha = function(param) cbind(a = c(1, 0, 0, 1) * exp(param["a"]), b = c(0, 0, 1, 0) * exp(param["b"])),
  beta = function(param) cbind(k = c(0, 1, 0, 0) * exp(param["k"]), d = c(1, 0, 1, 2) * exp(param["d"])),
  impact = impact_jac)
h <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = matrix(.01, 2, 2), impact = impact)
check <- function(model, jac, initial, dt = args$dt, type = args$type, mark = args$mark) {
  snap <- serialize(model, NULL)
  seed <- .Random.seed
  obs <- normalize(model, dt, type, mark, N0 = c(2, 3), Nc0 = c(4, 5))
  for (reduced in c(TRUE, FALSE)) {
    map <- build_map(model, reduced)
    j <- if (is.function(jac)) jac(map) else jac
    fun <- objective(model, obs, initial, map, jacobian = j)
    value_only <- objective(model, obs, initial, map, derivatives = FALSE)
    stopifnot(is.null(value_only$gradient))
    reference <- function(param) {
      candidate <- model
      values <- emhawkes:::restore_parameters(map, param)
      for (name in names(values)) slot(candidate, name) <- values[[name]]
      old <- options(emhawkes.engine = "r")
      on.exit(options(old))
      suppressMessages(as.numeric(logLik(candidate, dt, type, mark,
        N = obs$N, Nc = obs$Nc, N0 = obs$N0, Nc0 = obs$Nc0, lambda_component0 = initial)))
    }
    theta <- map$start
    direction <- rep(c(.2, -.3, .4), length.out = length(theta))
    for (point in list(theta, theta + .025 * direction)) {
      close(fun$loglik(point), reference(point), 1e-12)
      close(value_only$loglik(point), fun$loglik(point), 0)
      score <- fun$gradient(point)
      stopifnot(identical(names(score), names(theta)))
      for (step in c(1e-4, 2e-5)) close(score, central(reference, point, step), 4e-7)
      close(sum(score * direction), (reference(point + 2e-5 * direction) -
        reference(point - 2e-5 * direction)) / 4e-5, 4e-7)
      close(fun$gradient(rev(point)), score, 0)
    }
    close(fun$loglik(theta), reference(theta), 1e-12)
  }
  stopifnot(identical(snap, serialize(model, NULL)), identical(seed, .Random.seed))
}
check(h, jac, .1)
check(h, jac, NULL)
h@type_col_map <- list(c(2, 1), 2)
check(h, jac, matrix(c(.1, -.01, .05, .2), 2))

# A nonsquare mapping and a constant Jacobian with a fixed signed jump.
linear <- new("hspec", mu = c(.8, .9), alpha = matrix(.03, 2, 3), beta = matrix(c(.5, 1, 2, 3, 4, 5), 2),
  type_col_map = list(c(3, 1), c(2, 3)), impact = function(param = c(q = .1), n, mark) {
    matrix(c(param["q"], -.01, param["q"], 0, 0, param["q"]), 2)
  })
linear_jac <- function(map) {
  j <- matrix(0, 6, length(map$start), dimnames = list(NULL, names(map$start)))
  j[, "q"] <- c(1, 0, 1, 0, 0, 1)
  list(impact = j)
}
check(linear, linear_jac, .1)
fixed_impact <- new("hspec", mu = .8, alpha = .1, beta = 1.2,
  impact = function(n, mark, ...) matrix(.02 * mark[n]))
fixed_jac <- function(map) list(impact = matrix(0, 1, length(map$start), dimnames = list(NULL, names(map$start))))
# Use a univariate observation directly for this case.
obs <- normalize(fixed_impact, args$dt, mark = args$mark)
map <- build_map(fixed_impact, TRUE)
fun <- objective(fixed_impact, obs, .1, map, jacobian = fixed_jac(map))
reference <- function(param) {
  candidate <- fixed_impact
  values <- emhawkes:::restore_parameters(map, param)
  for (name in names(values)) slot(candidate, name) <- values[[name]]
  suppressMessages(as.numeric(logLik(candidate, args$dt, mark = args$mark, lambda_component0 = .1)))
}
close(fun$gradient(map$start), central(reference, map$start), 3e-7)
empty <- objective(fixed_impact, normalize(fixed_impact, 0), 0, map, jacobian = fixed_jac(map))
close(empty$loglik(map$start), 0, 0)
close(empty$gradient(map$start), c(0, 0, 0), 0)

# Malformed derivatives and forbidden latent-state access fail explicitly.
make <- function(model = h, j = jac) objective(model, normalize(model, args$dt, args$type, args$mark),
  .1, build_map(model, TRUE), jacobian = j)
fails(make(j = jac[names(jac) != "impact"]), "requires jacobian.impact")
bad <- jac
bad$impact <- matrix(0, 1, 1)
fails(make(j = bad), "one row")
bad$impact <- function(param, ...) matrix(0, 4, length(param), dimnames = list(NULL, rep("wrong", length(param))))
fails(make(j = bad)$gradient(build_map(h, TRUE)$start), "column names")
bad$impact <- function(param, ...) matrix(Inf, 4, length(param), dimnames = list(NULL, names(param)))
fails(make(j = bad)$gradient(build_map(h, TRUE)$start), "finite")
bad_model <- h
bad_model@impact <- function(lambda_component_n, ...) stop("must not be called")
fails(make(bad_model), "latent intensity")
bad_model@impact <- function(n, ...) { dots <- list(...); matrix(dots$lambda_component_n, 2, 2) }
fails(make(bad_model)$loglik(build_map(bad_model, TRUE)$start), "latent intensity")
bad_model@impact <- function(n, ...) matrix(..6, 2, 2) # lambda's existing ... position
fails(make(bad_model)$loglik(build_map(bad_model, TRUE)$start), "latent intensity")
bad_model@impact <- function(n, N) matrix(1, 1, 1)
fails(make(bad_model)$loglik(build_map(bad_model, TRUE)$start), "kernel dimensions")
bad_model@impact <- function(n, N) matrix(Inf, 2, 2)
invalid <- make(bad_model)
stopifnot(identical(as.numeric(invalid$loglik(build_map(bad_model, TRUE)$start)), -Inf))

# Cache reuse must not rerun either callback for the same named parameter vector.
counts <- c(value = 0L, jacobian = 0L)
counted <- h
counted@impact <- function(param = c(q = .15, k = -.7), n, N, Nc, mark, type, inter_arrival, mu, alpha, beta, ...) {
  counts["value"] <<- counts["value"] + 1L
  impact(param, n, N, Nc, mark, type, inter_arrival, mu, alpha, beta)
}
cj <- jac
cj$impact <- function(param, n, N, Nc, mark, type, inter_arrival, mu, alpha, beta, ...) {
  counts["jacobian"] <<- counts["jacobian"] + 1L
  impact_jac(param, n, N, Nc, mark, type, inter_arrival, mu, alpha, beta)
}
fun <- make(counted, cj)
point <- build_map(counted, TRUE)$start
invisible(fun$loglik(point))
stopifnot(counts["value"] == length(args$dt) - 1L, counts["jacobian"] == 0L)
invisible(fun$gradient(rev(point)))
seen <- counts
invisible(fun$loglik(point))
invisible(fun$gradient(point))
stopifnot(identical(counts, seen), all(seen == length(args$dt) - 1L))
invisible(fun$loglik(point + .001))
stopifnot(counts["value"] == 2L * seen["value"], counts["jacobian"] == seen["jacobian"])
invisible(fun$gradient(point + .001))
stopifnot(all(counts == 2L * seen))
# One free impact amplitude has a strictly concave likelihood with an
# independently computed score, root, and information matrix.
amplitude <- new("hspec",
  mu = function(param = numeric(0)) matrix(.4),
  alpha = function(param = numeric(0)) matrix(0),
  beta = function(param = numeric(0)) matrix(1.5),
  impact = function(param = c(a = .3), n, N, ...) matrix(param["a"] * (1 + .2 * (N[n, 1] %% 2))))
amplitude_jac <- list(impact = function(param, n, N, ...) matrix(1 + .2 * (N[n, 1] %% 2),
  dimnames = list(NULL, names(param))))
set.seed(9462)
path <- hsim(amplitude, size = 301, lambda_component0 = .1, N0 = 3, Nc0 = 4)
times <- cumsum(path$inter_arrival)
weights <- 1 + .2 * ((3 + seq_along(times) - 1) %% 2)
z <- vapply(seq_along(times)[-1], function(n) {
  earlier <- seq_len(n - 2L) + 1L
  sum(weights[earlier] * exp(-1.5 * (times[n] - times[earlier])))
}, numeric(1))
base <- .4 + .1 * exp(-1.5 * times[-1])
integral <- sum(weights[-1] * -expm1(-1.5 * (tail(times, 1) - times[-1])) / 1.5)
score <- function(a) sum(z / (base + a * z)) - integral
root <- uniroot(score, c(1e-8, 10), tol = 1e-12)$root
expected_variance <- 1 / sum((z / (base + root * z))^2)
fit_one <- function(mode, ...) suppressMessages(hfit(amplitude, path$inter_arrival,
  lambda_component0 = .1, N0 = 3, Nc0 = 4, gradient = mode,
  jacobian = amplitude_jac, method = "NR", ...))
numeric_fit <- fit_one("numeric")
analytic_fit <- fit_one("analytic", final_hessian = "stable")
close(coef(analytic_fit), root, 2e-7)
close(coef(numeric_fit), root, 2e-7)
close(vcov(analytic_fit), expected_variance, 2e-6)
close(logLik(analytic_fit), logLik(numeric_fit), 1e-10)
saved <- tempfile(fileext = ".rds")
saveRDS(analytic_fit, saved)
reloaded <- readRDS(saved)
unlink(saved)
stopifnot(is.function(fitted_hspec(reloaded)@impact),
          identical(attr(reloaded, "emhawkes")$jacobian, amplitude_jac))
close(logLik(reloaded), suppressMessages(do.call(logLik, c(list(object = fitted_hspec(reloaded),
  inter_arrival = path$inter_arrival), fit_initial_state(reloaded)))), 1e-12)
fixed_model <- amplitude
formals(fixed_model@impact)$param <- c(a = .3, unused = .2)
fixed_jac <- list(impact = function(param, n, N, ...) {
  out <- matrix(0, 1, length(param), dimnames = list(NULL, names(param)))
  out[, "a"] <- 1 + .2 * (N[n, 1] %% 2)
  out
})
fixed_fit <- suppressMessages(hfit(fixed_model, path$inter_arrival,
  lambda_component0 = .1, N0 = 3, Nc0 = 4, fixed = c(FALSE, TRUE), method = "NR",
  gradient = "analytic", jacobian = fixed_jac, final_hessian = "stable"))
close(coef(fixed_fit), c(root, .2), 2e-7)
close(diag(vcov(fixed_fit)), c(expected_variance, 0), 2e-6)
# A supplied total gradient continues to override the impact adapter.
options(emhawkes.engine = "r")
overridden <- suppressMessages(hfit(amplitude, path$inter_arrival, lambda_component0 = .1,
  N0 = 3, Nc0 = 4, gradient = "analytic", grad = function(param) score(param[1]), method = "NR"))
close(coef(overridden), root, 2e-7)
options(previous)
invisible(compiler::enableJIT(old_jit))

}
