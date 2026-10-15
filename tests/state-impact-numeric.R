library(emhawkes)
old_options <- options(emhawkes.engine = "r")
old_jit <- compiler::enableJIT(0)
map_for <- emhawkes:::build_parameter_map
normalize <- emhawkes:::normalize_observation
objective <- emhawkes:::make_hawkes_objective
close <- function(x, y, tolerance = 1e-10) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}
dt <- c(0, .3, 0, .7, .2, 1.1, .05)
type <- c(0, 1, 2, 2, 1, 2, 1)
mark <- c(0, 1, 3, .5, 2, 1, 2.5)
records <- list()
mu <- function(param = c(q = .2), n, Nc, ...) {
  previous <- if (n == 1L) Nc[1, ] else Nc[n - 1L, ]
  matrix(exp(param["q"]) * c(.5, .6) + .01 * previous, 2)
}
impact <- function(param = c(q = .2), n, N, Nc, mark, type, inter_arrival,
                   lambda_component_n, lambda, lambda_component, mu, alpha, beta,
                   extra = .001, ...) {
  # Both full histories include the current pre-jump row. Future rows are zero.
  stopifnot(identical(colnames(lambda), c("lambda1", "lambda2")),
    identical(colnames(lambda_component), c("lambda11", "lambda12", "lambda13", "lambda21", "lambda22", "lambda23")),
    is.function(mu))
  close(matrix(lambda_component[n, ], 2, byrow = TRUE), lambda_component_n, 0)
  close(lambda[n, ], as.vector(mu(n = n, Nc = Nc)) + rowSums(lambda_component_n), 0)
  if (n < length(mark)) {
    stopifnot(all(lambda[(n + 1):length(mark), ] == 0), all(lambda_component[(n + 1):length(mark), ] == 0))
  }
  records[[length(records) + 1L]] <<- list(n = n, current = lambda_component_n,
    lambda = lambda, components = lambda_component, local = param, N = N, Nc = Nc)
  value <- .02 * param["q"] + extra * (Nc[n, 1] - N[n - 1, 2]) - .03 * lambda_component_n +
    matrix(c(-.01, .01, .02, 0, -.02, .01), 2)
  # Local writes and retained arguments must not corrupt recursive storage.
  lambda_component_n[1, 1] <- 999
  lambda_component[n, 1] <- 999
  lambda[n, 1] <- 999
  value
}
h <- new("hspec", mu = mu, alpha = matrix(c(.08, -.02, .12, .06, .03, .09), 2),
  beta = matrix(c(.7, 1.1, 2.3, .9, 3.1, .9), 2), eta = matrix(.01, 2, 3),
  impact = impact, type_col_map = list(c(3, 1), c(3, 2)))
initial <- matrix(c(.11, -.03, .07, .13, -.02, .08), 2)
obs <- normalize(h, dt, type, mark, N0 = c(2, 3), Nc0 = c(4, 5))
reference <- function(model, point, initial = initial, reduced = TRUE) {
  map <- map_for(model, reduced)
  candidate <- model
  for (name in names(map$components)) slot(candidate, name) <- emhawkes:::restore_parameters(map, point)[[name]]
  suppressMessages(logLik(candidate, dt, type, mark, N0 = c(2, 3), Nc0 = c(4, 5), lambda_component0 = initial))
}
set.seed(519)
seed <- .Random.seed
snap <- serialize(h, NULL)
for (reduced in c(TRUE, FALSE)) {
  map <- map_for(h, reduced)
  for (point in list(map$start, map$start + .001)) {
    records <- list()
    ref <- reference(h, point, initial, reduced)
    expected <- records
    records <- list()
    fun <- objective(h, obs, initial, map, derivatives = FALSE)
    close(fun$loglik(point), ref, 1e-12)
    stopifnot(identical(records, expected)) # Includes retained snapshots and names.
    count <- length(records)
    close(fun$loglik(rev(point)), ref, 1e-12)
    stopifnot(length(records) == count)
  }
}
stopifnot(identical(seed, .Random.seed), identical(snap, serialize(h, NULL)))

# Independent full-history reference with feedback and an observed baseline.
independent <- function() {
  times <- cumsum(dt)
  horizon <- tail(times, 1)
  W <- function(t) -expm1(-h@beta * t) / h@beta
  total <- -sum(initial * W(horizon))
  jumps <- list()
  for (n in 2:length(dt)) {
    baseline <- mu(n = n, Nc = obs$Nc)
    state <- initial * exp(-h@beta * times[n])
    if (n > 2) for (k in 2:(n - 1)) state <- state + jumps[[k - 1]] * exp(-h@beta * (times[n] - times[k]))
    total <- total + log(baseline[type[n]] + sum(state[type[n], ])) - sum(baseline * dt[n])
    extra <- .02 * .2 + .001 * (obs$Nc[n, 1] - obs$N[n - 1, 2]) - .03 * state + matrix(c(-.01, .01, .02, 0, -.02, .01), 2)
    jump <- h@alpha + h@eta * (mark[n] - 1) + extra
    unused <- setdiff(1:3, h@type_col_map[[type[n]]])
    jump[, unused] <- 0
    jumps[[n - 1]] <- jump
    total <- total - sum(jump * W(horizon - times[n]))
  }
  total
}
map <- map_for(h, TRUE)
fun <- objective(h, obs, initial, map, derivatives = FALSE)
close(fun$loglik(map$start), independent(), 1e-12)

# Public hfit dispatch: keep only q free to make the short fit well specified.
a <- list(object = h, inter_arrival = dt, type = type, mark = mark,
  N0 = c(2, 3), Nc0 = c(4, 5), lambda_component0 = initial,
  fixed = names(map$start) != "q", control = list(iterlim = 3))
fit_with <- function(engine) {
  options(emhawkes.engine = engine)
  suppressWarnings(suppressMessages(do.call(hfit, c(a, list(jacobian = quote(stop("unused Jacobian")))))))
}
r <- fit_with("r")
cpp <- fit_with("cpp")
close(coef(cpp), coef(r), 1e-8)
close(logLik(cpp), logLik(r), 1e-11)
close(cpp$hessian, r$hessian, 1e-3)
stopifnot(attr(cpp, "emhawkes")$gradient == "numeric", is.null(attr(cpp, "emhawkes")$jacobian))
options(emhawkes.engine = "r")
close(suppressMessages(do.call(logLik, c(list(object = fitted_hspec(cpp), inter_arrival = dt, type = type, mark = mark), fit_initial_state(cpp)))), logLik(cpp), 1e-11)

# Current-state-only callbacks do not allocate histories. Also exercise automatic
# initialization and observed baselines without impacts.
unary <- new("hspec", mu = function(param = c(m = .4), n, N, ...) matrix(param + .01 * N[max(1, n - 1), 1]),
  alpha = .1, beta = 1.2, impact = function(n, lambda_component_n, ...) .02 * lambda_component_n)
for (model in list(unary, { x <- unary; x@impact <- NULL; x })) {
  m <- map_for(model, TRUE)
  f <- objective(model, normalize(model, dt), NULL, m, derivatives = FALSE)
  close(f$loglik(m$start), suppressMessages(logLik(model, dt)), 1e-12)
  empty <- objective(model, normalize(model, 0), 0, m, derivatives = FALSE)
  close(empty$loglik(m$start), 0, 0)
}

# Retaining lazy state arguments and forcing GC must preserve their old values.
lazy <- list()
retained <- unary
retained@impact <- function(n, lambda_component, ...) {
  lazy[[length(lazy) + 1L]] <<- function() list(n = n, state = lambda_component)
  invisible(gc(FALSE))
  matrix(0L)
}
m <- map_for(retained, TRUE)
lazy <- list()
ref <- suppressMessages(logLik(retained, dt, lambda_component0 = .1))
expected <- lapply(lazy, function(f) f())
lazy <- list()
f <- objective(retained, normalize(retained, dt), .1, m, derivatives = FALSE)
close(f$loglik(m$start), ref, 1e-12)
stopifnot(identical(lapply(lazy, function(f) f()), expected))
empty_model <- unary
empty_model@mu <- function(n, ...) { if (n != 1L) stop("unexpected terminal baseline"); matrix(.4) }
m <- map_for(empty_model, TRUE)
f <- objective(empty_model, normalize(empty_model, 0), 0, m, derivatives = FALSE)
close(f$loglik(m$start), 0, 0)

# Explicit limitations and native validation cannot produce a wrong analytic score.
fails(objective(h, obs, initial, map, jacobian = list(), derivatives = TRUE), "full intensity histories")
bad <- unary
bad@mu <- function(param = c(m = .4), lambda) matrix(param)
fails(objective(bad, normalize(bad, dt), .1, map_for(bad, TRUE), derivatives = FALSE), "latent intensity")
bad <- unary
bad@impact <- function(lambda_component_n) stop("callback sentinel")
f <- objective(bad, normalize(bad, dt), .1, map_for(bad, TRUE), derivatives = FALSE)
fails(f$loglik(map_for(bad, TRUE)$start), "callback sentinel")
for (value in list(0, matrix(0, 2), structure(matrix(0), class = c("custom", "matrix", "array")))) {
  bad@impact <- local({ v <- value; function(lambda_component_n) v })
  m <- map_for(bad, TRUE)
  f <- objective(bad, normalize(bad, dt), .1, m, derivatives = FALSE)
  fails(f$loglik(m$start), "plain numeric matrix")
}
bad@impact <- function(lambda_component_n) matrix(Inf)
m <- map_for(bad, TRUE)
f <- objective(bad, normalize(bad, dt), .1, m, derivatives = FALSE)
stopifnot(identical(as.numeric(f$loglik(m$start)), -Inf))
options(old_options)
compiler::enableJIT(old_jit)
