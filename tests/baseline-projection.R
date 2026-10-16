library(emhawkes)
old <- options(emhawkes.engine = "auto")
close <- function(x, y, tol = 1e-12) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
fails <- function(x, pattern) {
  e <- tryCatch(force(x), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
pack <- emhawkes:::hawkes_jacobian_path_cpp
native <- emhawkes:::hawkes_loglik_cpp

# Compare projected scores with the original full Jacobian path, independently
# of the R planner. Column order changes at each event; callbacks retain n lazily.
for (d in c(1L, 2L, 3L, 4L, 8L)) {
  dt <- c(0, .3, 0, .7, .2, 1.1, .05, .4)
  type <- c(0L, (seq_len(length(dt) - 1L) - 1L) %% d + 1L)
  p <- c(a = .3, b = -.2, c = .5, z = 0)
  value <- function(n) outer(seq_len(d), seq_along(p), function(i, j) sin(i + j + n) / 10)
  calls <- integer(); retained <- list()
  jac <- function(param, n) {
    calls <<- c(calls, n)
    retained[[length(retained) + 1L]] <<- function() n
    J <- value(n); colnames(J) <- names(param)
    J[, if (n %% 2L) rev(names(param)) else names(param), drop = FALSE]
  }
  result <- pack(jac, list(param = p, n = 1L), d, dt, type, names(p))
  gc(FALSE)
  stopifnot(identical(calls, seq_along(dt)), identical(vapply(retained, function(f) f(), 1L), seq_along(dt)))
  dense <- vapply(seq_along(dt), function(n) as.vector(value(n)), numeric(d * length(p)))
  if (d > 2L) {
    stopifnot(is.list(result))
    close(result$initial, value(1), 0)
    close(result$event[, -1L], vapply(seq_len(length(dt) - 1L) + 1L,
      function(n) value(n)[type[n], ], numeric(length(p))), 0)
    restored <- result$sum * 2^result$exponent
    close(restored[, -1L], vapply(seq_len(length(dt) - 1L) + 1L,
      function(n) colSums(value(n)), numeric(length(p))), 1e-14)
  } else stopifnot(identical(result, dense))
  k <- d + 1L; m <- d * k; np <- length(p)
  sensitivity <- list(mu = value(1), alpha = matrix(.003, m, np),
    beta = matrix(.002, m, np), eta = matrix(0, m, np),
    initial = matrix(.001, m, np), names = names(p), baseline = dense, stateful = FALSE)
  args <- list(dt = dt, type = type, mark = rep(1, length(dt)), mu = rep(.8, d),
    alpha = matrix(.03, d, k), beta = matrix(seq(.7, 2.3, length.out = m), d),
    eta = NULL, initial = matrix(.1, d, k), mapping = lapply(seq_len(d), function(i) c(i - 1L, k - 1L)),
    infer = FALSE, gradient = TRUE, baseline = matrix(.8, d, length(dt)), sensitivity = sensitivity)
  reference <- do.call(native, args)
  args$sensitivity$baseline <- result
  actual <- do.call(native, args)
  stopifnot(identical(actual$loglikelihood, reference$loglikelihood))
  close(actual$gradient, reference$gradient, 1e-13)
  # Empty and short paths never evaluate an extra Jacobian index.
  for (size in 1:3) {
    calls <- integer()
    short <- pack(jac, list(param = p, n = 1L), d, dt[seq_len(size)], type[seq_len(size)], names(p))
    stopifnot(identical(calls, seq_len(size)))
    if (size == 1L) stopifnot(is.matrix(short))
  }
}
# Scaling must retain subnormal derivatives and sums larger than double range.
p <- c(small = 0, large = 0)
jac <- function(param, n) matrix(rep(c(5e-324, 1e308), each = 4), 4,
  dimnames = list(NULL, names(param)))
z <- pack(jac, list(param = p, n = 1L), 4L, c(0, rep(1, 7)), c(0L, rep(1L, 7)), names(p))
stopifnot(z$sum[1, 2] > 0, z$exponent[1, 2] < -1000,
  z$sum[2, 2] >= .5, z$exponent[2, 2] > 1024)
bad <- function(param, n) {J <- jac(param, n); if (n == 3L) J[1] <- NA_real_; J}
fails(pack(bad, list(param = p, n = 1L), 4L, rep(0, 8), c(0L, rep(1L, 7)), names(p)), "finite")

# Closed-form Poisson scores also check the integral reconstruction at extreme
# scales. Intermediate row sums exceed double range in the large case.
for (scale in c("small", "large")) {
  derivative <- if (scale == "small") 5e-324 else 1e308
  interval <- if (scale == "small") 1e308 else 1e-308
  level <- if (scale == "small") 1e-308 else 1e308
  dt_extreme <- c(0, rep(interval, 7))
  jac_extreme <- local({v <- derivative; function(param, n) matrix(v, 4, 1, dimnames = list(NULL, "x"))})
  projected <- pack(jac_extreme, list(param = c(x = 0), n = 1L), 4L,
    dt_extreme, c(0L, rep(1L, 7)), "x")
  sens <- list(mu = matrix(derivative, 4, 1), alpha = matrix(0, 16, 1),
    beta = matrix(0, 16, 1), eta = matrix(0, 16, 1), initial = matrix(0, 16, 1),
    names = "x", baseline = projected, stateful = FALSE)
  actual <- native(dt_extreme, c(0L, rep(1L, 7)), rep(1, 8), rep(level, 4),
    matrix(0, 4, 4), matrix(1, 4, 4), NULL, matrix(0, 4, 4),
    lapply(0:3, function(i) i), FALSE, gradient = TRUE, baseline = matrix(level, 4, 8), sensitivity = sens)
  expected <- 7 * (derivative / level - 4 * (interval * derivative))
  stopifnot(is.finite(actual$gradient), abs(actual$gradient - expected) <= abs(expected) * 1e-13)
}

# Public observed-baseline models exercise shared parameters, automatic initial
# states, current-state feedback, cold scores and the numerical R reference.
mu <- function(param = c(m = -.8, a = -2), n, N, alpha, ...) {
  matrix(exp(param["m"]) * (1:4) / 4 * (1 + N[max(1L, n - 1L), ] / 20) + diag(alpha) / 10, 4)
}
alpha <- function(param = c(a = -2)) diag(exp(param), 4) + matrix(.01, 4, 4)
h <- new("hspec", mu = mu, alpha = alpha, beta = matrix(seq(.7, 2.2, length.out = 16), 4),
  impact = function(param = c(q = .02), lambda_component_n) matrix(param["q"] + .01 * mean(lambda_component_n), 4, 4))
j <- list(mu = function(param, n, N, ...) {
  J <- matrix(0, 4, length(param), dimnames = list(NULL, names(param)))
  J[, "m"] <- exp(param["m"]) * (1:4) / 4 * (1 + N[max(1L, n - 1L), ] / 20)
  J[, "a"] <- exp(param["a"]) / 10
  J[, rev(names(param)), drop = FALSE]
}, alpha = function(param) matrix(as.vector(diag(exp(param), 4)), dimnames = list(NULL, "a")),
impact = list(param = function(param) {
  J <- matrix(0, 16, length(param), dimnames = list(NULL, names(param))); J[, "q"] <- 1; J
}, state = matrix(.01 / 16, 16, 16)))
dt <- c(0, .2, .5, 0, .7, .3, 1, .1)
type <- c(0, 1, 4, 2, 3, 2, 4, 1)
obs <- emhawkes:::normalize_observation(h, dt, type)
for (reduced in c(TRUE, FALSE)) for (initial in list(NULL, .1)) {
  map <- emhawkes:::build_parameter_map(h, reduced)
  f <- emhawkes:::make_hawkes_objective(h, obs, initial, map, jacobian = j)
  reference <- function(p) {
    model <- h; values <- emhawkes:::restore_parameters(map, p)
    for (name in names(values)) slot(model, name) <- values[[name]]
    op <- options(emhawkes.engine = "r"); on.exit(options(op))
    as.numeric(suppressMessages(logLik(model, dt, type, lambda_component0 = initial)))
  }
  for (p in list(map$start, map$start + .01 * cos(seq_along(map$start)))) {
    g <- f$score(p); close(f$loglik(p), reference(p), 1e-12)
    close(f$gradient(p), g, 0)
    for (step in c(1e-4, 2e-5)) {
      fd <- vapply(seq_along(p), function(k) {
        a <- b <- p; delta <- step * max(1, abs(p[k])); a[k] <- a[k] + delta; b[k] <- b[k] - delta
        (reference(a) - reference(b)) / (2 * delta)
      }, 0)
      close(g, fd, 5e-7)
    }
  }
}
# The complete mu-Jacobian sequence still precedes every impact replay.
events <- character(); saved <- list()
j$mu <- local({base <- j$mu; function(param, n, N, ...) {
  events <<- c(events, paste0("M", n)); saved[[length(saved) + 1L]] <<- function() n
  base(param, n, N)
}})
h@impact <- function(param = c(q = .02), n, lambda_component_n) {
  events <<- c(events, paste0("I", n)); matrix(param["q"] + .01 * mean(lambda_component_n), 4, 4)
}
j$impact$param <- local({base <- j$impact$param; function(param, n) {
  events <<- c(events, paste0("P", n)); base(param)
}})
j$impact$state <- function(n) {events <<- c(events, paste0("S", n)); matrix(.01 / 16, 16, 16)}
map <- emhawkes:::build_parameter_map(h, TRUE)
f <- emhawkes:::make_hawkes_objective(h, obs, .1, map, jacobian = j)
invisible(f$score(map$start))
expected <- c(paste0("M", seq_along(dt)), as.vector(rbind(paste0("I", 2:length(dt)),
  paste0("P", 2:length(dt)), paste0("S", 2:length(dt)))))
stopifnot(identical(events, expected), identical(vapply(saved, function(f) f(), 1L), seq_along(dt)))
invisible(f$loglik(map$start)); invisible(f$gradient(map$start))
stopifnot(identical(events, expected))
options(old)
