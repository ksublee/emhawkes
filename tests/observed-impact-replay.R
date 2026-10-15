library(emhawkes)
old_jit <- compiler::enableJIT(0L)
old <- options(emhawkes.engine = "auto")
evaluate <- function(engine, object, dt, infer = FALSE, initial = .1, ...) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(logLik(object, dt, infer = infer, lambda_component0 = initial, ...))
}
close <- function(a, b, tolerance = 1e-10) stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance)))
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), ignore.case = TRUE))
}
native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
calls <- 0L
assignInNamespace("hawkes_loglik_cpp", function(...) { calls <<- calls + 1L; native(...) }, "emhawkes")
# Independent shot-noise sums fix the intensity, compensator and residuals.
h <- new("hspec", mu = .4, alpha = .05, beta = 1.3, eta = .01,
  impact = function(n, N, Nc, mark) matrix(.03 * N[n, 1] + .02 * Nc[n, 1] + .01 * (mark[n] - 1)))
dt <- c(0, .2, .3, .4, .1)
marks <- c(0, 1, 2, 1, 3)
times <- cumsum(dt)
jumps <- .05 + .01 * (marks[-1] - 1) + .03 * (2 + seq_len(4)) +
  .02 * (3 + cumsum(marks[-1])) + .01 * (marks[-1] - 1)
pre <- vapply(2:5, function(n) {
  .1 * exp(-1.3 * times[n]) + if (n > 2) sum(jumps[seq_len(n - 2)] * exp(-1.3 * (times[n] - times[2:(n - 1)]))) else 0
}, 0)
compensator <- function(t) .4 * t + .1 * (-expm1(-1.3 * t)) / 1.3 +
  sum(jumps * (-expm1(-1.3 * pmax(0, t - times[-1])))) / 1.3
expected <- sum(log(.4 + pre)) - compensator(tail(times, 1))
for (engine in c("r", "auto", "cpp")) {
  calls <- 0L
  value <- evaluate(engine, h, dt, infer = TRUE, mark = marks, N0 = 2, Nc0 = 3)
  close(as.numeric(value$loglikelihood), expected)
  close(as.vector(value$lambda[-1, ]), .4 + pre)
  close(as.vector(value$rambda_component[-1, ]), pre + jumps)
  close(value$res_process1, diff(vapply(times[-1], compensator, 0)))
  stopifnot(calls == if (engine == "r") 0L else 1L)
}
# Same default names in different functions remain independent during replay.
independent <- new("hspec", mu = .3,
  alpha = function(param = c(a = .15)) matrix(param), beta = 1.2,
  impact = function(param = c(a = .04)) matrix(param))
combined <- new("hspec", mu = .3, alpha = .19, beta = 1.2)
for (infer in c(FALSE, TRUE)) close(evaluate("cpp", independent, dt, infer),
                                  evaluate("r", combined, dt, infer))
# Matrix order, overlapping/reordered mappings, marks and supplied initial counts.
mu <- function(param = c(m = .3)) matrix(c(param, 1.2 * param), 2)
alpha <- function(param = c(a = .04)) matrix(param, 2, 4)
beta <- function(param = c(b = 1.1)) matrix(param + seq(0, .7, length.out = 8), 2)
impact <- function(param = c(a = .005), n, N, Nc, mark, type, inter_arrival, mu, alpha, beta, ...) {
  alpha * param + matrix(param * (sum(N[n, ]) + .1 * sum(Nc[n, ]) + mark[n] +
    type[n] * inter_arrival[n] + sum(mu) + sum(beta)) / 100, 2, 4)
}
h <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = matrix(.01, 2, 4),
  impact = impact, type_col_map = list(c(3L, 1L), c(4L, 1L)))
a <- list(type = c(0, 1, 2, 2, 1), mark = marks, N0 = c(2, 3), Nc0 = c(4, 5))
snapshot <- serialize(h, NULL)
for (infer in c(FALSE, TRUE)) {
  reference <- do.call(evaluate, c(list("r", h, dt, infer), a))
  for (engine in c("auto", "cpp")) {
    calls <- 0L
    close(do.call(evaluate, c(list(engine, h, dt, infer), a)), reference)
    stopifnot(calls == 1L)
  }
}
stopifnot(identical(snapshot, serialize(h, NULL)))
# Empty paths do not call impact, and integer matrices are supported.
empty <- new("hspec", mu = .3, alpha = 0, beta = 1,
  impact = function(n) stop("impact must not run on an empty path"))
close(evaluate("cpp", empty, 0, TRUE), evaluate("r", empty, 0, TRUE))
empty@impact <- function() matrix(1L)
close(evaluate("cpp", empty, dt, TRUE), evaluate("r", empty, dt, TRUE))
# Hidden state access through dots triggers R before native invocation. Preserve
# the already-computed automatic initial state on this fallback.
hidden <- new("hspec", mu = function(param = c(m = .3)) matrix(param), alpha = .1, beta = 1.2,
  impact = function(n, ...) { if (n >= 3) list(...)$lambda_component_n * .01 else matrix(.02) })
for (initial in list(.1, NULL)) {
  calls <- 0L
  close(evaluate("auto", hidden, dt, TRUE, initial), evaluate("r", hidden, dt, TRUE, initial))
  stopifnot(calls == 0L)
}
fails(evaluate("cpp", hidden, dt), "latent intensity")
explicit <- hidden
explicit@impact <- function(lambda_component_n) .01 * lambda_component_n
calls <- 0L
close(evaluate("auto", explicit, dt, TRUE), evaluate("r", explicit, dt, TRUE))
stopifnot(calls == 1L)
close(evaluate("cpp", explicit, dt), evaluate("r", explicit, dt))
# Nonfinite unselected entries can be ignored by R; preserve that behavior
# instead of turning a finite R likelihood into native -Inf.
unused <- new("hspec", mu = c(.3, .4), alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2),
  impact = function(n, type) { x <- matrix(NA_real_, 2, 2); x[, type[n]] <- .01; x })
calls <- 0L
close(evaluate("auto", unused, dt, TRUE, type = c(0, 1, 2, 1, 2)),
      evaluate("r", unused, dt, TRUE, type = c(0, 1, 2, 1, 2)))
stopifnot(calls == 0L)
fails(evaluate("cpp", unused, dt, type = c(0, 1, 2, 1, 2)), "finite plain numeric")
failures <- 0L
failed <- empty
failed@impact <- function(n) { failures <<- failures + 1L; stop("observed replay sentinel") }
fails(evaluate("auto", failed, dt), "observed replay sentinel")
stopifnot(failures == 1L)
# Replay must not draw simulation marks or advance RNG state.
empty@rmark <- function(...) stop("do not draw marks")
empty@rresidual <- function(...) stop("do not draw residuals")
set.seed(56); seed <- .Random.seed
invisible(evaluate("cpp", empty, dt, TRUE))
stopifnot(identical(seed, .Random.seed))
assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
options(old)
compiler::enableJIT(old_jit)
