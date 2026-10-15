library(emhawkes)
old <- options(emhawkes.engine = "auto")
old_jit <- compiler::enableJIT(0L)
evaluate <- function(engine, object, dt, infer = TRUE, initial = .1, ...) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(logLik(object, dt, infer = infer, lambda_component0 = initial, ...))
}
close <- function(a, b, tolerance = 1e-11) stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance)))
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), ignore.case = TRUE))
}
native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
native_calls <- 0L
assignInNamespace("hawkes_loglik_cpp", function(...) { native_calls <<- native_calls + 1L; native(...) }, "emhawkes")
# An independent piecewise-constant baseline fixes pre/post intensity, likelihood
# and residuals, including the terminal next baseline after the final event.
mu_calls <- integer()
mu <- function(param = c(m = .2), n, N, ...) {
  mu_calls <<- c(mu_calls, n)
  matrix(param + .01 * N[max(1L, n - 1L), ], 2)
}
h <- new("hspec", mu = mu, alpha = matrix(0, 2, 2), beta = matrix(c(.7, 1.1, 1.7, 2), 2))
dt <- c(0, .2, .3, .4, .1, .5)
types <- c(0, 1, 2, 1, 1, 2)
counts <- cbind(2 + cumsum(types == 1), 3 + cumsum(types == 2))
pre <- t(vapply(seq_along(dt), function(n) .2 + .01 * counts[max(1L, n - 1L), ], numeric(2)))
post <- rbind(pre[1, ], .2 + .01 * counts[-1, , drop = FALSE])
phi <- pre * dt
expected <- sum(log(pre[cbind(2:length(dt), types[-1])])) - sum(phi)
for (initial in list(0, NULL)) for (infer in c(FALSE, TRUE)) for (engine in c("r", "auto", "cpp")) {
  mu_calls <- integer(); native_calls <- 0L
  value <- evaluate(engine, h, dt, infer, initial, type = types, N0 = c(2, 3))
  close(as.numeric(if (infer) value$loglikelihood else value), expected)
  stopifnot(identical(as.integer(mu_calls), c(1L, if (is.null(initial)) 1L, seq.int(2L, length(dt) + 1L))))
  stopifnot(native_calls == if (engine == "r") 0L else 1L)
  if (infer) {
    close(unname(value$lambda), unname(pre))
    close(unname(value$rambda), unname(post))
    stopifnot(all(value$lambda_component == 0), all(value$rambda_component == 0))
    for (k in 1:2) close(value[[paste0("res_process", k)]], diff(cumsum(phi[, k])[which(types == k)]))
  }
}
# State feedback, marks, maps and named/attributed kernels. Same-name defaults
# in mu, beta and impact intentionally have different values.
records <- list()
mu <- function(param = c(q = .2), n, Nc, ...) {
  previous <- Nc[max(1L, n - 1L), ]
  matrix(exp(param) * c(.5, .6) + .01 * previous, 2)
}
beta <- function(param = c(q = .9)) {
  x <- matrix(param + c(.7, 1.1, 2.3, .9, 3.1, .9), 2,
    dimnames = list(c("up", "down"), c("A", "B", "C")))
  attr(x, "label") <- "decay"
  x
}
impact <- function(param = c(q = .03), n, N, Nc, lambda_component_n, lambda, lambda_component, mu, alpha, beta, ...) {
  stopifnot(is.function(mu), identical(dimnames(lambda_component_n), dimnames(beta)),
    identical(attr(lambda_component_n, "label"), "decay"),
    identical(colnames(lambda), c("lambda1", "lambda2")),
    identical(colnames(lambda_component), c("lambda11", "lambda12", "lambda13", "lambda21", "lambda22", "lambda23")))
  close(as.vector(matrix(lambda_component[n, ], 2, byrow = TRUE)), as.vector(lambda_component_n), 0)
  close(as.numeric(lambda[n, ]), as.vector(mu(n = n, Nc = Nc)) + as.vector(rowSums(lambda_component_n)), 0)
  if (n < length(N[, 1])) stopifnot(all(lambda[(n + 1L):nrow(lambda), ] == 0),
    all(lambda_component[(n + 1L):nrow(lambda_component), ] == 0))
  records[[length(records) + 1L]] <<- list(n = n, current = lambda_component_n,
    lambda = lambda, components = lambda_component, param = param)
  value <- .02 * param + .001 * (Nc[n, 1] - N[n - 1L, 2]) - .03 * lambda_component_n +
    .001 * lambda_component_n["up", "A"]
  lambda_component_n[1, 1] <- 999
  lambda_component[n, 1] <- 999
  lambda[n, 1] <- 999
  alpha[1, 1] <- 999
  beta[1, 1] <- 999
  value
}
h <- new("hspec", mu = mu, alpha = matrix(c(.08, -.02, .12, .06, .03, .09), 2),
  beta = beta, eta = matrix(.01, 2, 3), impact = impact,
  type_col_map = list(c(3, 1), c(3, 2)))
dt <- c(0, .3, 0, .7, .2, 1.1, .05)
types <- c(0, 1, 2, 2, 1, 2, 1)
marks <- c(0, 1, 3, .5, 2, 1, 2.5)
initial <- matrix(c(.11, -.03, .07, .13, -.02, .08), 2)
for (infer in c(FALSE, TRUE)) {
  records <- list()
  reference <- evaluate("r", h, dt, infer, initial, type = types, mark = marks, N0 = c(2, 3), Nc0 = c(4, 5))
  expected_records <- records
  for (engine in c("auto", "cpp")) {
    records <- list(); native_calls <- 0L
    close(evaluate(engine, h, dt, infer, initial, type = types, mark = marks, N0 = c(2, 3), Nc0 = c(4, 5)), reference)
    stopifnot(identical(records, expected_records), native_calls == 1L)
  }
}
# The R block-kernel cache and a named numeric beta have the same semantics.
long_dt <- c(0, rep(dt[-1], length.out = 39L))
long_types <- c(0, rep(types[-1], length.out = 39L))
long_marks <- c(0, rep(marks[-1], length.out = 39L))
named_numeric <- h
named_numeric@beta <- emhawkes:::evalf(h@beta)
model_before <- serialize(named_numeric, NULL)
reference <- evaluate("r", named_numeric, long_dt, TRUE, initial,
  type = long_types, mark = long_marks, N0 = c(2, 3), Nc0 = c(4, 5))
close(evaluate("cpp", named_numeric, long_dt, TRUE, initial,
  type = long_types, mark = long_marks, N0 = c(2, 3), Nc0 = c(4, 5)), reference)
stopifnot(identical(model_before, serialize(named_numeric, NULL)))
# Lazy captures through dots receive full histories in inference and NULL
# histories without inference, matching R when only current state is declared.
saved <- list()
lazy <- new("hspec", mu = .3, alpha = .1, beta = 1.2,
  impact = function(n, lambda_component_n, ...) {
    saved[[length(saved) + 1L]] <<- function() list(n = n, current = lambda_component_n, dots = list(...))
    invisible(gc(FALSE))
    matrix(.01)
  })
for (infer in c(FALSE, TRUE)) {
  saved <- list(); reference <- evaluate("r", lazy, dt, infer)
  expected_snapshots <- lapply(saved, function(f) f())
  saved <- list(); actual <- evaluate("cpp", lazy, dt, infer)
  close(actual, reference)
  stopifnot(identical(lapply(saved, function(f) f()), expected_snapshots))
}
# Terminal callback errors are preserved even if only the likelihood is wanted.
terminal <- new("hspec", mu = function(n, inter_arrival, ...) {
  if (n > length(inter_arrival)) stop("terminal baseline sentinel")
  matrix(.3)
}, alpha = .1, beta = 1.2)
for (engine in c("r", "auto", "cpp")) fails(evaluate(engine, terminal, dt, FALSE), "terminal baseline sentinel")
terminal@mu <- function(n, inter_arrival, ...) {
  if (n > length(inter_arrival)) return(NULL)
  matrix(.3)
}
close(evaluate("cpp", terminal, dt, FALSE), evaluate("r", terminal, dt, FALSE))
empty <- terminal
empty@mu <- function(n, ...) { if (n != 1L) stop("no terminal call on empty path"); matrix(.3) }
empty@impact <- function(lambda_component_n) stop("no impact on empty path")
close(evaluate("cpp", empty, 0), evaluate("r", empty, 0))
# Unavailable baseline state reads remain on R; ordinary callback errors do not
# become fallback. Nonfinite unselected impact entries preserve finite R results.
latent <- lazy
latent@mu <- function(n, ...) matrix(.3 + .01 * sum(list(...)$lambda_component_n))
close(evaluate("auto", latent, dt), evaluate("r", latent, dt))
fails(evaluate("cpp", latent, dt), "latent intensity")
unused <- new("hspec", mu = c(.3, .4), alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2),
  impact = function(n, type, lambda_component_n) { x <- matrix(NA_real_, 2, 2); x[, type[n]] <- .01; x })
for (infer in c(FALSE, TRUE)) close(evaluate("auto", unused, dt, infer, type = types),
                                  evaluate("r", unused, dt, infer, type = types))
fails(evaluate("cpp", unused, dt, type = types), "finite numeric")
failures <- 0L
failed <- lazy
failed@impact <- function(lambda_component_n) { failures <<- failures + 1L; stop("state replay sentinel") }
fails(evaluate("auto", failed, dt), "state replay sentinel")
stopifnot(failures == 1L)
# Native baseline dimensions must protect the terminal post-event lookup.
a <- list(dt = c(0, .2, .3), type = c(0L, 1L, 1L), mark = c(0, 1, 1), mu = .3,
  alpha = matrix(0), beta = matrix(1), eta = NULL, initial = matrix(0), mapping = list(0L), infer = TRUE)
fails(do.call(native, c(a, list(baseline = matrix(.3, 1, 3)))), "terminal")
value <- do.call(native, c(a, list(baseline = matrix(c(.3, .3, .4, .5), 1))))
close(as.vector(value$lambda), c(.3, .3, .4))
close(as.vector(value$rambda), c(.3, .4, .5))
assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
options(old)
compiler::enableJIT(old_jit)
