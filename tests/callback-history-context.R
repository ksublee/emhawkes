library(emhawkes)
close <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-12)))
old <- options(emhawkes.engine = "r")
dt <- c(0, .2, .3, .4)
marks <- c(0, 2, 1, 3)
dots <- function(...) {
  x <- list(...)
  matrix(.1 * sum(x$lambda) + .01 * sum(x$lambda_component) +
           .001 * x$N[x$n, 1] + .002 * x$Nc[x$n, 1])
}
named <- function(n, lambda, lambda_component, N, Nc) {
  matrix(.1 * sum(lambda) + .01 * sum(lambda_component) + .001 * N[n, 1] + .002 * Nc[n, 1])
}
h <- new("hspec", mu = 1, alpha = 0, beta = 1, impact = dots)
explicit <- h
explicit@impact <- named
evaluate <- function(object, infer, engine) {
  options(emhawkes.engine = engine)
  out <- logLik(object, dt, mark = marks, lambda_component0 = 0,
                N0 = 2, Nc0 = 3, infer = infer)
  if (infer) out$loglikelihood else out
}
reference <- evaluate(explicit, FALSE, "r")
for (engine in c("r", "auto")) for (infer in c(FALSE, TRUE))
  close(evaluate(h, infer, engine), reference)

# A baseline also receives counts and histories through dots. Its first call
# sees NULL histories; later calls see the same state regardless of infer.
baseline_dots <- function(n, ...) {
  x <- list(...)
  matrix(1 + .01 * sum(x$lambda) + .001 * x$N[max(1, n - 1), 1])
}
baseline_named <- function(n, lambda, N) matrix(1 + .01 * sum(lambda) + .001 * N[max(1, n - 1), 1])
h <- new("hspec", mu = baseline_dots, alpha = 0, beta = 1)
explicit <- h
explicit@mu <- baseline_named
reference <- evaluate(explicit, FALSE, "r")
for (engine in c("r", "auto")) for (infer in c(FALSE, TRUE))
  close(evaluate(h, infer, engine), reference)

# Mark likelihoods must not change when inference output is requested either.
h <- new("hspec", mu = 1, alpha = .1, beta = 1,
  dmark = function(...) exp(-.01 * sum(list(...)$lambda)))
reference <- evaluate(h, TRUE, "r")
for (engine in c("r", "auto")) close(evaluate(h, FALSE, engine), reference)

# Numerical fitting uses the sequential native adapter when current state is
# named. Its dots must carry the same full history as standalone likelihoods.
h <- new("hspec", mu = 1, alpha = 0, beta = 1,
  impact = function(lambda_component_n, ...) matrix(.1 * sum(list(...)$lambda) +
                                                   .01 * sum(lambda_component_n)))
reference <- evaluate(h, TRUE, "r")
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  fit <- suppressMessages(hfit(h, dt, mark = marks, N0 = 2, Nc0 = 3,
    lambda_component0 = 0, optimizer = "stats", final_hessian = "optimizer",
    fixed = rep(TRUE, 3)))
  close(logLik(fit), reference)
}
options(old)
