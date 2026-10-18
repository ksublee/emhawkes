library(emhawkes)
close <- function(actual, expected) stopifnot(
  isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)))
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

# Independent compound Poisson processes have variance
# T * (mu[1] * E[mark1^2] + mu[2] * E[mark2^2]).
poisson <- new("hspec", mu = c(.2, .3), alpha = matrix(0, 2, 2),
  beta = matrix(c(1, 2), 2, 2))
dt <- c(0, .2, .3, .4, .5)
types <- c(0, 1, 2, 1, 2)
old <- options(emhawkes.engine = "auto", warn = 2)
for (engine in c("r", "cpp")) for (dep in c(FALSE, TRUE)) for (horizon in c(0, .5, 3)) {
  options(emhawkes.engine = engine)
  value <- function(h = poisson, mark = NULL, type = types) hvol(h, horizon = horizon,
    inter_arrival = dt, type = type, mark = mark, dependence = dep, lambda_component0 = 0)
  close(value(), sqrt(.5 * horizon))
  close(value(mark = c(0, 1, 2, 3, 4)), sqrt((.2 * 5 + .3 * 10) * horizon))
  explicit <- poisson
  explicit@eta <- matrix(0, 2, 2)
  close(value(explicit, c(0, 1, 1, 1, 1)), value())
  # The initial row is not an event, even if its type/mark are nonzero.
  close(value(type = c(1, types[-1]), mark = c(100, 1, 1, 1, 1)), value())
}

# Unit-mark Hawkes counts have long-run covariance
# (I-G)^-1 diag(mean intensity) (I-G)^-T, independent of this moment recursion.
for (i in seq_len(5)) {
  set.seed(300 + i)
  h <- new("hspec", mu = runif(2, .2, .8), alpha = matrix(runif(4, .01, .1), 2),
    beta = matrix(runif(2, .8, 2), 2, 2))
  resolvent <- solve(diag(2) - h@alpha / h@beta)
  means <- resolvent %*% h@mu
  covariance <- resolvent %*% diag(as.numeric(means)) %*% t(resolvent)
  reference <- sqrt(c(1, -1) %*% covariance %*% c(1, -1))
  for (dep in c(FALSE, TRUE)) close(hvol(h, inter_arrival = dt, type = types,
    dependence = dep, lambda_component0 = 0), reference)
}

args <- list(inter_arrival = dt, type = types, lambda_component0 = 0)
check <- function(h = poisson, ...) do.call(hvol, c(list(object = h), args, list(...)))
for (horizon in list(-1, NA_real_, Inf, c(1, 2))) fails(check(horizon = horizon), "horizon")
fails(check(dependence = NA), "dependence")
fails(hvol(poisson, inter_arrival = c(0, .2), type = c(0, 1)), "each type")
unequal <- poisson
unequal@beta[1, 2] <- 3
fails(check(unequal), "equal decay")
mapped <- poisson
mapped@type_col_map <- list(2L, 1L)
fails(check(mapped), "one kernel")
mapped@type_col_map <- list(1L, 2L)
close(check(mapped), check())
callback <- poisson
callback@impact <- function(...) stop("unsupported callback must not run")
fails(check(callback), "custom impact")
residual <- set_residual(poisson, c(rate = 1), dresidual = dexp, presidual = pexp)
fails(check(residual), "residual distribution")
univariate <- new("hspec", mu = .2, alpha = .1, beta = 1)
fails(check(univariate), "bivariate")
for (a in c(1, 1.1, 2)) for (dep in c(FALSE, TRUE)) {
  unstable <- poisson
  unstable@alpha <- diag(a, 2)
  fails(check(unstable, dependence = dep), "stationary mean-intensity")
}
marked_unstable <- poisson
marked_unstable@alpha <- matrix(.1, 2, 2)
marked_unstable@eta <- matrix(.3, 2, 2)
fails(hvol(marked_unstable, inter_arrival = dt, type = types,
           mark = c(0, 3, 3, 3, 3), lambda_component0 = 0), "stationary mean-intensity")
options(old)
