library(emhawkes)
close <- function(actual, expected, tolerance = 1e-8) stopifnot(
  isTRUE(all.equal(actual, expected, tolerance = tolerance)))
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
h <- new("hspec", mu = 0, alpha = 0, beta = 1)
for (excitation in c(0, .1, 1, 50)) stopifnot(identical(expected_tau(h, excitation), Inf))
for (mu in c(1e-8, 1e-6, .5, 1e4)) close(expected_tau(h, 0, mu = mu), 1 / mu)

# For ordinary Hawkes residuals, survival is a Poisson mixture of exponential
# survivals. Summing its independent closed-form moments checks two time scales.
for (mu in c(1e-6, .2, 2, 1e4)) for (beta in c(.1, 1, 1e3)) for (a in c(.1, 1, 20)) {
  counts <- 0:200
  reference <- sum(dpois(counts, a) / (mu + beta * counts))
  close(expected_tau(h, a * beta, mu = mu, beta = beta, tol = 1e-9), reference)
}
for (mu in c(1e-300, 1e-320)) for (a in c(1, 100)) {
  reference <- sum(dpois(0:250, a) / (mu + 0:250))
  close(expected_tau(h, a, mu = mu, tol = 1e-9), reference)
}
reference <- sum(dpois(0:200, 1) / (.2 + 0:200))
for (unit in c(1e-100, 1e-20, 1, 1e20, 1e100)) close(
  expected_tau(h, 1/unit, mu = .2/unit, beta = 1/unit, tol = 1e-9) / unit, reference)
close(expected_tau(h, 0, mu = .5, max_upper = 3),
  2 * (1 - exp(-1.5) * (1 + 1.5)))
close(expected_tau(h, 1, max_upper = 0), 0)

exponential <- set_residual(h, c(rate = 2), dresidual = dexp, presidual = pexp)
stopifnot(identical(expected_tau(exponential, 1), Inf))
close(expected_tau(exponential, 0, mu = .5), 1)
close(expected_tau(exponential, 0, mu = 1e-6, tol = 1e-9), 5e5)
no_cdf <- set_residual(h, c(rate = 2), dresidual = dexp)
fails(expected_tau(no_cdf, 1), "presidual")
# A bounded residual law can guarantee an event despite a zero baseline.
bounded <- set_residual(h, c(min = 0, max = .5), dresidual = dunif, presidual = punif)
close(expected_tau(bounded, 1, tol = 1e-9), 1 - log(2))
for (excitation in list(c(1, 2), -1, NA_real_)) fails(expected_tau(h, excitation), "scalars")
fails(expected_tau(h, 0, type = 2), "type")
unequal <- new("hspec", mu = c(.2, .3), alpha = matrix(.01, 2, 2),
  beta = matrix(c(1, 2, 3, 2), 2))
fails(expected_tau(unequal, .1), "equal decay")

# Independent exponential/Gamma first moments check integration when the time
# scale is far from one. A finite limit returns E[T * I(T <= upper)].
for (mu in c(1e-12, 1e-6, 1, 1e4, 1e6, 1e12)) {
  upper <- 5 / mu
  close(expected_tau(h, 0, mu = mu, max_upper = upper, tol = 1e-9) * mu,
        pgamma(5, shape = 2))
  close(expected_tau(no_cdf, 0, mu = mu, tol = 1e-9) * mu, .5)
  close(expected_tau(no_cdf, 0, mu = mu, max_upper = upper, tol = 1e-9) * mu,
        .5 * pgamma(10, shape = 2))
}
gamma_density <- set_residual(h, c(shape = 2, rate = 3), dresidual = dgamma)
for (mu in c(1e-8, 1, 1e8)) close(
  expected_tau(gamma_density, 0, mu = mu, tol = 1e-9) * mu, 2 / 3)
close(expected_tau(h, 0, mu = 1e4, max_upper = 1, tol = 1e-9) * 1e4, 1)
