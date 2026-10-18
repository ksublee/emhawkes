# Example helpers for Fosset, Bouchaud and Benzaquen (2020), Section 4.1.
# These local example functions are not exported package functions.

## ---- fosset-specification
fosset_spec <- function(nu_plus = 1, nu_minus = 3, a = 0.4, b = 2,
                        spread0 = 1) {
  values <- c(nu_plus, nu_minus, a, b, spread0)
  stopifnot(length(values) == 5L, all(is.finite(values)),
    nu_plus > 0, nu_minus > 0, a > 0, a < 1, b > 0,
    spread0 >= 1, spread0 == floor(spread0))
  # Type 1 widens the spread; type 2 narrows it. Read the previous row.
  mu <- function(param = c(log_nu_plus = log(nu_plus),
                           log_nu_minus = log(nu_minus)), n, N) {
    previous <- max(1L, n - 1L)
    spread <- spread0 + N[previous, 1] - N[previous, 2]
    matrix(c(exp(param[1]), exp(param[2]) * (spread >= 2)), 2, 1)
  }
  # The paper normalizes its kernel as a * b * exp(-b * t).
  # Functional components keep the other matrix entries fixed at zero.
  alpha <- function(param = c(logit_a = stats::qlogis(a), log_b = log(b))) {
    matrix(c(stats::plogis(param[1]) * exp(param[2]), 0, 0, 0), 2, 2)
  }
  beta <- function(param = c(log_b = log(b))) matrix(exp(param[1]), 2, 2)
  methods::new("hspec", mu = mu, alpha = alpha, beta = beta)
}

fosset_parameters <- function(theta) {
  c(nu_plus = unname(exp(theta["log_nu_plus"])),
    nu_minus = unname(exp(theta["log_nu_minus"])),
    a = unname(stats::plogis(theta["logit_a"])),
    b = unname(exp(theta["log_b"])))
}
