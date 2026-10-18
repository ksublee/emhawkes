# Fosset, Bouchaud and Benzaquen (2020), Section 4.1, equations (8a)--(8b).
# Local example helpers; these are not exported package functions.
# Event type 1 widens the spread; event type 2 narrows it.

fosset_spec <- function(nu_plus = 1, nu_minus = 3, a = 0.4, b = 2,
                        spread0 = 1) {
  values <- c(nu_plus, nu_minus, a, b, spread0)
  stopifnot(
    length(values) == 5L,
    all(is.finite(values)),
    nu_plus > 0,
    nu_minus > 0,
    a > 0, a < 1,
    b > 0,
    spread0 >= 1, spread0 == floor(spread0)
  )

  # 1. Baseline rates: narrowing is allowed only when the spread is open.
  baseline <- function(
    param = c(log_nu_plus = log(nu_plus), log_nu_minus = log(nu_minus)),
    n, N
  ) {
    widening_rate <- exp(param[["log_nu_plus"]])
    narrowing_rate <- exp(param[["log_nu_minus"]])

    # N stores cumulative counts. Row n - 1 is the state before event n.
    previous_event <- max(1L, n - 1L)
    widening_count <- N[previous_event, 1]
    narrowing_count <- N[previous_event, 2]
    spread_before_event <- spread0 + widening_count - narrowing_count
    narrowing_enabled <- spread_before_event >= 2

    matrix(
      c(widening_rate, narrowing_rate * narrowing_enabled),
      nrow = 2, ncol = 1
    )
  }

  # 2. Excitation jumps: only widening excites future widening.
  excitation <- function(
    param = c(logit_a = stats::qlogis(a), log_b = log(b))
  ) {
    branching_ratio <- stats::plogis(param[["logit_a"]])
    decay_rate <- exp(param[["log_b"]])

    # Rows are affected event types; columns are triggering event types.
    # The paper's kernel is a * b * exp(-b * t), so the jump is a * b.
    jumps <- matrix(0, nrow = 2, ncol = 2)
    jumps[1, 1] <- branching_ratio * decay_rate
    jumps
  }

  # 3. Decay rates: log_b is shared with the excitation function above.
  decay <- function(param = c(log_b = log(b))) {
    decay_rate <- exp(param[["log_b"]])
    matrix(decay_rate, nrow = 2, ncol = 2)
  }

  methods::new(
    "hspec",
    mu = baseline,
    alpha = excitation,
    beta = decay
  )
}

# Convert fitted log/logit coefficients back to the paper's parameters.
fosset_parameters <- function(theta) {
  c(
    nu_plus = exp(theta[["log_nu_plus"]]),
    nu_minus = exp(theta[["log_nu_minus"]]),
    a = stats::plogis(theta[["logit_a"]]),
    b = exp(theta[["log_b"]])
  )
}
