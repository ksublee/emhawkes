# Exact derivatives for the current-component Bid-ask example.
bid_ask_jacobians <- function(initial_ask_price, initial_bid_price, initial_level,
    initial_mid_price = (initial_ask_price + initial_bid_price) / 2) {
  force(initial_ask_price); force(initial_bid_price); force(initial_level); force(initial_mid_price)
  market_state <- function(n, Nc, after = FALSE) {
    if (n == 1L) return(list(level = initial_level, mid = initial_mid_price))
    k <- if (after) n else n - 1L
    ask <- initial_ask_price + Nc[k, 1] - Nc[k, 2]
    bid <- initial_bid_price + Nc[k, 3] - Nc[k, 4]
    list(level = initial_level + Nc[k, 1] - Nc[k, 2] - (Nc[k, 3] - Nc[k, 4]),
      mid = (ask + bid) / 2)
  }
  mu_jac <- function(param, n, Nc, ...) {
    s <- market_state(n, Nc)
    J <- matrix(0, 4, length(param), dimnames = list(NULL, names(param)))
    J[c(1, 4), 'mu1'] <- 1
    if (s$level > 0) J[c(2, 3), 'zeta1'] <- s$level / s$mid
    J
  }
  alpha_jac <- matrix(0, 16, 5, dimnames = list(NULL,
    c('alpha_s1', 'alpha_m', 'alpha_s2', 'alpha_w1', 'alpha_w2')))
  alpha_jac[c(1, 16), 'alpha_s1'] <- 1
  alpha_jac[c(5, 12), 'alpha_m'] <- 1
  alpha_jac[c(9, 8), 'alpha_s2'] <- 1
  alpha_jac[c(2, 15), 'alpha_w1'] <- 1
  alpha_jac[c(3, 14), 'alpha_w2'] <- 1
  impact_param_jac <- function(param, n, Nc, ...) {
    s <- market_state(n, Nc, after = TRUE)
    J <- matrix(0, 16, length(param), dimnames = list(NULL, names(param)))
    J[c(6, 7, 10, 11), 'xi'] <- s$level / s$mid
    J
  }
  impact_state_jac <- matrix(0, 16, 16)
  impact_state_jac[c(6, 10), c(2, 6, 10, 14)] <- -1
  impact_state_jac[c(7, 11), c(3, 7, 11, 15)] <- -1
  list(mu = mu_jac, alpha = alpha_jac,
    impact = list(param = impact_param_jac, state = impact_state_jac))
}
