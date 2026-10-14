# Shared fixtures for timing and previous-version comparisons.
observed_impact_cases <- function() {
  zero <- function(param = numeric(0)) matrix(0)
  impact_count <- function(param = c(a = .3, extra = .2), n, N, ...) {
    matrix(param["a"] + param["extra"] * (N[n, 1] %% 2))
  }
  jac_count <- function(param, n, N, ...) {
    out <- matrix(0, 1, length(param), dimnames = list(NULL, names(param)))
    out[, "a"] <- 1
    out[, "extra"] <- N[n, 1] %% 2
    out
  }
  impact_mark <- function(param = c(self = .3, cross = .1, gain = .07), n, N, Nc, mark, ...) {
    weight <- 1 + .2 * (N[n, 1] %% 2) + .1 * Nc[n - 1, 2] / (1 + N[n - 1, 2])
    weight * matrix(param[c("self", "cross", "cross", "self")], 2) + param["gain"] * (mark[n] - 1)
  }
  jac_mark <- function(param, n, N, Nc, mark, ...) {
    weight <- 1 + .2 * (N[n, 1] %% 2) + .1 * Nc[n - 1, 2] / (1 + N[n - 1, 2])
    out <- matrix(0, 4, length(param), dimnames = list(NULL, names(param)))
    out[, "self"] <- weight * c(1, 0, 0, 1)
    out[, "cross"] <- weight * c(0, 1, 1, 0)
    out[, "gain"] <- mark[n] - 1
    out
  }
  beta <- function(param = c(b = 1.3)) matrix(c(.9, 1.5, 1.2, 2) * param["b"], 2)
  impact_multi <- function(param = c(slow = .25, fast = .9), n, N, Nc, ...) {
    matrix(c(param["slow"] * (1 + .1 * Nc[n - 1, 1] / (1 + N[n - 1, 1])),
             param["fast"] * (1 + .2 * (N[n, 1] %% 2))), 1)
  }
  jac_multi <- function(param, n, N, Nc, ...) {
    out <- matrix(0, 2, length(param), dimnames = list(NULL, names(param)))
    out[1, "slow"] <- 1 + .1 * Nc[n - 1, 1] / (1 + N[n - 1, 1])
    out[2, "fast"] <- 1 + .2 * (N[n, 1] %% 2)
    out
  }
  unary <- new("hspec", mu = .35, alpha = zero, beta = 1.5, impact = impact_count)
  unary_start <- unary
  unary_start@mu[] <- .3
  unary_start@beta[] <- 1.2
  formals(unary_start@impact)$param <- c(a = .25, extra = .15)
  marked <- new("hspec", mu = c(.25, .35),
    alpha = function(param = numeric(0)) matrix(0, 2, 2), beta = beta,
    impact = impact_mark, rmark = function(n, ...) 1 + (n %% 3))
  marked_start <- marked
  marked_start@mu[] <- c(.2, .3)
  formals(marked_start@beta)$param <- c(b = 1)
  formals(marked_start@impact)$param <- c(self = .25, cross = .08, gain = .05)
  multi <- new("hspec", mu = .3, alpha = function(param = numeric(0)) matrix(0, 1, 2),
    beta = function(param = numeric(0)) matrix(c(.8, 5), 1),
    impact = impact_multi, type_col_map = list(1:2))
  multi_start <- multi
  multi_start@mu[] <- .25
  formals(multi_start@impact)$param <- c(slow = .2, fast = .75)
  list(univariate_counts = list(truth = unary, start = unary_start, seed = 9351L,
         initial = NULL, jacobian = list(impact = jac_count)),
       bivariate_counts_marks = list(truth = marked, start = marked_start, seed = 9352L,
         initial = .1, jacobian = list(beta = matrix(c(.9, 1.5, 1.2, 2), dimnames = list(NULL, "b")), impact = jac_mark)),
       multi_kernel_counts = list(truth = multi, start = multi_start, seed = 9353L,
         initial = .1, jacobian = list(impact = jac_multi)))
}
