library(emhawkes)

assert_close <- function(actual, expected, tolerance = 1e-9) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                            tolerance = tolerance)))
}

# Compute intensities from the complete event history, independently of the
# recursion in logLik(), and integrate them numerically between events.
reference_path <- function(object, inter_arrival, type, mark, initial) {
  arrival <- cumsum(inter_arrival)
  size <- length(arrival)
  dimens <- object@dimens
  beta <- object@beta
  mu <- as.vector(object@mu)

  intensity <- function(time, past) {
    component <- initial * exp(-beta * time)
    for (j in past) {
      columns <- if (is.null(object@type_col_map)) {
        type[j]
      } else {
        object@type_col_map[[type[j]]]
      }
      jump <- matrix(0, dimens, ncol(beta))
      jump[, columns] <- object@alpha[, columns]
      if (!is.null(object@eta)) {
        jump[, columns] <- jump[, columns] +
          object@eta[, columns] * (mark[j] - 1)
      }
      component <- component + jump * exp(-beta * (time - arrival[j]))
    }
    mu + rowSums(component)
  }

  lambda <- matrix(0, size, dimens)
  lambda[1, ] <- intensity(0, integer())
  integrated <- matrix(0, size, dimens)
  for (n in 2:size) {
    past <- if (n > 2) 2:(n - 1) else integer()
    lambda[n, ] <- intensity(arrival[n], past)
    for (k in seq_len(dimens)) {
      integrated[n, k] <- integrate(
        function(t) vapply(t, function(s) intensity(s, past)[k], numeric(1)),
        lower = arrival[n - 1], upper = arrival[n], rel.tol = 1e-11
      )$value
    }
  }

  event_log <- sum(log(lambda[cbind(2:size, type[-1])]))
  list(lambda = lambda, integrated = integrated, event_log = event_log,
       loglikelihood = event_log - sum(integrated))
}

check_inference <- function(object, inter_arrival, type, mark, initial) {
  reference <- reference_path(object, inter_arrival, type, mark, initial)
  args <- list(object = object, inter_arrival = inter_arrival, type = type,
               mark = mark, lambda_component0 = initial)
  likelihood <- do.call(logLik, args)
  inferred <- do.call(logLik, c(args, list(infer = TRUE)))

  assert_close(likelihood, reference$loglikelihood)
  assert_close(inferred$loglikelihood, reference$loglikelihood)
  assert_close(inferred$lambda, reference$lambda)

  for (k in seq_len(object@dimens)) {
    events <- which(type == k)
    residual <- vapply(seq_len(length(events) - 1), function(j) {
      sum(reference$integrated[(events[j] + 1):events[j + 1], k])
    }, numeric(1))
    assert_close(inferred[[paste0("res_process", k)]], residual)
  }

  # Exp(1) residuals must recover the ordinary Hawkes likelihood and residuals.
  exponential <- set_residual(object, param = c(rate = 1),
                             dresidual = dexp, presidual = pexp,
                             qresidual = qexp, rresidual = rexp)
  args$object <- exponential
  exponential_inferred <- do.call(logLik, c(args, list(infer = TRUE)))
  assert_close(do.call(logLik, args), reference$loglikelihood)
  assert_close(exponential_inferred$loglikelihood, reference$loglikelihood)
  for (k in seq_len(object@dimens)) {
    assert_close(exponential_inferred[[paste0("res_process", k)]],
                 inferred[[paste0("res_process", k)]])
    assert_close(exponential_inferred[[paste0("eres_process", k)]],
                 inferred[[paste0("res_process", k)]])
  }

  # Check a non-exponential residual law against the numerical integrals.
  gamma <- set_residual(object, param = c(shape = 2, rate = 2),
                       dresidual = dgamma, presidual = pgamma,
                       qresidual = qgamma, rresidual = rgamma)
  args$object <- gamma
  gamma_loglikelihood <- reference$event_log
  for (n in 2:length(type)) {
    for (k in seq_len(object@dimens)) {
      value <- reference$integrated[n, k]
      gamma_loglikelihood <- gamma_loglikelihood + if (type[n] == k) {
        dgamma(value, shape = 2, rate = 2, log = TRUE)
      } else {
        pgamma(value, shape = 2, rate = 2, lower.tail = FALSE, log.p = TRUE)
      }
    }
  }
  assert_close(do.call(logLik, args), gamma_loglikelihood)
  assert_close(do.call(logLik, c(args, list(infer = TRUE)))$loglikelihood,
               gamma_loglikelihood)
}

inter_arrival <- c(0, 0.4, 0.9, 0.2, 0.7, 1.1, 0.3)
type <- c(0, 1, 2, 1, 2, 1, 2)
mark <- c(0, 1, 2, 1, 3, 2, 1)
initial <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
alpha <- matrix(c(0.3, 0.1, 0.2, 0.3), nrow = 2)
beta <- matrix(c(1, 1.2, 1.5, 1.8), nrow = 2)

unequal <- new("hspec", mu = c(0.2, 0.3), alpha = alpha, beta = beta)
check_inference(unequal, inter_arrival, type, mark, initial)

# Each row may have its own common rate without requiring all rows to match.
equal_rows <- unequal
equal_rows@beta <- matrix(c(1, 1.2, 1, 1.2), nrow = 2)
check_inference(equal_rows, inter_arrival, type, mark, initial)

reverse <- unequal
reverse@beta <- beta[, 2:1]
check_inference(reverse, inter_arrival, type, mark, initial)

marked <- unequal
marked@eta <- matrix(0.05, nrow = 2, ncol = 2)
check_inference(marked, inter_arrival, type, mark, initial)

multi <- new("hspec", mu = c(0.2, 0.3),
             alpha = cbind(alpha, alpha / 3), beta = cbind(beta, beta * 3),
             type_col_map = list(c(1, 3), c(2, 4)))
check_inference(multi, inter_arrival, type, mark, cbind(initial, initial / 2))

univariate <- new("hspec", mu = 0.2, alpha = 0.3, beta = 1)
check_inference(univariate, inter_arrival, c(0, rep(1, 6)), mark, matrix(0.1))
