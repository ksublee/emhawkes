library(emhawkes)

close <- function(x, y, tol = 1e-9) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

# Splitting an exponential kernel leaves its total intensity unchanged.
h <- new("hspec", mu = .2, alpha = .3, beta = 2)
multi <- new("hspec", mu = .2, alpha = matrix(c(.1, .2), 1),
             beta = matrix(2, 1, 2), type_col_map = list(1:2))
dt <- c(0, .2, .7, .3)
single <- logLik(h, dt, lambda_component0 = matrix(.15), infer = TRUE)
split <- logLik(multi, dt, lambda_component0 = matrix(c(.05, .1), 1), infer = TRUE)
close(single$loglikelihood, split$loglikelihood)
close(single$lambda, split$lambda)
close(single$res_process1, split$res_process1)
explicit <- logLik(h, dt, type = c(0, 1, 1, 1), lambda_component0 = .15, infer = TRUE)
close(single$res_process1, explicit$res_process1)

# Kernel column layouts cannot change event types in a Poisson special case.
h <- new("hspec", mu = c(100, 1), alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2))
set.seed(45)
reference <- hsim(h, size = 300, lambda_component0 = matrix(0, 2, 2))
for (mapping in list(list(1:2, 3:4), list(c(2, 4), c(1, 3)))) {
  multi <- new("hspec", mu = c(100, 1), alpha = matrix(0, 2, 4),
               beta = matrix(1, 2, 4), type_col_map = mapping)
  set.seed(45)
  path <- hsim(multi, size = 300, lambda_component0 = matrix(0, 2, 4))
  close(path$type, reference$type)
  close(path$inter_arrival, reference$inter_arrival)
}

# The baseline jumps immediately after an observed event.
mu <- function(param = c(base = .1, increment = .1), n, N, ...) {
  counts <- if (n <= 1) c(0, 0) else N[n - 1, ]
  matrix(param[1] + param[2] * counts, 2, 1)
}
state <- new("hspec", mu = mu, alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2))
types <- c(0, 1, 2, 1, 2)
out <- logLik(state, c(0, rep(.1, 4)), type = types,
              lambda_component0 = matrix(0, 2, 2), infer = TRUE)
counts <- cbind(c(0, 1, 1, 2, 2), c(0, 0, 1, 1, 2))
close(out$rambda, .1 + .1 * counts)
close(out$lambda[-1, ], .1 + .1 * counts[-nrow(counts), ])
set.seed(4)
path <- hsim(state, size = 30, lambda_component0 = matrix(0, 2, 2))
restored <- logLik(state, path$inter_arrival, type = path$type,
                   lambda_component0 = matrix(0, 2, 2), infer = TRUE)
close(path$lambda, restored$lambda)
close(path$rambda, restored$rambda)

# Empty and singleton event types have no complete inter-event residual.
h <- new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 2), beta = matrix(1, 2, 2))
for (types in list(c(0, 1, 1), c(0, 1, 2), 0)) {
  dt <- c(0, rep(.5, length(types) - 1))
  out <- logLik(h, dt, type = types, lambda_component0 = matrix(.1, 2, 2), infer = TRUE)
  close(out$loglikelihood, logLik(h, dt, type = types, lambda_component0 = matrix(.1, 2, 2)))
  stopifnot(length(out$res_process2) == 0)
  legacy <- suppressWarnings(residual_process(2, dt, types, out$rambda_component,
    h@mu, h@beta))
  stopifnot(identical(legacy, numeric(0)))
}

# The initial row is not an event, even with a nonzero initial type.
out <- logLik(h, c(0, .2, .3), type = c(2, 2, 2), lambda_component0 = .1, infer = TRUE)
close(suppressWarnings(residual_process(2, c(0, .2, .3), c(2, 2, 2),
  out$rambda_component, h@mu, h@beta)), out$res_process2)

# A zero-baseline process can have no further events. No infinite row is added.
extinct <- new("hspec", mu = 0, alpha = matrix(0, 1, 2),
               beta = matrix(c(1, 2), 1), type_col_map = list(1:2))
set.seed(42)
path <- suppressWarnings(hsim(extinct, size = 5, lambda_component0 = matrix(.1, 1, 2)))
stopifnot(length(path$type) == 1, all(is.finite(as.matrix(path))))
stopifnot(nrow(path$lambda) == 1, nrow(path$N) == 1)
close(logLik(extinct, path$inter_arrival, type = path$type,
             lambda_component0 = matrix(.1, 1, 2)), 0)
invert <- getFromNamespace("invert_exp_hazard", "emhawkes")
waiting <- invert(.1, c(.1, .1), 0, c(1, 2))
close(integrate(function(t) .1 * exp(-t) + .1 * exp(-2 * t), 0, waiting)$value, .1)
stopifnot(is.infinite(invert(.2, c(.1, .1), 0, c(1, 2))))
close(invert(1, 0, 1e-6, 1), 1e6)

# The small-decay limit has constant intensity and integrated intensity one.
h <- new("hspec", mu = 0, alpha = 0, beta = 1e-18)
for (infer in c(FALSE, TRUE)) {
  out <- logLik(h, c(0, 1), type = c(0, 1), lambda_component0 = matrix(1), infer = infer)
  close(if (infer) out$loglikelihood else out, -1)
}

# Exponential residuals recover ordinary Hawkes even far into the tail.
h <- new("hspec", mu = c(1, 1), alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2))
exp_h <- set_residual(h, param = c(rate = 1), dresidual = dexp,
                     presidual = pexp, qresidual = qexp, rresidual = rexp)
dt <- c(0, 50, 1, 800, 1); types <- c(0, 1, 2, 1, 2)
ordinary <- logLik(h, dt, type = types, lambda_component0 = matrix(0, 2, 2), infer = TRUE)
flexible <- logLik(exp_h, dt, type = types, lambda_component0 = matrix(0, 2, 2), infer = TRUE)
close(ordinary$loglikelihood, flexible$loglikelihood)
for (k in 1:2) {
  close(ordinary[[paste0("res_process", k)]], flexible[[paste0("res_process", k)]])
  close(ordinary[[paste0("res_process", k)]], flexible[[paste0("eres_process", k)]])
}

# Old callbacks without log/tail options still work; parameters remain bindable.
legacy <- set_residual(h, param = c(rate = 2),
  dresidual = function(x, rate) dexp(x, rate),
  presidual = function(q, rate) pexp(q, rate),
  qresidual = function(p, rate) qexp(p, rate))
close(legacy@dresidual(.3), dexp(.3, 2))
close(legacy@presidual(.3), pexp(.3, 2))
close(legacy@qresidual(.3), qexp(.3, 2))
close(exp_h@dresidual(.3, param = c(rate = 2)), dexp(.3, 2))

# Constructor validation and multi-kernel branching matrix.
fails(new("hspec", mu = .2, alpha = .1, beta = -1))
fails(new("hspec", mu = -.2, alpha = .1, beta = 1))
fails(new("hspec", mu = .2, alpha = .1, beta = Inf))
fails(new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 4),
          beta = matrix(1, 2, 4), type_col_map = list(c(1, 9), c(2, 4))))
stable <- new("hspec", mu = .2, alpha = .1, beta = 1, stability_check = TRUE)
stable <- new("hspec", mu = .2, alpha = function(param = c(a = .1)) matrix(param),
              beta = 1, stability_check = TRUE)
stable <- new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 4),
              beta = matrix(1, 2, 4), type_col_map = list(c(1, 3), c(2, 4)),
              stability_check = TRUE)
warning_seen <- FALSE
withCallingHandlers(new("hspec", mu = .2, alpha = 1, beta = 1, stability_check = TRUE),
  warning = function(w) { warning_seen <<- TRUE; invokeRestart("muffleWarning") })
stopifnot(warning_seen)

# The deprecated initial-condition alias reaches the same calculation.
h <- new("hspec", mu = .2, alpha = .1, beta = 1)
old <- suppressWarnings(infer_lambda(h, c(0, .2, .3), lambda0 = matrix(.1)))
current <- logLik(h, c(0, .2, .3), lambda_component0 = matrix(.1), infer = TRUE)
close(old$lambda, current$lambda)
