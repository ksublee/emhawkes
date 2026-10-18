library(emhawkes)
source(system.file("examples", "fosset-spread.R", package = "emhawkes"))
close_fosset <- function(x, y, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
}

# Direct sums over past widening events implement paper equations (8a)--(8b).
# This oracle has no package recursions, callback helpers or kernel integrators.
fosset_reference <- function(theta, dt, type, spread0 = 1) {
  p <- fosset_parameters(theta); time <- cumsum(dt); size <- length(dt)
  count_plus <- cumsum(type == 1); count_minus <- cumsum(type == 2)
  spread <- spread0 + count_plus - count_minus
  pre_spread <- c(spread0, head(spread, -1))
  opening <- which(type == 1)
  pre <- post <- matrix(0, size, 2)
  for (i in seq_len(size)) {
    past <- opening[opening < i]
    excitation <- sum(p["a"] * p["b"] * exp(-p["b"] * (time[i] - time[past])))
    pre[i, ] <- c(p["nu_plus"] + excitation, p["nu_minus"] * (pre_spread[i] >= 2))
    post[i, ] <- c(pre[i, 1] + p["a"] * p["b"] * (type[i] == 1),
      p["nu_minus"] * (spread[i] >= 2))
  }
  horizon <- tail(time, 1)
  exposure <- sum(dt[-1] * (head(spread, -1) >= 2))
  integral_plus <- p["nu_plus"] * horizon +
    sum(p["a"] * (1 - exp(-p["b"] * (horizon - time[opening]))))
  loglik <- sum(log(pre[cbind(seq.int(2, size), type[-1])])) -
    integral_plus - p["nu_minus"] * exposure
  list(lambda = pre, rambda = post, loglik = unname(loglik), exposure = exposure)
}

old_fosset_options <- options(emhawkes.engine = "auto")
h <- fosset_spec()
theta <- c(log_nu_plus = 0, log_nu_minus = log(3), logit_a = qlogis(.4), log_b = log(2))
dt <- c(0, .2, .1, .4, .5, .15, .3, .2)
type <- c(0, 1, 2, 1, 1, 2, 2, 1)
for (spread0 in c(1, 3)) for (a in c(.2, .4, .7)) for (b in c(.5, 2, 3)) {
  candidate <- theta; candidate[c("logit_a", "log_b")] <- c(qlogis(a), log(b))
  object <- fosset_spec(a = a, b = b, spread0 = spread0)
  expected <- fosset_reference(candidate, dt, type, spread0)
  for (engine in c("r", "cpp")) {
    options(emhawkes.engine = engine)
    actual <- logLik(object, dt, type, lambda_component0 = 0, infer = TRUE)
    close_fosset(actual$lambda, expected$lambda)
    close_fosset(actual$rambda, expected$rambda)
    close_fosset(actual$loglikelihood, expected$loglik)
    close_fosset(logLik(object, dt, type, lambda_component0 = 0), expected$loglik)
  }
}
# A narrowing event at one tick has zero probability, without clipping its rate.
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  stopifnot(is.infinite(logLik(h, c(0, .1), c(0, 2), lambda_component0 = 0)),
    logLik(h, c(0, .1), c(0, 2), lambda_component0 = 0) < 0)
}

# R simulation respects the boundary and replays identically through both engines.
set.seed(20261018)
path <- hsim(h, size = 401, lambda_component0 = 0); seed <- .Random.seed
spread <- 1 + path$N[, 1] - path$N[, 2]
pre_spread <- c(1, head(spread, -1))
stopifnot(all(spread >= 1), any(spread == 1),
  all(pre_spread[path$type == 2] >= 2),
  all(path$lambda[pre_spread == 1, 2] == 0),
  all(path$rambda[spread == 1, 2] == 0))
close_fosset(path$lambda[, 2], 3 * (pre_spread >= 2))
close_fosset(path$rambda[, 2], 3 * (spread >= 2))
reference <- fosset_reference(theta, path$inter_arrival, path$type)
close_fosset(path$lambda, reference$lambda)
close_fosset(path$rambda, reference$rambda)
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  replay <- logLik(h, path$inter_arrival, path$type, lambda_component0 = 0, infer = TRUE)
  close_fosset(replay$lambda, path$lambda); close_fosset(replay$rambda, path$rambda)
  close_fosset(replay$loglikelihood, reference$loglik)
  stopifnot(identical(seed, .Random.seed))
}

# With the widening parameters fixed, the narrowing MLE is known exactly.
# The negative log-likelihood curvature in log(nu_minus) is the event count.
options(emhawkes.engine = "auto")
fit <- suppressMessages(hfit(fosset_spec(nu_minus = 2), path$inter_arrival, path$type,
  lambda_component0 = 0, optimizer = "stats", method = "BFGS",
  fixed = c("log_nu_plus", "logit_a", "log_b"),
  control = list(maxit = 300, reltol = 1e-12)))
count <- sum(path$type == 2)
stopifnot(fit$code == 0, fit_diagnostics(fit)$status == "ok")
close_fosset(fosset_parameters(coef(fit))["nu_minus"], count / reference$exposure, 1e-5)
close_fosset(vcov(fit)["log_nu_minus", "log_nu_minus"], 1 / count, 1e-4)
initial <- fit_initial_state(fit)
close_fosset(initial$lambda_component0, matrix(0, 2, 2))
restored <- fitted_hspec(fit)
close_fosset(do.call(logLik, c(list(object = restored, inter_arrival = path$inter_arrival,
  type = path$type), initial)), logLik(fit))
options(old_fosset_options)
