# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

close <- function(x, y, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance)))
}
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
fit_quiet <- function(...) suppressMessages(hfit(...))

# A genuine likelihood reproduces its fitted value and keeps maxLik methods.
h <- new("hspec", mu = .3, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .5, 1)
original <- serialize(h, NULL)
fit <- fit_quiet(h, dt, lambda_component0 = 0, N0 = 3, Nc0 = 7,
                 fixed = c(FALSE, TRUE, TRUE))
stopifnot(identical(class(fit), c("maxLik", "maxim")), fit$code == 0)
estimated <- fitted_hspec(fit)
initial <- fit_initial_state(fit)
close(initial, list(lambda_component0 = matrix(0), N0 = 3, Nc0 = 7))
close(estimated@mu, matrix(2), tolerance = 1e-5)
close(as.numeric(logLik(fit)), as.numeric(do.call(logLik,
  c(list(object = estimated, inter_arrival = dt), initial))))
stopifnot(is.matrix(vcov(fit)), inherits(summary(fit), "summary.maxLik"),
          identical(original, serialize(h, NULL)))
set.seed(28)
simulated <- do.call(hsim, c(list(object = estimated, size = 8), initial))
set.seed(28)
manual <- hsim(new("hspec", mu = coef(fit)["mu1"], alpha = 0, beta = 1),
               size = 8, lambda_component0 = 0, N0 = 3, Nc0 = 7)
close(as.matrix(simulated), as.matrix(manual))

# Omitted intensity remains a parameter-dependent initialization policy.
h <- new("hspec", mu = .3, alpha = .2, beta = 1)
fit <- fit_quiet(h, dt, fixed = c(FALSE, TRUE, TRUE))
stopifnot(identical(fit_initial_state(fit)$lambda_component0, "mean"))
close(as.numeric(logLik(fit)), as.numeric(suppressMessages(do.call(logLik,
  c(list(object = fitted_hspec(fit), inter_arrival = dt), fit_initial_state(fit))))))

# Non-square matrices, reordered kernel mappings, and both sharing settings.
alpha <- matrix(c(.1, .2, .1, .3, .2, .4, .3, .4), 2, byrow = TRUE)
h <- new("hspec", mu = c(.2, .3), alpha = alpha, beta = matrix(1.5, 2, 4),
         eta = matrix(.01, 2, 4), type_col_map = list(c(3, 1), c(4, 2)))
objective <- function(object) {
  -sum((object@mu - c(.5, .6))^2) - sum((object@alpha - .25)^2) -
    sum((object@beta - 2)^2) - sum((object@eta - .02)^2)
}

# An ordinary multi-kernel fit replays the same marked observation likelihood.
dt_multi <- c(0, .2, .3, .1)
type_multi <- c(0, 1, 2, 1)
mark_multi <- c(0, 2, 1, 3)
fit <- fit_quiet(h, dt_multi, type_multi, mark_multi,
  lambda_component0 = .1, N0 = c(4, 5), Nc0 = c(9, 10), iterlim = 0)
close(as.numeric(logLik(fit)), as.numeric(do.call(logLik,
  c(list(object = fitted_hspec(fit), inter_arrival = dt_multi,
         type = type_multi, mark = mark_multi), fit_initial_state(fit)))))
for (reduced in c(TRUE, FALSE)) {
  original <- serialize(h, NULL)
  fit <- fit_quiet(h, mylogLik = objective, reduced = reduced, lambda_component0 = .1)
  metadata <- attr(fit, "emhawkes")
  estimated <- fitted_hspec(fit)
  close(estimated@alpha, matrix(.25, 2, 4), tolerance = 1e-5)
  close(estimated@eta, matrix(.02, 2, 4), tolerance = 1e-5)
  stopifnot(identical(estimated@type_col_map, h@type_col_map),
            identical(metadata$parameter_map$start,
              emhawkes:::build_parameter_map(h, reduced)$start),
            identical(metadata$reduced, reduced), isTRUE(metadata$custom_loglik),
            identical(original, serialize(h, NULL)))
  close(objective(estimated), as.numeric(logLik(fit)))
  close(fit_initial_state(fit)$lambda_component0, matrix(.1, 2, 4))
  # A new order of named coefficients still uses the saved map.
  reordered <- fit
  reordered$estimate <- rev(reordered$estimate)
  close(fitted_hspec(reordered)@alpha, estimated@alpha)
}

# Function environments, shared names, fixed zeros, and generator defaults.
make_model <- function() {
  multiplier <- 2
  new("hspec", mu = .3,
    alpha = function(param = c(shared = .1)) matrix(multiplier * param["shared"]),
    beta = 2,
    impact = function(param = c(shared = .1), n) matrix(param["shared"] * (n > 1)),
    dmark = function(param = c(p = .4), n, mark) param["p"],
    rmark = function(param = c(p = .4, separate = 9), n) param["p"] + param["separate"])
}
h <- make_model()
objective <- function(object) {
  -(object@mu - .5)^2 - (object@alpha() - .6)^2 -
    (object@beta - 2)^2 - (object@impact(n = 2) - .3)^2 -
    (object@dmark(n = 2, mark = c(0, 1)) - .7)^2
}
original <- serialize(h, NULL)
fit <- fit_quiet(h, mylogLik = objective)
estimated <- fitted_hspec(fit)
close(estimated@alpha(), matrix(.6), tolerance = 1e-5)
close(estimated@impact(n = 2), matrix(.3), tolerance = 1e-5)
close(estimated@rmark(n = 2), c(p = 9.7), tolerance = 1e-5)
stopifnot(identical(environment(estimated@alpha), environment(h@alpha)),
          identical(environment(estimated@rmark), environment(h@rmark)),
          identical(original, serialize(h, NULL)),
          sum(names(coef(fit)) == "shared") == 1L)
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
saved <- readRDS(path)
unlink(path)
close(fitted_hspec(saved)@alpha(), estimated@alpha())
close(fitted_hspec(saved)@rmark(n = 2), estimated@rmark(n = 2))
h@alpha <- function(param = numeric(0)) matrix(0)
fit <- fit_quiet(h, mylogLik = function(object) -(object@mu - .5)^2 -
  (object@beta - 2)^2 - (object@impact(n = 2) - .3)^2 -
  (object@dmark(n = 2, mark = c(0, 1)) - .7)^2)
close(fitted_hspec(fit)@alpha(), matrix(0))

# Function-defined dimensions are recomputed as in the optimizer's constructor.
mu_shape <- function(param = c(level = .2)) {
  matrix(param[1], if (param[1] < .5) 1L else 2L)
}
alpha_shape <- function(param = c(level = .2)) {
  n <- if (param[1] < .5) 1L else 2L
  matrix(0, n, n)
}
beta_shape <- function(param = c(level = .2)) {
  n <- if (param[1] < .5) 1L else 2L
  matrix(1, n, n)
}
shape <- new("hspec", mu = mu_shape, alpha = alpha_shape, beta = beta_shape)
objective_shape <- function(object) -sum((object@mu() - 1)^2)
shape_fit <- fit_quiet(shape, mylogLik = objective_shape)
restored_shape <- fitted_hspec(shape_fit)
stopifnot(shape@dimens == 1, restored_shape@dimens == 2)
close(objective_shape(restored_shape), as.numeric(logLik(shape_fit)))

# Residual generators share named distribution estimates without losing closures.
h <- new("hspec", mu = .3, alpha = 0, beta = 1)
h <- set_residual(h, param = c(rate = 1), dresidual = dexp,
                  presidual = pexp, qresidual = qexp, rresidual = rexp)
objective <- function(object) -(object@mu - .3)^2 - object@alpha^2 -
  (object@beta - 1)^2 - (object@qresidual(.5) - qexp(.5, rate = 2))^2
fit <- fit_quiet(h, mylogLik = objective)
estimated <- fitted_hspec(fit)
for (name in c("dresidual", "presidual", "qresidual", "rresidual")) {
  close(eval(formals(methods::slot(estimated, name))$param), c(rate = 2), tolerance = 1e-4)
}
set.seed(29)
sample <- estimated@rresidual(5)
set.seed(29)
close(sample, rexp(5, rate = coef(fit)["rate"]))

# Initial counts inferred from supplied paths survive a state-dependent fit.
h <- new("hspec", mu = function(param = c(base = .2), n, Nc) {
  matrix(param[1] + .1 * Nc[max(1, n - 1), 1])
}, alpha = 0, beta = 1)
marks <- c(0, 2, 1, 3)
dt <- c(0, .2, .3, .4)
fit <- fit_quiet(h, dt, mark = marks, N = matrix(c(4, 5, 6, 7)),
  Nc = matrix(c(5, 7, 8, 11)), lambda_component0 = 0, fixed = c(FALSE, TRUE, TRUE))
initial <- fit_initial_state(fit)
close(initial$N0, 4)
close(initial$Nc0, 5)
close(as.numeric(logLik(fit)), as.numeric(do.call(logLik,
  c(list(object = fitted_hspec(fit), inter_arrival = dt, mark = marks), initial))))

# Metadata survives serialization; accessors do not mutate saved values.
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
saved <- readRDS(path)
unlink(path)
close(fit_initial_state(saved), initial)
close(as.numeric(logLik(saved)), as.numeric(do.call(logLik,
  c(list(object = fitted_hspec(saved), inter_arrival = dt, mark = marks), initial))))
initial$N0 <- 99
close(fit_initial_state(fit)$N0, 4)
changed <- fitted_hspec(fit)
changed@beta[] <- 99
close(fitted_hspec(fit)@beta, matrix(1))

# Finite non-converged fits can be inspected; invalid or older fits fail clearly.
early <- fit_quiet(h, dt, mark = marks, lambda_component0 = 0, iterlim = 0)
stopifnot(methods::is(fitted_hspec(early), "hspec"))
bad <- fit
attr(bad, "emhawkes") <- NULL
fails(fitted_hspec(bad), "Fit metadata is unavailable")
fails(fit_initial_state(bad), "Fit metadata is unavailable")
bad <- fit
bad$estimate[1] <- NA_real_
fails(fitted_hspec(bad), "not finite")
bad <- fit
names(bad$estimate)[1] <- "unknown"
fails(fitted_hspec(bad), "saved parameter map")

}
