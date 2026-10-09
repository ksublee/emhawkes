library(emhawkes)

close <- function(x, y, tol = 1e-9) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
warns <- function(expr, pattern) {
  seen <- FALSE
  result <- withCallingHandlers(force(expr), warning = function(w) {
    stopifnot(grepl(pattern, conditionMessage(w), fixed = TRUE))
    seen <<- TRUE
    invokeRestart("muffleWarning")
  })
  stopifnot(seen)
  result
}
build_map <- getFromNamespace("build_parameter_map", "emhawkes")
restore <- getFromNamespace("restore_parameters", "emhawkes")

# Known coefficient order and locations, including non-square kernel matrices.
alpha <- matrix(c(.1, .2, .1, .3, .2, .4, .3, .4), 2, byrow = TRUE)
h <- new("hspec", mu = c(.2, .3), alpha = alpha, beta = matrix(2, 2, 4),
         type_col_map = list(c(1, 3), c(2, 4)))
reduced <- build_map(h)
stopifnot(identical(names(reduced$start),
  c("mu1", "mu2", "alpha1.1", "alpha1.2", "alpha1.4", "alpha2.2", "beta1.1")))
pr <- reduced$start
pr[3:6] <- c(10, 20, 30, 40)
decoded <- restore(reduced, pr)
close(decoded$alpha, matrix(c(10, 20, 10, 30, 20, 40, 30, 40), 2, byrow = TRUE))
close(restore(reduced, rev(pr))$alpha, decoded$alpha)
full <- build_map(h, reduced = FALSE)
stopifnot(identical(names(full$start)[3:10],
  c("alpha1.1", "alpha2.1", "alpha1.2", "alpha2.2", "alpha1.3", "alpha2.3", "alpha1.4", "alpha2.4")))
pr <- full$start
pr[3:10] <- 11:18
close(restore(full, pr)$alpha, matrix(11:18, 2, 4))
close(restore(full, full$start)$alpha, alpha)

# Function defaults keep their environments, and shared names remain one parameter.
make_alpha <- function() {
  multiplier <- 2
  function(param = c(shared = .1)) matrix(multiplier * param["shared"])
}
alpha_fun <- make_alpha()
impact <- function(param = c(shared = .1), n) matrix(param["shared"] * (n > 1))
fun_h <- new("hspec", mu = .2, alpha = alpha_fun, beta = 2, impact = impact)
map <- build_map(fun_h)
stopifnot(sum(names(map$start) == "shared") == 1)
pr <- map$start
pr["shared"] <- .3
decoded <- restore(map, pr)
close(decoded$alpha(), .6)
close(decoded$impact(n = 2), .3)
close(alpha_fun(), .2)
stopifnot(identical(environment(decoded$alpha), environment(alpha_fun)))
conflict <- fun_h
conflict@impact <- function(param = c(shared = .9), n) matrix(param["shared"])
conflict_map <- warns(build_map(conflict), "different starting values")
close(conflict_map$start["shared"], .1)
bad <- fun_h
for (param in list(c(.1, .2), c(a = .1, a = .2), c(a = NA_real_), c(a = Inf))) {
  bad@impact <- impact
  formals(bad@impact)$param <- param
  fails(build_map(bad), "unique, nonempty names")
}

# A parameterless component is fixed; an explicit empty param is also valid.
fixed <- fun_h
fixed@impact <- function(n) matrix(0)
stopifnot(!("shared" %in% names(build_map(fixed)$components$impact$param)))
fixed@impact <- function(param = numeric(0), n) matrix(0)
close(restore(build_map(fixed), build_map(fixed)$start)$impact(n = 2), 0)
index_only <- new("hspec", mu = function(n) matrix(.2 + .1 * (n > 2)), alpha = 0, beta = 1)
index_out <- logLik(index_only, c(0, .2, .3), infer = TRUE, showWarning = FALSE)
close(index_out$lambda[-1, ], c(.2, .3))
set.seed(16)
index_path <- hsim(index_only, size = 5)
index_restored <- logLik(index_only, index_path$inter_arrival, infer = TRUE, showWarning = FALSE)
close(index_path$lambda, index_restored$lambda)

# Custom likelihoods may omit observations, use defaults, and accept dots.
objective <- function(object, offset = 7, ...) {
  offset - sum((object@mu - c(.4, .6))^2) - sum((object@alpha - .2)^2) - sum((object@beta - 1.5)^2)
}
fit <- hfit(h, mylogLik = objective, iterlim = 100)
stopifnot(inherits(fit, "maxLik"), identical(names(coef(fit)), names(reduced$start)))
close(as.numeric(logLik(fit)), 7, tol = 1e-7)
close(coef(fit)[1:2], c(.4, .6), tol = 1e-5)
fails(hfit(h, mylogLik = function(object, unavailable) 0), "unavailable")

# State-dependent Poisson model: initial marked counts must reach hfit's likelihood.
mu <- function(param = c(base = .2), n, Nc) matrix(param[1] + .1 * Nc[max(1, n - 1), 1])
state <- new("hspec", mu = mu, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .4)
types <- c(0, 1, 1, 1)
marks <- c(0, 2, 1, 3)
N0 <- 4
Nc0 <- 5
N <- matrix(c(4, 5, 6, 7), ncol = 1)
Nc <- matrix(c(5, 7, 8, 11), ncol = 1)
baseline <- .2 + .1 * c(5, 7, 8)
expected <- sum(log(baseline) - baseline * dt[-1])
auto <- logLik(state, dt, type = types, mark = marks, N0 = N0, Nc0 = Nc0,
               lambda_component0 = 0, infer = TRUE)
explicit <- logLik(state, dt, type = types, mark = marks, N = N, Nc = Nc,
                   lambda_component0 = matrix(0), infer = TRUE)
close(auto$loglikelihood, expected)
close(explicit$loglikelihood, expected)
close(auto$lambda[-1, ], baseline)
close(auto$rambda, .2 + .1 * Nc)
fit <- hfit(state, dt, type = types, mark = marks, N0 = N0, Nc0 = Nc0,
            lambda_component0 = 0, fixed = c(FALSE, TRUE, TRUE), iterlim = 0)
close(as.numeric(logLik(fit)), expected)
stopifnot(identical(names(coef(fit)), c("base", "alpha1", "beta1")))
estimated <- hfit(state, dt, type = types, mark = marks, N0 = N0, Nc0 = Nc0,
                  lambda_component0 = 0, fixed = c(FALSE, TRUE, TRUE))
optimal_base <- uniroot(function(base) sum(1 / (base + .1 * c(5, 7, 8))) - sum(dt),
                        c(.01, 10), tol = 1e-10)$root
close(coef(estimated)["base"], optimal_base, tol = 1e-5)
optimal_rates <- optimal_base + .1 * c(5, 7, 8)
close(as.numeric(logLik(estimated)), sum(log(optimal_rates) - optimal_rates * dt[-1]), tol = 1e-7)
custom <- function(object, inter_arrival, type, mark, Nc, Nc0, lambda_component0, multiplier = 1, ...) {
  stopifnot(Nc0 == 5, Nc[1, 1] == 5)
  multiplier * logLik(object, inter_arrival, type, mark, Nc = Nc, Nc0 = Nc0,
                      lambda_component0 = lambda_component0)
}
custom_fit <- hfit(state, dt, type = types, mark = marks, Nc0 = Nc0,
                   lambda_component0 = 0, mylogLik = custom,
                   fixed = c(FALSE, TRUE, TRUE), iterlim = 0)
close(as.numeric(logLik(custom_fit)), expected)

# Component functions can declare only the arguments they use, without dots.
set.seed(17)
path <- hsim(state, size = 20, N0 = N0, Nc0 = Nc0, lambda_component0 = 0)
inferred <- logLik(state, path$inter_arrival, type = path$type, mark = path$mark,
                   N0 = N0, Nc0 = Nc0, lambda_component0 = 0, infer = TRUE)
close(path$lambda, inferred$lambda)
close(path$rambda, inferred$rambda)
mark_h <- new("hspec", mu = .5, alpha = .1, beta = 1,
  impact = function(param = c(extra = .01), n, mark) matrix(param[1] * mark[n]),
  rmark = function(n) 2,
  dmark = function(param = c(p = .5), n, N, Nc, lambda) {
    stopifnot(N[n, 1] == n - 1, Nc[n, 1] == 2 * (n - 1), is.matrix(lambda))
    param[1]
  })
set.seed(18)
path <- hsim(mark_h, size = 10, lambda_component0 = .1)
stopifnot(all(path$mark[-1] == 2))
out <- logLik(mark_h, path$inter_arrival, type = path$type, mark = path$mark,
              lambda_component0 = .1, infer = TRUE)
close(out$lambda, path$lambda)
close(out$rambda, path$rambda)
plain <- mark_h
plain@dmark <- NULL
close(logLik(mark_h, path$inter_arrival, type = path$type, mark = path$mark,
             lambda_component0 = .1),
      logLik(plain, path$inter_arrival, type = path$type, mark = path$mark,
             lambda_component0 = .1) + 9 * log(.5))

# Default initialization remains dependent on each proposed parameter value.
poisson <- new("hspec", mu = .3, alpha = .2, beta = 1)
from_default <- function(object, inter_arrival, ...) {
  logLik(object, inter_arrival, showWarning = FALSE)
}
fit <- hfit(poisson, dt, fixed = c(FALSE, TRUE, TRUE), iterlim = 0)
close(as.numeric(logLik(fit)), logLik(poisson, dt, showWarning = FALSE))
fit2 <- hfit(poisson, dt, mylogLik = from_default,
             fixed = c(FALSE, TRUE, TRUE), iterlim = 0)
close(as.numeric(logLik(fit2)), as.numeric(logLik(fit)))
state_default <- suppressMessages(logLik(state, dt, mark = marks, Nc0 = 5, infer = TRUE))
stopifnot(all(is.finite(state_default$lambda)))

# Invalid observations and callback contracts fail before optimization calls them.
called <- FALSE
should_not_run <- function(object, ...) { called <<- TRUE; 0 }
fails(hfit(poisson, c(0, -1), mylogLik = should_not_run), "inter_arrival")
stopifnot(!called)
fails(logLik(state, dt, type = c(0, 1, 2, 1)), "type")
fails(logLik(state, dt, mark = c(0, 1)), "mark")
fails(logLik(state, dt, Nc = matrix(0, 4, 2)), "one column per dimension")
fails(logLik(state, dt, Nc = matrix(c(5, 6, 8, 9), 4), mark = marks), "inconsistent")
fails(logLik(state, dt, Nc = Nc, Nc0 = 0, mark = marks), "inconsistent")
fails(logLik(state, dt, N0 = c(1, 2)), "N0")
fails(hsim(state, size = 2, Nc0 = Inf), "Nc0")
unsupported <- state
unsupported@mu <- function(param = c(base = .2), covariate) matrix(param[1])
fails(logLik(unsupported, dt, lambda_component0 = 0), "covariate")
fails(hsim(unsupported, size = 2, lambda_component0 = 0), "covariate")

# Scalar initial excitation and a full matrix describe the same initial state.
set.seed(19)
scalar <- hsim(h, size = 20, lambda_component0 = .1, N0 = 3, Nc0 = 4)
set.seed(19)
matrix_path <- hsim(h, size = 20, lambda_component0 = matrix(.1, 2, 4),
                    N0 = c(3, 3), Nc0 = c(4, 4))
close(as.matrix(scalar), as.matrix(matrix_path))
close(logLik(h, scalar$inter_arrival, type = scalar$type, lambda_component0 = .1),
      logLik(h, scalar$inter_arrival, type = scalar$type, lambda_component0 = matrix(.1, 2, 4)))
legacy <- warns(hsim(h, size = 1, lambda_component0 = c(.1, .2)), "deprecated")
close(legacy$rambda_component[1, ], as.vector(t(matrix(rep(c(.1, .2), 4), 2, 4))))
fails(logLik(h, c(0, .1), type = c(0, 1), lambda_component0 = matrix(.1, 4, 2)), "dimensions")

# Only the initial row is omitted; tied event times retain their input order.
close(logLik(state, 0, mark = 0, Nc0 = 5, lambda_component0 = 0), 0)
tied <- logLik(state, c(0, 0, .2), mark = c(0, 2, 1), Nc0 = 5,
               lambda_component0 = 0, infer = TRUE)
close(tied$lambda[-1, ], c(.7, .9))
close(tied$loglikelihood, log(.7) + log(.9) - .9 * .2)
