# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

fit_engine <- function(engine, ...) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(hfit(...))
}
set.seed(105)
truth <- new("hspec", mu = c(.3, .25), alpha = matrix(c(.15, .05, .04, .1), 2),
              beta = matrix(1.5, 2, 2))
path <- hsim(truth, size = 701, lambda_component0 = .1)
start <- truth
start@mu[] <- .2
map <- getFromNamespace("build_parameter_map", "emhawkes")
for (reduced in c(TRUE, FALSE)) {
  parameters <- map(start, reduced)$start
  fixed <- startsWith(names(parameters), "beta")
  fit <- function(engine) fit_engine(engine, start, path$inter_arrival, path$type,
    lambda_component0 = .1, reduced = reduced, fixed = fixed,
    constraints = list(ineqA = diag(length(parameters)), ineqB = rep(-1e-6, length(parameters))),
    control = list(iterlim = 200, reltol = 1e-9))
  reference <- fit("r")
  native <- fit("cpp")
  stopifnot(inherits(native, "maxLik"), reference$code == native$code,
            identical(names(coef(reference)), names(coef(native))),
            abs(as.numeric(logLik(reference)) - as.numeric(logLik(native))) < 1e-7,
            max(abs(coef(reference) - coef(native))) < 1e-4)
}

# Custom derivatives still pass through to maxLik. The Poisson special case
# has a closed-form optimum, independent of the Hawkes recurrence.
calls <- c(gradient = 0L, hessian = 0L)
gradient <- function(param) {
  calls["gradient"] <<- calls["gradient"] + 1L
  c(5 / param[1] - 3, 0, 0)
}
hessian <- function(param) {
  calls["hessian"] <<- calls["hessian"] + 1L
  diag(c(-5 / param[1]^2, 0, 0))
}
poisson <- new("hspec", mu = 1, alpha = 0, beta = 1)
fit <- fit_engine("cpp", poisson, c(0, .2, .3, .5, 1, 1), lambda_component0 = 0,
                   fixed = c(FALSE, TRUE, TRUE), grad = gradient, hess = hessian, method = "NR")
stopifnot(all(calls > 0), abs(coef(fit)[1] - 5 / 3) < 1e-7)

}
