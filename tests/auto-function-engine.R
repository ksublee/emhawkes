library(emhawkes)
old_options <- options(emhawkes.engine = "auto")
native <- emhawkes:::hawkes_loglik_cpp
calls <- 0L
assignInNamespace("hawkes_loglik_cpp", function(...) {
  calls <<- calls + 1L
  native(...)
}, "emhawkes")
fit <- function(model, engine = "auto") {
  options(emhawkes.engine = engine)
  suppressWarnings(suppressMessages(hfit(model, c(0, .2, .4, .3, .8),
    lambda_component0 = .1, method = "BFGS", control = list(iterlim = 3))))
}
compare <- function(model, supported = TRUE, probes = 0L) {
  calls <<- 0L
  auto <- fit(model)
  stopifnot(if (supported) calls > 0L else calls == probes)
  reference <- fit(model, "r")
  for (name in c("estimate", "maximum", "gradient", "hessian", "code", "iterations"))
    stopifnot(isTRUE(all.equal(auto[[name]], reference[[name]], tolerance = 1e-7)))
  if (supported) {
    compiled <- fit(model, "cpp")
    for (name in c("estimate", "maximum", "gradient", "hessian", "code", "iterations"))
      stopifnot(identical(auto[[name]], compiled[[name]]))
  }
}
h <- new("hspec", mu = function(param = c(rate = .3)) matrix(param), alpha = .05, beta = 1.2)
compare(h)
h@impact <- function(param = c(jump = .01), n, N, ...) matrix(param * (1 + N[n, 1] %% 2))
compare(h)
h@mu <- function(param = c(rate = .3), n, Nc, ...) matrix(param + .001 * Nc[max(1, n - 1), 1])
h@impact <- function(param = c(jump = .01), lambda_component_n, ...) matrix(param - .02 * lambda_component_n)
compare(h)
h@mu <- function(param = c(rate = .3), lambda_component_n, ...) matrix(param)
compare(h, FALSE)
h@mu <- function(param = c(rate = .3)) matrix(param)
h@dmark <- function(...) 1
compare(h, FALSE)
h@dmark <- NULL
h@mu <- function(param = c(rate = .3)) structure(matrix(param), class = c("custom", "matrix", "array"))
compare(h, FALSE)
h@mu <- function(param = c(rate = .3)) matrix(param)
h@impact <- function(param = c(jump = .01), n, ...) {
  structure(matrix(param), class = c("custom", "matrix", "array"))
}
compare(h, FALSE)
h@impact <- function(param = c(jump = .01), lambda_component_n, ...) {
  structure(matrix(param), class = c("custom", "matrix", "array"))
}
compare(h, FALSE, probes = 1L)
# An unnamed latent argument switches permanently to R.
h@mu <- function(param = c(rate = .3)) matrix(param)
h@impact <- function(param = c(jump = .01), ...) matrix(param + 0 * sum(..7))
unsupported <- emhawkes:::stop_cpp_unsupported
fallbacks <- 0L
assignInNamespace("stop_cpp_unsupported", function(message) {
  fallbacks <<- fallbacks + 1L
  unsupported(message)
}, "emhawkes")
calls <- 0L
auto <- fit(h)
stopifnot(calls == 0L, fallbacks == 1L)
assignInNamespace("stop_cpp_unsupported", unsupported, "emhawkes")
reference <- fit(h, "r")
stopifnot(identical(auto$estimate, reference$estimate), identical(auto$maximum, reference$maximum))
# The same capability fallback covers latent baseline reads through dots.
h@impact <- NULL
h@mu <- function(param = c(rate = .3), n, ...) matrix(param + 0 * sum(..6))
compare(h, FALSE)
h@mu <- function(param = c(rate = .3)) matrix(param)
# User errors propagate rather than causing an R retry.
failures <- 0L
h@impact <- function(param = c(jump = .01), n, ...) {
  failures <<- failures + 1L
  stop("callback failure")
}
error <- tryCatch(fit(h), error = identity)
stopifnot(inherits(error, "error"), grepl("callback failure", conditionMessage(error)), failures == 1L)
assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
options(old_options)
