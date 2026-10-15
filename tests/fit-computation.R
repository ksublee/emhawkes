library(emhawkes)
previous <- options(emhawkes.engine = "auto")
h <- new("hspec", mu = 1, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .5)
fixed <- c(FALSE, TRUE, TRUE)
fit <- function(object = h, engine = "auto", ...) {
  old <- options(emhawkes.engine = engine)
  on.exit(options(old))
  suppressWarnings(suppressMessages(hfit(object, dt, lambda_component0 = 0,
    fixed = fixed, ...)))
}
computation <- function(x) fit_diagnostics(x)$computation
base <- fit()
info <- computation(base)
stopifnot(info$recorded, info$engine$requested == "auto", identical(info$engine$used, "cpp"),
  is.null(info$engine$fallback_reason), info$gradient$used == "numeric",
  info$hessian$requested == "maxLik", info$hessian$used == "maxLik_numeric")
for (engine in c("r", "cpp")) {
  current <- fit(engine = engine)
  stopifnot(computation(current)$engine$requested == engine,
    identical(computation(current)$engine$used, engine),
    is.null(computation(current)$engine$fallback_reason))
  for (field in c("estimate", "maximum", "gradient", "hessian", "code", "iterations"))
    stopifnot(identical(current[[field]], base[[field]]))
}
analytic <- fit(gradient = "analytic")
stopifnot(computation(analytic)$gradient$used == "analytic",
  identical(computation(analytic)$engine$used, "cpp"))
grad <- function(param) c(3 / param[1] - 1, 0, 0)
hess <- function(param) diag(c(-3 / param[1]^2, 0, 0))
user <- fit(gradient = "analytic", grad = grad, hess = hess)
stopifnot(computation(user)$gradient$requested == "analytic",
  computation(user)$gradient$used == "user", computation(user)$hessian$used == "user")
# Unrelated Hessian attributes must not make recording fail.
tagged_hess <- function(param) {
  value <- hess(param)
  attr(value, "type") <- function() NULL
  value
}
tagged <- fit(hess = tagged_hess)
stopifnot(computation(tagged)$hessian$used == "user")
stable <- fit(final_hessian = "stable")
stopifnot(computation(stable)$hessian$requested == "stable",
  computation(stable)$hessian$used == "scaled_central_Richardson",
  computation(stable)$gradient$used == "numeric")
for (flag in list(FALSE, 0)) {
  disabled <- fit(finalHessian = flag, final_hessian = "stable")
  stopifnot(computation(disabled)$hessian$used == "disabled",
    computation(disabled)$hessian$reason == "finalHessian_false",
    identical(computation(disabled)$hessian$finalHessian, flag))
}
marked <- h
marked@dmark <- function(...) 1
r <- fit(marked)
stopifnot(identical(computation(r)$engine$used, "r"),
  grepl("mark density", computation(r)$engine$fallback_reason))
latent <- h
latent@mu <- function(param = c(rate = 1), lambda_component_n, ...) matrix(param)
r <- fit(latent)
stopifnot(identical(computation(r)$engine$used, "r"),
  grepl("latent", computation(r)$engine$fallback_reason))
classed <- h
classed@mu <- function(param = c(rate = 1)) structure(matrix(param), class = c("custom", "matrix", "array"))
r <- fit(classed)
stopifnot(identical(computation(r)$engine$used, "r"),
  grepl("plain numeric", computation(r)$engine$fallback_reason))
implicit <- h
implicit@impact <- function(n, ...) matrix(0 * sum(..6))
r <- fit(implicit)
stopifnot(identical(computation(r)$engine$used, "r"),
  grepl("latent intensity", computation(r)$engine$fallback_reason))
# A partial native pass is recorded even if its callback output needs R.
sequential <- h
sequential@impact <- function(lambda_component_n, ...) {
  structure(matrix(0), class = c("custom", "matrix", "array"))
}
mixed <- fit(sequential)
stopifnot(identical(computation(mixed)$engine$used, c("cpp", "r")),
  grepl("plain numeric", computation(mixed)$engine$fallback_reason))
# Native likelihoods can also succeed before a later proposal triggers fallback.
late <- h
late@impact <- function(n, mu, ...) {
  value <- matrix(0)
  if (mu[1] > 2) class(value) <- c("custom", "matrix", "array")
  value
}
mixed <- fit(late)
stopifnot(identical(computation(mixed)$engine$used, c("cpp", "r")),
  grepl("plain numeric", computation(mixed)$engine$fallback_reason),
  isTRUE(all.equal(coef(mixed), coef(base), tolerance = 1e-6)))
residual <- h
residual@dresidual <- function(x) dexp(x)
residual@presidual <- function(q) pexp(q)
residual@qresidual <- function(p) qexp(p)
r <- fit(residual, final_hessian = "stable")
stopifnot(identical(computation(r)$engine$used, "r"),
  grepl("residual", computation(r)$engine$fallback_reason),
  computation(r)$hessian$used == "maxLik_numeric",
  computation(r)$hessian$reason == "residual_likelihood")
# Custom likelihood attributes are recorded without evaluating them extra times.
custom <- function(object, param) {
  value <- sum(dnorm(c(.1, .3, .7), param[1], log = TRUE))
  attr(value, "gradient") <- c(1.1 - 3 * param[1], 0, 0)
  attr(value, "hessian") <- diag(c(-3, 0, 0))
  value
}
x <- fit(mylogLik = custom, engine = "not_an_engine", final_hessian = "stable")
stopifnot(is.null(computation(x)$engine$requested), identical(computation(x)$engine$used, "custom"),
  computation(x)$gradient$used == "attribute", computation(x)$hessian$used == "attribute",
  computation(x)$hessian$reason == "custom_likelihood")
x <- fit(mylogLik = custom, grad = grad, hess = hess)
stopifnot(is.na(computation(x)$gradient$used), is.na(computation(x)$hessian$used),
  grepl("precedence", computation(x)$gradient$reason))
by_event <- function(object, param) dnorm(c(.1, .3, .7), param[1], log = TRUE)
x <- fit(mylogLik = by_event, method = "BHHH")
stopifnot(computation(x)$hessian$used == "BHHH")
# Saved metadata is independent of current options and later objective calls.
base_record <- computation(base)
options(emhawkes.engine = "r")
invisible(maxLik::objectiveFn(base)(coef(base)))
stopifnot(identical(base_record, computation(base)))
before <- computation(mixed)
options(emhawkes.engine = "r")
invisible(maxLik::objectiveFn(mixed)(coef(mixed)))
stopifnot(identical(before, computation(mixed)))
path <- tempfile(fileext = ".rds")
saveRDS(mixed, path)
stopifnot(identical(computation(readRDS(path)), before))
unlink(path)
old <- base
attr(old, "emhawkes")$diagnostics$computation <- NULL
bytes <- serialize(old, NULL)
unknown <- computation(old)
stopifnot(!unknown$recorded, is.na(unknown$engine$used),
  unknown$gradient$used == "numeric", is.na(unknown$hessian$used),
  identical(bytes, serialize(old, NULL)))
options(previous)
