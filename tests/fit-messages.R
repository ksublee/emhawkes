library(emhawkes)

capture_messages <- function(expr) {
  messages <- character(0)
  value <- withCallingHandlers(force(expr), message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = messages)
}

h <- new("hspec", mu = c(.2, .2), alpha = matrix(c(.1, .2, .2, .1), 2),
         beta = matrix(1, 2, 2))
evaluations <- 0L
objective <- function(object) {
  evaluations <<- evaluations + 1L
  -sum((object@mu - .2)^2) - sum((object@alpha - h@alpha)^2) - sum((object@beta - 1)^2)
}
out <- capture_messages(hfit(h, mylogLik = objective))
stopifnot(evaluations > 1L, length(out$messages) == 1L,
          grepl("within each matrix", out$messages, fixed = TRUE),
          grepl("mu, alpha, beta", out$messages, fixed = TRUE),
          grepl("reduced = FALSE", out$messages, fixed = TRUE),
          inherits(out$value, "maxLik"), length(coef(out$value)) == 4L)

# Separate matrix estimates do not produce a sharing message.
out <- capture_messages(hfit(h, reduced = FALSE, mylogLik = objective))
stopifnot(length(out$messages) == 0L, length(coef(out$value)) == 10L)
out <- capture_messages(suppressMessages(hfit(h, mylogLik = objective)))
stopifnot(length(out$messages) == 0L, length(coef(out$value)) == 4L)

# No message is needed when all entries within each matrix are distinct.
distinct <- new("hspec", mu = c(.2, .3), alpha = matrix(c(.1, .2, .3, .4), 2),
                beta = matrix(1:4, 2))
out <- capture_messages(hfit(distinct, mylogLik = function(object) {
  -sum((object@mu - distinct@mu)^2) - sum((object@alpha - distinct@alpha)^2) -
    sum((object@beta - distinct@beta)^2)
}))
stopifnot(length(out$messages) == 0L)

# Named function parameters are a different sharing rule; reduced=FALSE does not separate them.
functions <- new("hspec", mu = .2, alpha = .1, beta = 1,
  impact = function(param = c(shared = .1), n) matrix(param[1]),
  dmark = function(param = c(shared = .1), n) param[1])
out <- capture_messages(hfit(functions, mylogLik = function(object) {
  -sum((object@mu - .2)^2) - sum((object@alpha - .1)^2) - sum((object@beta - 1)^2) -
    (object@impact(n = 2) - .1)^2
}))
stopifnot(length(out$messages) == 0L, sum(names(coef(out$value)) == "shared") == 1L)

# Invalid data fail before an estimation-start message is emitted.
out <- capture_messages(tryCatch(hfit(h, c(0, -1), type = c(0, 1)), error = identity))
stopifnot(inherits(out$value, "error"), length(out$messages) == 0L)

# Failed default numerical curvature gives one hint, without changing the fit.
poisson <- new("hspec", mu = 1, alpha = 0, beta = 1)
dt <- c(0, .5, 1.5, .25)
fixed <- c(FALSE, TRUE, TRUE)
valid <- suppressMessages(hfit(poisson, dt, lambda_component0 = 0, fixed = fixed))
failed <- valid
failed$hessian[1, 1] <- -failed$hessian[1, 1]
stopifnot(is.nan(suppressWarnings(miscTools::stdEr(failed)[1])))
finalize <- getFromNamespace("finalize_hawkes_fit", "emhawkes")
check_message <- function(fit = failed, object = poisson, hess = NULL,
                          additional = list(), method = "BFGS") {
  capture_messages(finalize(fit, object, constraints = NULL, hess = hess,
                             additional = additional, method = method,
                             final_hessian = "maxLik"))
}
out <- check_message()
stopifnot(length(out$messages) == 1L,
          grepl('final_hessian = "auto"', out$messages, fixed = TRUE),
          identical(out$value$fit, failed),
          identical(out$value$diagnostics$status, "not_checked"))
out <- capture_messages(suppressMessages(check_message()))
stopifnot(length(out$messages) == 0L)
stopifnot(length(check_message(valid)$messages) == 0L,
          length(check_message(hess = function(param) failed$hessian)$messages) == 0L,
          length(check_message(method = "BHHH")$messages) == 0L,
          length(check_message(additional = list(finalHessian = "bhhh"))$messages) == 0L,
          length(check_message(additional = list(bhhhHessian = TRUE))$messages) == 0L)
disabled <- suppressMessages(hfit(poisson, dt, lambda_component0 = 0,
                                   fixed = fixed, finalHessian = FALSE))
stopifnot(length(check_message(disabled)$messages) == 0L)
residual <- new("hspec", mu = 1, alpha = 0, beta = 1,
  dresidual = function(x) dexp(x), presidual = function(q) pexp(q),
  qresidual = function(p) qexp(p))
stopifnot(length(check_message(object = residual)$messages) == 0L)
# Custom objectives cannot use the Hawkes stabilization and retain their messages.
custom <- capture_messages(hfit(poisson, fixed = fixed,
  mylogLik = function(object) -(object@mu - 1)^2,
  hess = function(param) diag(c(2, 0, 0))))
stopifnot(is.nan(suppressWarnings(miscTools::stdEr(custom$value)[1])),
          length(custom$messages) == 0L)
