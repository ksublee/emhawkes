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
