library(emhawkes)

prepare <- getFromNamespace("prepare_component_call", "emhawkes")
plan_for <- getFromNamespace("callback_plan", "emhawkes")
legacy <- getFromNamespace("call_component", "emhawkes")
close <- function(actual, expected, tolerance = 1e-12) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
}
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
context <- list(n = 2L, mark = c(0, 2), type = c(0, 1),
                inter_arrival = c(0, .2), N = matrix(c(3, 4)),
                Nc = matrix(c(5, 7)), lambda = matrix(c(.8, .9)),
                lambda_component = matrix(c(.3, .4)),
                lambda_component_n = matrix(.4), mu = matrix(.5),
                alpha = matrix(.1), beta = matrix(1))

# Argument matching, optional defaults and missing() match the previous dispatcher.
callbacks <- list(
  function() 17,
  function(beta, n, mark) list(beta = beta, n = n, mark = mark),
  function(n = 99, optional = 7) {
    list(n = n, optional = optional,
         missing_n = missing(n), missing_optional = missing(optional))
  },
  function(n, unused = stop("do not force optional defaults")) n,
  function(...) list(...),
  function(n, ..., optional = 7) {
    list(n = n, first = ..1, second = ..2, dots = list(...),
         optional = optional, missing_optional = missing(optional))
  },
  function(beta, ..., n) list(beta = beta, n = n, dots = list(...))
)
for (callback in callbacks) {
  plan <- plan_for(callback, names(context), "test")
  prepared <- prepare(plan)
  actual <- do.call(prepared, c(list(FUN = callback), context))
  expected <- legacy(callback, plan, context)
  stopifnot(identical(actual, expected))
}
callback <- function(n, offset = 4) n + offset
prepared <- prepare(plan_for(callback, names(context), "test"))
replacement <- callback
formals(replacement)$offset <- 8
stopifnot(do.call(prepared, c(list(FUN = callback), context)) == 6,
          do.call(prepared, c(list(FUN = replacement), context)) == 10)
callback <- function(n) stop("callback failure", call. = FALSE)
plan <- plan_for(callback, names(context), "test")
fails(do.call(prepare(plan), c(list(FUN = callback), context)), "callback failure")
fails(legacy(callback, plan, context), "callback failure")

# Callbacks may retain promises in closures without evaluating arguments themselves.
# Later loop iterations must not change the values visible to those closures.
snapshots <- function(prepared, dots = FALSE) {
  saved <- list()
  callback <- if (dots) function(...) {
    saved[[length(saved) + 1L]] <<- function() list(...)
    0
  } else function(n, lambda_component_n) {
    saved[[length(saved) + 1L]] <<- function() list(n, lambda_component_n)
    0
  }
  plan <- plan_for(callback, names(context), "test")
  invoke <- prepare(plan)
  current <- context
  for (n in 2:4) {
    lambda_component_n <- matrix(n / 10)
    if (prepared) {
      invoke(callback, n = n, mark = context$mark, type = context$type,
             inter_arrival = context$inter_arrival, N = context$N, Nc = context$Nc,
             lambda = context$lambda, lambda_component = context$lambda_component,
             lambda_component_n = lambda_component_n, mu = context$mu,
             alpha = context$alpha, beta = context$beta)
    } else {
      current$n <- n
      current$lambda_component_n <- lambda_component_n
      legacy(callback, plan, current)
    }
  }
  n <- 99
  lambda_component_n[] <- 99
  lapply(saved, function(f) f())
}
for (dots in c(FALSE, TRUE)) {
  stopifnot(identical(snapshots(TRUE, dots), snapshots(FALSE, dots)))
}

# Independent recursion checks callback order and all changing event context.
# The example has unequal rates, overlapping mapped columns and a tied event time.
dt <- c(0, .2, 0, .4)
types <- c(0, 1, 2, 1)
marks <- c(0, 2, 3, 1)
counts <- cbind(c(2, 3, 3, 4), c(4, 4, 5, 5))
marked_counts <- cbind(c(5, 7, 7, 8), c(6, 6, 9, 9))
alpha <- matrix(c(.06, .04, .03, .02, .05, .01), 2, 3)
beta <- matrix(c(.8, 1.1, 1.7, 2.2, 3.1, 4), 2, 3)
eta <- matrix(c(.01, .02, -.005, .01, .02, -.005), 2, 3)
initial <- matrix(c(.07, .08, .09, .1, .11, .12), 2, 3)
mapping <- list(c(3, 1), c(2, 3))
base_at <- function(n) matrix(c(.3, .4) + .01 * counts[max(1, n - 1), ] +
                               .002 * marked_counts[max(1, n - 1), ])
impact_at <- function(n) matrix(.005 * (counts[n, ] + marks[n]), 2, 3)
density_at <- function(n) exp(-.03 * marks[n])

expected_trace <- list()
record_expected <- function(label, n, total, components, current) {
  expected_trace[[length(expected_trace) + 1L]] <<-
    list(label = label, n = n, lambda = total,
         lambda_component = components, lambda_component_n = current)
}
record_expected("mu", 1, NULL, NULL, NULL)
expected_lambda <- matrix(0, length(dt), 2)
expected_components <- matrix(0, length(dt), 6)
expected_lambda[1, ] <- as.vector(base_at(1)) + rowSums(initial)
expected_components[1, ] <- as.vector(t(initial))
record_expected("mu", 2, expected_lambda, expected_components, initial)
expected_right <- expected_lambda
expected_right_components <- expected_components
right <- initial
expected_loglik <- 0
for (n in 2:length(dt)) {
  current <- right * exp(-beta * dt[n])
  baseline <- base_at(n)
  integrated <- right * (-expm1(-beta * dt[n])) / beta
  expected_lambda[n, ] <- as.vector(baseline) + rowSums(current)
  expected_components[n, ] <- as.vector(t(current))
  expected_loglik <- expected_loglik + log(expected_lambda[n, types[n]]) -
    sum(integrated) - sum(baseline) * dt[n] + log(density_at(n))
  record_expected("dmark", n, expected_lambda, expected_components, current)
  record_expected("impact", n, expected_lambda, expected_components, current)
  selected <- mapping[[types[n]]]
  right <- current
  right[, selected] <- right[, selected] + alpha[, selected] +
    eta[, selected] * (marks[n] - 1) + impact_at(n)[, selected]
  record_expected("mu", n + 1, expected_lambda, expected_components, current)
  expected_right[n, ] <- as.vector(base_at(n + 1)) + rowSums(right)
  expected_right_components[n, ] <- as.vector(t(right))
}

observed_trace <- list()
record_actual <- function(label, n, N, Nc, lambda, lambda_component, lambda_component_n) {
  close(unname(N), counts)
  close(unname(Nc), marked_counts)
  observed_trace[[length(observed_trace) + 1L]] <<-
    list(label = label, n = n, lambda = if (is.null(lambda)) NULL else unname(lambda),
         lambda_component = if (is.null(lambda_component)) NULL else unname(lambda_component),
         lambda_component_n = lambda_component_n)
}
mu <- function(n, N, Nc, lambda, lambda_component, lambda_component_n) {
  record_actual("mu", n, N, Nc, lambda, lambda_component, lambda_component_n)
  base_at(n)
}
impact <- function(n, N, Nc, lambda, lambda_component, lambda_component_n,
                   mark, type, inter_arrival, alpha, beta) {
  close(mark, marks)
  close(type, types)
  close(inter_arrival, dt)
  record_actual("impact", n, N, Nc, lambda, lambda_component, lambda_component_n)
  impact_at(n)
}
dmark <- function(n, N, Nc, lambda, lambda_component, lambda_component_n) {
  record_actual("dmark", n, N, Nc, lambda, lambda_component, lambda_component_n)
  density_at(n)
}
h <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = eta,
         impact = impact, dmark = dmark, type_col_map = mapping)
previous <- options(emhawkes.engine = "r")
for (infer in c(FALSE, TRUE)) {
  observed_trace <- list()
  actual <- logLik(h, dt, type = types, mark = marks, N0 = counts[1, ],
                   Nc0 = marked_counts[1, ], lambda_component0 = initial, infer = infer)
  close(observed_trace, expected_trace)
  close(as.numeric(if (infer) actual$loglikelihood else actual), expected_loglik)
  if (infer) {
    close(unname(actual$lambda), expected_lambda)
    close(unname(actual$lambda_component), expected_components)
    close(unname(actual$rambda), expected_right)
    close(unname(actual$rambda_component), expected_right_components)
  }
}
options(previous)


# Cached numeric kernels span multiple blocks without changing mapped updates.
# With 64 components, 1,029 events cross a 65,536-entry block boundary.
evaluate <- function(engine, object, intervals, types = NULL, initial = 0, infer = TRUE) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  logLik(object, intervals, type = types, lambda_component0 = initial, infer = infer)
}
kernel_names <- list(c("target1", "target2"), paste0("kernel", 1:32))
block_beta <- matrix(seq(.8, 3, length.out = 64), 2, 32, dimnames = kernel_names)
block_h <- new("hspec", mu = c(.3, .4),
               alpha = matrix(seq(.001, .003, length.out = 64), 2, 32),
               beta = block_beta, type_col_map = list(seq(1, 31, 2), seq(2, 32, 2)))
block_dt <- c(0, rep(c(.01, .04, 0, .07), length.out = 1029))
block_type <- c(0, rep(1:2, length.out = 1029))
block_reference <- evaluate("cpp", block_h, block_dt, block_type, initial = .02)
close(evaluate("r", block_h, block_dt, block_type, initial = .02), block_reference)
block_h@impact <- function(n, beta, lambda_component_n, ...) {
  stopifnot(identical(dimnames(beta), kernel_names),
            identical(dimnames(lambda_component_n), kernel_names))
  matrix(0, 2, 32)
}
close(evaluate("r", block_h, block_dt, block_type, initial = .02), block_reference)
empty_callback_h <- block_h
empty_callback_h@impact <- function() matrix(0, 2, 32)
close(evaluate("r", empty_callback_h, block_dt, block_type, initial = .02),
      block_reference)

# An underflowed beta * dt still has integral dt, including in cached blocks.
tiny_h <- new("hspec", mu = 0, alpha = 0, beta = 1e-300)
tiny_dt <- c(0, rep(c(0, 1e-100), 32))
tiny_expected <- 64 * log(1e100) - 32
for (infer in c(FALSE, TRUE)) {
  actual <- evaluate("r", tiny_h, tiny_dt, initial = 1e100, infer = infer)
  close(as.numeric(if (infer) actual$loglikelihood else actual), tiny_expected)
}

# Zero components contribute zero even at the largest supported interval scales.
# Their weights may round to infinity near the edge of the floating-point range.
zero_h <- new("hspec", mu = 1e-308, alpha = 0, beta = 5e-324)
large_dt <- c(0, rep(.Machine$double.xmax, 32))
zero_expected <- 32 * (log(1e-308) - 1e-308 * .Machine$double.xmax)
for (infer in c(FALSE, TRUE)) {
  actual <- evaluate("r", zero_h, large_dt, infer = infer)
  close(as.numeric(if (infer) actual$loglikelihood else actual), zero_expected)
}

# Classed numeric inputs retain their R dispatch; they cannot be flattened into
# the plain numeric cache. Both the beta matrix and scalar interval are exercised.
math_calls <- character()
Math.dispatch_kernel <- function(x, ...) {
  math_calls <<- c(math_calls, .Generic)
  NextMethod()
}
Ops.dispatch_kernel <- function(e1, e2) NextMethod()
class_h <- block_h
class_h@beta <- function(param = numeric(0)) {
  structure(block_beta, class = c("dispatch_kernel", "matrix", "array"))
}
short_dt <- block_dt[1:40]
short_type <- block_type[1:40]
short_reference <- evaluate("r", block_h, short_dt, short_type, initial = .02)
close(evaluate("r", class_h, short_dt, short_type, initial = .02), short_reference)
stopifnot("expm1" %in% math_calls)

interval_kernel_calls <- 0L
Ops.dispatch_interval <- function(e1, e2) {
  if (.Generic == "*" && !missing(e2) &&
      identical(dim(e1), c(2L, 32L)) && length(e2) == 1L) {
    interval_kernel_calls <<- interval_kernel_calls + 1L
  }
  NextMethod()
}
`[.dispatch_interval` <- function(x, i, ...) {
  structure(NextMethod("["), class = "dispatch_interval")
}
classed_dt <- structure(short_dt, class = "dispatch_interval")
close(evaluate("r", block_h, classed_dt, short_type, initial = .02), short_reference)
stopifnot(interval_kernel_calls >= length(short_dt) - 1L)


# Invalid proposals retain their special values without log-domain warnings.
without_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  stopifnot(length(warnings) == 0L)
  value
}
for (rate in c(0, -1, NA_real_, NaN)) {
  special_h <- new("hspec", mu = function(n) matrix(if (n == 1) 1 else rate),
                   alpha = 0, beta = 1)
  actual <- without_warnings(evaluate("r", special_h, c(0, .25), infer = FALSE))
  expected <- suppressWarnings(log(rate) - rate * .25)
  stopifnot(identical(as.numeric(actual), expected))
}

# Custom Math methods can warn even for positive intensities. They still run,
# and their log warnings are suppressed just as in the previous likelihood.
custom_log_calls <- 0L
Math.dispatch_rate <- function(x, ...) {
  if (.Generic == "log") {
    custom_log_calls <<- custom_log_calls + 1L
    warning("custom log warning")
    return(unclass(NextMethod()))
  }
  NextMethod()
}
Ops.dispatch_rate <- function(e1, e2) NextMethod()
`[.dispatch_rate` <- function(x, i, ...) {
  structure(NextMethod("["), class = "dispatch_rate")
}
custom_h <- new("hspec", mu = function(n) {
  structure(matrix(.4), class = c("dispatch_rate", "matrix", "array"))
}, alpha = 0, beta = 1)
actual <- without_warnings(evaluate("r", custom_h, c(0, .25), infer = FALSE))
close(as.numeric(actual), log(.4) - .4 * .25)
stopifnot(custom_log_calls == 1L)
