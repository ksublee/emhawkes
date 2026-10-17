# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

capture_conditions <- function(expr) {
  messages <- warnings <- list()
  value <- withCallingHandlers(tryCatch(force(expr), error = identity),
    message = function(m) { messages[[length(messages) + 1L]] <<- m; invokeRestart("muffleMessage") },
    warning = function(w) { warnings[[length(warnings) + 1L]] <<- w; invokeRestart("muffleWarning") })
  list(value = value, messages = messages, warnings = warnings)
}
text <- function(conditions) vapply(conditions, conditionMessage, character(1))
has <- function(string, ...) all(vapply(list(...), grepl, logical(1), x = string, fixed = TRUE))
old <- options(emhawkes.engine = "auto")

# Mapped numeric models now support mean initialization before RNG draws.
h <- new("hspec", mu = c(.2, .3), alpha = matrix(.01, 2, 4), beta = matrix(1, 2, 4),
  type_col_map = list(c(1, 3), c(2, 4)))
set.seed(508)
rng <- .Random.seed
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  out <- capture_conditions(hsim(h, 1))
  stopifnot(!inherits(out$value, "error"), identical(rng, .Random.seed),
    length(out$messages) == 1L, has(text(out$messages), "auto -> mean"), !length(out$warnings))
}
out <- capture_conditions(logLik(h, c(0, .2), c(0, 1), lambda_component0 = matrix(.1, 4, 2)))
stopifnot(has(conditionMessage(out$value), "2 x 4", "received 4 x 2"))
out <- capture_conditions(new("hspec", mu = c(.2, .3), alpha = matrix(.1, 3, 3), beta = matrix(1, 3, 3)))
stopifnot(has(conditionMessage(out$value), "rows", "mu=2", "alpha=3"))
out <- capture_conditions(new("hspec", mu = .2, alpha = matrix(.1, 1, 2), beta = matrix(1, 1, 3)))
stopifnot(has(conditionMessage(out$value), "columns", "alpha=2", "beta=3"))
evaluations <- 0L
out <- capture_conditions(new("hspec", mu = function(param = c(m = .2)) {
  evaluations <<- evaluations + 1L; param
}, alpha = .1, beta = 1))
stopifnot(has(conditionMessage(out$value), "mu(param)", "matrix"), evaluations == 2L)
out <- capture_conditions(new("hspec", mu = .2, alpha = 1, beta = 1, stability_check = TRUE))
stopifnot(length(out$warnings) == 1L, has(text(out$warnings), "spectral radius = 1", "must be < 1"))
out <- capture_conditions(new("hspec", mu = .2, alpha = .1, beta = 1, eta = .01, stability_check = TRUE))
stopifnot(length(out$warnings) == 1L, has(text(out$warnings), "Stability was not checked"))

# Singular starting parameters select zero and explain the actual state once.
singular <- new("hspec", mu = c(.2, .3), alpha = diag(2), beta = matrix(1, 2, 2))
args <- list(inter_arrival = c(0, .2, .3), type = c(0, 1, 2), infer = TRUE)
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  automatic <- capture_conditions(do.call(logLik, c(list(object = singular), args)))
  explicit <- do.call(logLik, c(list(object = singular, lambda_component0 = 0), args))
  stopifnot(identical(automatic$value, explicit), length(automatic$messages) == 1L,
    !length(automatic$warnings),
    has(text(automatic$messages), "auto -> zero", "lambda_component0", "explicitly is recommended"))
}
functions <- singular
functions@alpha <- function(param = c(a = 1)) diag(param[1], 2)
map <- emhawkes:::build_parameter_map(functions)
observation <- emhawkes:::normalize_observation(functions, args$inter_arrival, args$type)
objective <- emhawkes:::make_hawkes_objective(functions, observation, NULL, map, FALSE,
  derivatives = FALSE)
out <- capture_conditions(objective$loglik(map$start))
stopifnot(!length(out$warnings),
  identical(as.numeric(out$value), as.numeric(explicit$loglikelihood)))

# Verbose custom or invalid fits must not claim automatic Hawkes initialization.
poisson <- new("hspec", mu = 1, alpha = 0, beta = 1)
printed <- capture.output(out <- capture_conditions(hfit(poisson, verbose = TRUE,
  fixed = c(FALSE, TRUE, TRUE), iterlim = 0, mylogLik = function(object) -(object@mu - 1)^2)))
stopifnot(!length(out$messages), inherits(out$value, "maxLik"))
out <- capture_conditions(hfit(poisson, c(0, -1), verbose = TRUE))
stopifnot(inherits(out$value, "error"), !length(out$messages))
printed <- capture.output(out <- capture_conditions(hfit(poisson, c(0, .2, .3), verbose = TRUE,
  fixed = c(FALSE, TRUE, TRUE), iterlim = 0)))
stopifnot(length(out$messages) == 1L, has(text(out$messages), "auto -> mean"))

# Deprecations keep their warning class and point to the replacement API.
out <- capture_conditions(hsim(poisson, 1, lambda0 = 0))
stopifnot(length(out$warnings) == 1L, is.null(conditionCall(out$warnings[[1]])),
  has(text(out$warnings), "deprecated", "lambda0", "lambda_component0"),
  identical(out$value, hsim(poisson, 1, lambda_component0 = 0)))

# Early stopping reports returned rows, including the initial row. The message,
# result and RNG state must match R for all three native simulation loops.
r_sim <- emhawkes:::simulate_hawkes_r
models <- list(new("hspec", mu = 1e-308, alpha = 0, beta = 1),
  new("hspec", mu = 1e-308, alpha = matrix(0, 1, 2), beta = matrix(c(1, 2), 1, 2),
    type_col_map = list(1:2)),
  new("hspec", mu = 1e-308, alpha = 0, beta = 1, rmark = function() stop("No mark should be drawn")))
for (model in models) {
  set.seed(2); reference <- capture_conditions(r_sim(model, 5, lambda_component0 = 0)); rng <- .Random.seed
  set.seed(2); native <- capture_conditions(hsim(model, 5, lambda_component0 = 0))
  stopifnot(identical(reference$value, native$value), identical(rng, .Random.seed),
    identical(text(reference$warnings), text(native$warnings)),
    has(tail(text(native$warnings), 1), "Returned 1 of 5", "initial row"))
}
# Warnings execute R code, including handlers that draw, nest simulation or
# abort. Their draws and conditions must survive protected native unwinding.
saved_kind <- RNGkind()
for (normal in c("Inversion", "Box-Muller")) for (model in models) {
  for (mode in c("error", "draw", "handler_error")) {
    run <- function(simulate) {
      RNGkind("Mersenne-Twister", normal)
      set.seed(2)
      prior <- options(warn = if (mode == "error") 2 else 0)
      on.exit(options(prior))
      draws <- numeric()
      nested <- NULL
      value <- tryCatch(withCallingHandlers(simulate(model, 5, lambda_component0 = 0),
        warning = function(w) {
          if (mode == "error") return(invisible(NULL))
          draws <<- c(runif(2), rnorm(1))
          nested <<- simulate(poisson, 4, lambda_component0 = 0)
          gc()
          if (mode == "handler_error") stop(structure(
            list(message = "warning handler failure", call = NULL),
            class = c("warning_handler_failure", "error", "condition")))
          invokeRestart("muffleWarning")
        }), error = identity)
      following <- simulate(poisson, 4, lambda_component0 = 0)
      next_normal <- rnorm(1)
      list(value = value, draws = draws, nested = nested, following = following,
        next_normal = next_normal, rng = .Random.seed)
    }
    reference <- run(r_sim)
    native <- run(hsim)
    stopifnot(identical(reference, native))
    if (mode == "handler_error") stopifnot(inherits(native$value, "warning_handler_failure"))
  }
}
do.call(RNGkind, as.list(saved_kind))
precision <- options(digits = 3, scipen = -2)
set.seed(2); reference <- capture_conditions(r_sim(models[[1]], 100001, lambda_component0 = 0)); rng <- .Random.seed
set.seed(2); native <- capture_conditions(hsim(models[[1]], 100001, lambda_component0 = 0))
stopifnot(identical(reference$value, native$value), identical(rng, .Random.seed),
  identical(text(reference$warnings), text(native$warnings)),
  has(tail(text(native$warnings), 1), "Returned 1 of 100001"))
options(precision)
overflow <- models[[2]]
overflow@mu[] <- 1e-320
set.seed(509); reference <- capture_conditions(r_sim(overflow, 5, lambda_component0 = 0)); rng <- .Random.seed
set.seed(509); native <- capture_conditions(hsim(overflow, 5, lambda_component0 = 0))
stopifnot(inherits(reference$value, "error"), inherits(native$value, "error"),
  identical(conditionMessage(reference$value), conditionMessage(native$value)),
  has(conditionMessage(native$value), "NaN during thinning", "rescale time units"),
  identical(rng, .Random.seed), identical(text(reference$warnings), text(native$warnings)))
stops <- new("hspec", mu = function(n) matrix(if (n <= 2) .5 else 0), alpha = 0, beta = 1)
out <- capture_conditions(hsim(stops, 5, lambda_component0 = 0))
stopifnot(length(out$value$arrival) == 2L, has(text(out$warnings), "Returned 2 of 5"))

# Curvature reason codes stay machine-readable; warnings also explain them.
out <- capture_conditions(hfit(poisson, c(0, .2, .3), lambda_component0 = 0,
  fixed = c(FALSE, TRUE, FALSE), final_hessian = "stable"))
reason <- fit_diagnostics(out$value)$reason
stopifnot(reason %in% c("not_positive_definite", "rank_deficient"),
  any(grepl("information matrix", text(out$warnings), fixed = TRUE)),
  any(grepl(reason, text(out$warnings), fixed = TRUE)))
options(old)

}
