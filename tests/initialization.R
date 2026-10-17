library(emhawkes)

capture_initial <- function(expr) {
  messages <- character(0)
  value <- withCallingHandlers(force(expr), message = function(m) {
    messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  list(value = value, messages = messages)
}
fails_initial <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
close_initial <- function(a, b, tolerance = 1e-10)
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance)))
old_initial_options <- options(emhawkes.engine = "auto")
u <- new("hspec", mu = 1, alpha = .4, beta = 1)
set.seed(814)
rng <- .Random.seed
one <- capture_initial(hsim(u, 1))
stopifnot(identical(rng, .Random.seed), length(one$messages) == 1L,
          grepl("auto -> mean", one$messages))
close_initial(one$value$lambda_component[1, ], 2 / 3)
close_initial(one$value$lambda[1, ], 5 / 3)
close_initial(capture_initial(hsim(u, 1, lambda_component0 = NULL))$value$lambda, 5 / 3)

# Verify stationary balance independently, including overlapping mappings.
models <- list(u,
  new("hspec", mu = c(.7, 1), alpha = matrix(c(.12, .03, .04, .1), 2), beta = matrix(1.5, 2, 2)),
  new("hspec", mu = c(.7, 1), alpha = matrix(seq(.01, .08, .01), 2, 4), beta = matrix(1.5, 2, 4),
      type_col_map = list(c(1, 2, 3), c(2, 4))))
for (h in models) {
  out <- capture_initial(hsim(h, 1))$value
  components <- matrix(out$lambda_component[1, ], h@dimens, byrow = TRUE)
  rates <- as.vector(out$lambda[1, ])
  mapping <- if (is.null(h@type_col_map)) as.list(seq_len(h@dimens)) else h@type_col_map
  expected <- matrix(0, nrow(h@beta), ncol(h@beta))
  for (j in seq_len(h@dimens)) expected[, mapping[[j]]] <-
    expected[, mapping[[j]], drop = FALSE] + h@alpha[, mapping[[j]], drop = FALSE] / h@beta[, mapping[[j]], drop = FALSE] * rates[j]
  close_initial(components, expected)
  for (seed in c(1, 29, 704)) {
    set.seed(seed); auto <- capture_initial(hsim(h, 41))$value; rng <- .Random.seed
    set.seed(seed); explicit <- hsim(h, 41, lambda_component0 = components)
    stopifnot(identical(auto, explicit), identical(rng, .Random.seed))
    set.seed(seed); reference <- capture_initial(emhawkes:::simulate_hawkes_r(h, 41))$value
    stopifnot(identical(auto, reference), identical(rng, .Random.seed))
    for (engine in c("r", "cpp")) {
      options(emhawkes.engine = engine)
      replay <- capture_initial(logLik(h, auto$inter_arrival, auto$type, auto$mark, infer = TRUE))$value
      stopifnot(identical(auto$lambda, replay$lambda), identical(auto$rambda, replay$rambda))
    }
  }
}

# Unknown means and unstable starts select zero once, without drawing RNG.
unknown <- new("hspec", mu = function(n) matrix(.8), alpha = .2, beta = 1)
unstable <- new("hspec", mu = 1, alpha = 1, beta = 1)
for (h in list(unknown, unstable, new("hspec", mu = 1, alpha = .2, beta = 1, eta = .1))) {
  set.seed(14); rng <- .Random.seed
  out <- capture_initial(hsim(h, 1))
  stopifnot(identical(rng, .Random.seed), all(out$value$lambda_component == 0),
    length(out$messages) == 1L, grepl("auto -> zero", out$messages),
    grepl("explicitly is recommended", out$messages))
  fails_initial(hsim(h, 1, lambda_component0 = "mean"), "Mean initialization is unavailable")
}
zero <- capture_initial(hsim(u, 1, lambda_component0 = "zero"))
stopifnot(all(zero$value$lambda_component == 0), grepl("components: zero", zero$messages))
fails_initial(hsim(u, 1, lambda_component0 = "invalid"), "lambda_component0")

# Selecting a zero state must not evaluate beta an extra time.
beta_calls <- 0L
callback_beta <- u
callback_beta@beta <- function(param = c(b = 1)) {
  beta_calls <<- beta_calls + 1L
  matrix(param[1])
}
for (engine in c("r", "auto", "cpp")) {
  options(emhawkes.engine = engine)
  for (initial in list("auto", 0)) {
    beta_calls <- 0L
    suppressMessages(logLik(callback_beta, c(0, .2, .4), lambda_component0 = initial))
    stopifnot(beta_calls == 1L)
  }
}

# Mean initial-state score includes the full chain rule for mapped kernels.
for (h in models) {
  obs <- emhawkes:::normalize_observation(h, c(0, .2, .4, .1, .3),
    type = c(0, rep(seq_len(h@dimens), length.out = 4)))
  map <- emhawkes:::build_parameter_map(h, FALSE)
  objective <- emhawkes:::make_hawkes_objective(h, obs, NULL, map)
  p <- map$start
  score <- objective$gradient(p)
  for (step in c(1e-5, 3e-6)) {
    numeric <- vapply(seq_along(p), function(j) {
      above <- below <- p; above[j] <- above[j] + step; below[j] <- below[j] - step
      (objective$loglik(above) - objective$loglik(below)) / (2 * step)
    }, numeric(1))
    close_initial(score, numeric, 2e-7)
  }
}

# Selection is fixed for a fit; saved mean replay is parameter dependent.
large <- new("hspec", mu = 1, alpha = .9, beta = 1)
large_map <- emhawkes:::build_parameter_map(large, FALSE)
large_objective <- emhawkes:::make_hawkes_objective(large,
  emhawkes:::normalize_observation(large, c(0, .2, .4)), NULL, large_map)
proposal <- large_map$start; proposal["mu1"] <- 1e308
stopifnot(identical(as.numeric(large_objective$loglik(proposal)), -Inf),
  all(is.na(large_objective$gradient(proposal))))
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  fit <- capture_initial(hfit(u, c(0, .2, .4, .1), optimizer = "stats", fixed = rep(TRUE, 3)))
  stopifnot(sum(grepl("Initial intensity components", fit$messages)) == 1L,
    identical(fit_initial_state(fit$value)$lambda_component0, "mean"))
  p <- coef(fit$value); p["alpha1"] <- p["beta1"]
  stopifnot(identical(as.numeric(fit$value$objectiveFn(p)), -Inf))
  replay <- suppressMessages(do.call(logLik, c(list(object = fitted_hspec(fit$value),
    inter_arrival = c(0, .2, .4, .1)), fit_initial_state(fit$value))))
  close_initial(replay, logLik(fit$value))
  # Pre-policy metadata used the legacy half-mean rule. Recover its actual state.
  legacy <- fit$value
  metadata <- attr(legacy, "emhawkes")
  metadata$initialization <- NULL
  metadata$initial_state$lambda_component0 <- NULL
  attr(legacy, "emhawkes") <- metadata
  restored <- fit_initial_state(legacy)
  close_initial(restored$lambda_component0, 1 / 3)
  legacy$maximum <- as.numeric(suppressMessages(logLik(u, c(0, .2, .4, .1), lambda_component0 = 1 / 3)))
  legacy_replay <- suppressMessages(do.call(logLik, c(list(object = fitted_hspec(legacy),
    inter_arrival = c(0, .2, .4, .1)), fit_initial_state(legacy))))
  close_initial(legacy_replay, logLik(legacy), 0)
}
options(old_initial_options)

# Numerical constrained optimization must never probe an invalid negative rate.
near <- new("hspec", mu = 1e-7, alpha = 0, beta = 1)
strict <- function(object) {
  if (object@mu[1] <= 0) stop("outside positive domain")
  -(object@mu[1] - .2)^2
}
near_fit <- hfit(near, optimizer = "stats", mylogLik = strict,
  fixed = c(FALSE, TRUE, TRUE),
  constraints = list(ineqA = matrix(c(1, 0, 0), 1), ineqB = 0),
  control = list(maxit = 300))
close_initial(coef(near_fit)["mu1"], .2, 1e-5)
