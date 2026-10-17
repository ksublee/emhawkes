# Small fixtures for GC torture and native-memory diagnostics, not benchmarks.
native_safety_cases <- function() {
  ns <- asNamespace("emhawkes")
  native <- function(name) get(name, envir = ns)
  numeric_case <- function(equal) {
    beta <- if (equal) matrix(1, 1, 2) else matrix(c(1, 2), 1)
    native(if (equal) "hawkes_simulate_equal_cpp" else "hawkes_simulate_cpp")(
      .4, matrix(.01, 1, 2), beta, matrix(.02, 1, 2), list(0:1),
      matrix(0, 5, 1), matrix(0, 5, 1), numeric(5), numeric(5), c(0, rep(1, 4)))
  }
  callbacks <- function() {
    saved <- list()
    template <- structure(matrix(.02, 2, 4), dimnames = list(c("a", "b"), letters[1:4]), unit = "rate")
    impact <- function(lambda_component_n, event) {
      saved[[length(saved) + 1L]] <<- function() lambda_component_n
      gc(FALSE)
      value <- matrix(0, 2, 4)
      columns <- list(c(2L, 1L), c(4L, 3L))[[event]]
      value[, columns] <- .01 * lambda_component_n[, columns]
      value
    }
    result <- native("hawkes_simulate_callbacks_cpp")(c(.3, .5), matrix(.01, 2, 4),
      matrix(1, 2, 4), matrix(.02, 2, 4), list(c(1L, 0L), c(3L, 2L)),
      matrix(0, 5, 2), matrix(0, 5, 2), numeric(5), numeric(5), c(0, rep(1, 4)),
      matrix(.001, 2, 4), function() { gc(FALSE); 1 + runif(1) }, impact, template)
    list(result = result, snapshots = lapply(saved, function(f) f()))
  }
  baseline <- function() {
    saved <- list()
    callback <- function(n) {
      saved[[length(saved) + 1L]] <<- function() n
      gc(FALSE)
      matrix(c(.3, .5) + n / 1000, 2)
    }
    value <- native("hawkes_baseline_cpp")(callback, list(n = 1L), 2L, 4L,
      initial = NULL, infer = TRUE, integer_n = TRUE)
    indices <- lapply(saved, function(f) f())
    stopifnot(identical(indices, as.list(seq_len(5))))
    list(value = value, n = indices)
  }
  jacobian <- function() {
    saved <- list()
    callback <- function(n) {
      saved[[length(saved) + 1L]] <<- function() n
      gc(FALSE)
      matrix(seq_len(8) / 100 + n / 1000, 4, 2,
        dimnames = list(NULL, c("level", "jump")))
    }
    value <- native("hawkes_jacobian_path_cpp")(callback, list(n = 1L), 4L,
      c(0, .1, .2, .3, .4), c(0L, 1L, 2L, 3L, 4L), c("level", "jump"))
    indices <- lapply(saved, function(f) f())
    stopifnot(identical(indices, as.list(seq_len(5))))
    list(value = value, n = indices)
  }
  likelihood <- function(observed) {
    saved <- list()
    h <- new("hspec", mu = c(.3, .5), alpha = matrix(.01, 2, 2), beta = matrix(1, 2, 2),
      impact = function(lambda_component_n) {
        saved[[length(saved) + 1L]] <<- function() lambda_component_n
        gc(FALSE)
        .01 * lambda_component_n
      })
    if (observed) h@mu <- function(n) { gc(FALSE); matrix(c(.3, .5) + n / 1000, 2) }
    value <- logLik(h, c(0, .1, .2, .3, .4), c(0, 1, 2, 1, 2),
      lambda_component0 = .02, infer = TRUE)
    list(value = value, snapshots = lapply(saved, function(f) f()))
  }
  gradient <- function() {
    native("hawkes_loglik_cpp")(c(0, .1, .2, .3), c(0L, 1L, 2L, 1L), c(0, 1, 1, 1),
      c(.3, .5), matrix(.01, 2, 2), matrix(1, 2, 2), NULL, matrix(.02, 2, 2),
      list(0L, 1L), infer = TRUE, gradient = TRUE)
  }
  warning_path <- function() {
    old <- options(warn = 0)
    on.exit(options(old))
    h <- new("hspec", mu = 1e-308, alpha = 0, beta = 1)
    nested <- NULL
    draws <- numeric()
    set.seed(2)
    value <- withCallingHandlers(hsim(h, 3, lambda_component0 = 0), warning = function(w) {
      draws <<- runif(2)
      nested <<- hsim(new("hspec", mu = .4, alpha = .01, beta = 1), 3, lambda_component0 = 0)
      gc(FALSE)
      invokeRestart("muffleWarning")
    })
    options(warn = 2)
    set.seed(2)
    reference <- tryCatch(emhawkes:::simulate_hawkes_r(h, 3, lambda_component0 = 0), error = identity)
    expected_seed <- .Random.seed
    set.seed(2)
    error <- tryCatch(hsim(h, 3, lambda_component0 = 0), error = identity)
    stopifnot(inherits(error, "error"), identical(error, reference), identical(expected_seed, .Random.seed))
    list(value = value, nested = nested, draws = draws, error = conditionMessage(error), following = runif(2))
  }
  list(simulation_equal = function() numeric_case(TRUE),
    simulation_general = function() numeric_case(FALSE), simulation_callbacks = callbacks,
    simulation_warning = warning_path,
    observed_baseline = baseline, baseline_jacobian = jacobian,
    state_likelihood = function() likelihood(FALSE),
    observed_likelihood = function() likelihood(TRUE), numeric_gradient = gradient)
}
