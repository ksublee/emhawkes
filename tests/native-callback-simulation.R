library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")
capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  list(value = value, warnings = warnings, seed = if (exists(".Random.seed", .GlobalEnv)) .Random.seed)
}
compare <- function(h, size = 31L, initial = .02, seed = 7361L) {
  set.seed(seed)
  old <- capture(r_sim(h, size, lambda_component0 = initial, N0 = 2, Nc0 = 5))
  set.seed(seed)
  new <- capture(hsim(h, size, lambda_component0 = initial, N0 = 2, Nc0 = 5))
  if (!identical(old, new)) stop("Callback simulation mismatch at size=", size,
    ", seed=", seed, ": ", paste(all.equal(old, new, tolerance = 0), collapse = "; "))
  new$value
}
models <- list()
for (d in c(1L, 2L, 4L)) for (k in unique(c(d, 2L * d))) local({
  h <- new("hspec", mu = seq(.3, .6, length.out = d),
    alpha = matrix(.03 / k, d, k), beta = matrix(rep(seq(.8, 1.7, length.out = d), k), d),
    type_col_map = if (k == d) NULL else lapply(seq_len(d), function(i) rev(c(i, i + d))))
  variants <- list(
    current = function(lambda_component_n, param = c(scale = .02), unused = stop("unused"))
      param[1] * lambda_component_n,
    hidden_rng = function(lambda_component_n) .01 * lambda_component_n + runif(1, 0, .001),
    no_args = function() matrix(0, d, k))
  for (name in names(variants)) {
    x <- h; x@impact <- variants[[name]]
    models[[paste(d, k, name)]] <<- x
  }
  h@rmark <- function(param = c(p = .7), unused = stop("unused")) rgeom(1, param[1]) + 1
  h@eta <- matrix(.002 / k, d, k)
  models[[paste(d, k, "marked")]] <<- h
  h@impact <- variants$current
  models[[paste(d, k, "combined")]] <<- h
})
models$signed <- new("hspec", mu = c(.2, .3), alpha = matrix(0, 2, 2),
  beta = matrix(1, 2, 2), impact = function(lambda_component_n) {
    matrix(rep(.1 - rowSums(lambda_component_n), 2), 2)
  })
saved_rng <- RNGkind()
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
  RNGkind(kind, "Inversion", "Rejection")
  for (h in models) for (size in c(1L, 2L, 31L, 251L)) for (seed in c(1L, 62L, 7361L))
    compare(h, size, seed = seed)
}
do.call(RNGkind, as.list(saved_rng))

# Resolved functions use the numeric loop, including unequal decays. Preserve
# evaluation order/count and the original hspec; do not rebind parameter defaults.
for (automatic in c(FALSE, TRUE)) local({
  calls <- character()
  h <- new("hspec", mu = .4, alpha = .03, beta = 1)
  h@mu <- function(param = c(level = .4)) { calls <<- c(calls, "mu"); matrix(param[1]) }
  h@alpha <- function(param = c(jump = .03)) { calls <<- c(calls, "alpha"); matrix(param[1]) }
  h@beta <- function(param = c(rate = 1)) { calls <<- c(calls, "beta"); matrix(param[1]) }
  initial <- if (automatic) NULL else .02
  set.seed(95); calls <- character()
  old <- r_sim(h, 31, lambda_component0 = initial); count <- calls; rng <- .Random.seed
  set.seed(95); calls <- character()
  new <- hsim(h, 31, lambda_component0 = initial)
  stopifnot(identical(old, new), identical(count, calls), identical(rng, .Random.seed),
    identical(calls, rep(c("mu", "alpha", "beta"), if (automatic) 2L else 1L)),
    identical(new$hspec, h))
})
h <- models[["2 4 combined"]]
h@beta <- function(param = c(rate = 1)) matrix(seq(.8, 2.2, length.out = 8), 2)
h@impact <- h@rmark <- h@eta <- NULL
compare(h)

# Dimnames and custom attributes on beta must reach each fresh decayed matrix.
local({
  saved <- list()
  h <- models[["2 4 current"]]
  dimnames(h@beta) <- list(c("a", "b"), letters[1:4]); attr(h@beta, "unit") <- "rate"
  h@impact <- function(lambda_component_n, unused = stop("unused default")) {
    stopifnot(!missing(lambda_component_n), missing(unused))
    saved[[length(saved) + 1L]] <<- function() lambda_component_n
    value <- .02 * lambda_component_n
    lambda_component_n[] <- 123 # Local modification must not change the kernel.
    value
  }
  set.seed(24); old <- r_sim(h, 31, lambda_component0 = .02)
  snapshots <- lapply(saved, function(f) f()); saved <- list(); rng <- .Random.seed
  set.seed(24); new <- hsim(h, 31, lambda_component0 = .02)
  stopifnot(identical(old, new), identical(rng, .Random.seed),
    identical(snapshots, lapply(saved, function(f) f())))
  # Keep unmodified promises too, so reusing a matrix for later events is caught.
  saved <- list()
  h@impact <- function(lambda_component_n) {
    saved[[length(saved) + 1L]] <<- function() lambda_component_n
    lambda_component_n * .02
  }
  set.seed(25); old <- r_sim(h, 31, lambda_component0 = .02)
  snapshots <- lapply(saved, function(f) f()); saved <- list(); rng <- .Random.seed
  set.seed(25); new <- hsim(h, 31, lambda_component0 = .02)
  stopifnot(identical(old, new), identical(rng, .Random.seed),
    identical(snapshots, lapply(saved, function(f) f())),
    !identical(snapshots[[1]], snapshots[[2]]))
})

# R warnings, replacement recycling and errors retain every consumed draw.
for (callback in list(function() { runif(1); stop("mark failed", call. = FALSE) },
                      function() { runif(1); warning("mark warning"); 2 },
                      function() c(2, 3), function() NULL, function() "bad")) {
  h <- models[["1 1 marked"]]; h@rmark <- callback
  compare(h)
}
for (callback in list(
  function(lambda_component_n) { runif(1); stop("impact failed", call. = FALSE) },
  function(lambda_component_n) { runif(1); warning("impact warning"); lambda_component_n * .01 },
  function(lambda_component_n) numeric(0), function(lambda_component_n) matrix("bad"))) {
  h <- models[["1 1 current"]]; h@impact <- callback
  compare(h)
}
local({
  h <- models[["1 1 current"]]
  h@impact <- function(lambda_component_n) { runif(1); warning("converted warning"); lambda_component_n }
  op <- options(warn = 2); on.exit(options(op))
  set.seed(7); old <- tryCatch(r_sim(h, 5, lambda_component0 = .02), error = conditionMessage)
  rng <- .Random.seed
  set.seed(7); new <- tryCatch(hsim(h, 5, lambda_component0 = .02), error = conditionMessage)
  stopifnot(identical(old, new), identical(rng, .Random.seed))
})

# Eta-only and parameter functions with random draws retain the R evaluation path
# during setup and exactly the same event draws afterwards.
h <- new("hspec", mu = .4, alpha = .03, beta = 1, eta = .01)
compare(h)
h@mu <- function(param = c(level = .4)) matrix(param[1] + runif(1, 0, .01))
h@alpha <- function(param = c(jump = .03)) matrix(param[1] + runif(1, 0, .001))
h@beta <- function(param = c(rate = 1)) matrix(param[1] + runif(1, 0, .01))
h@eta <- function(param = c(scale = .01)) matrix(param[1] + runif(1, 0, .001))
for (initial in list(NULL, .02)) compare(h, initial = initial)
normal_settings <- RNGkind()
for (normal in c("Inversion", "Box-Muller", "Kinderman-Ramage")) {
  RNGkind(normal.kind = normal)
  h <- models[["1 1 marked"]]
  h@rmark <- function() 1 + abs(rnorm(1))
  compare(h)
}
do.call(RNGkind, as.list(normal_settings))
h@rmark <- function() 1 + abs(.Random.seed[2] %% 3)
compare(h)

# A callback can use nested native simulations as well as ordinary R RNG calls.
# Rcpp's default nested RNG scope must not leave the outer callback on stale state.
local({
  inner <- new("hspec", mu = .3, alpha = 0, beta = 1)
  h <- models[["1 1 marked"]]
  h@rmark <- function() {
    runif(1)
    path <- hsim(inner, 3, lambda_component0 = 0)
    1 + path$arrival[3] + runif(1)
  }
  compare(h)
  h@rmark <- function() {
    runif(1); hsim(inner, 3, lambda_component0 = 0); runif(1)
    stop("nested failure", call. = FALSE)
  }
  compare(h)
})

# Initial-only paths never call event functions or create a seed.
h <- models[["1 1 combined"]]
h@rmark <- function() stop("must not call mark")
h@impact <- function(lambda_component_n) stop("must not call impact")
if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
hsim(h, 1, lambda_component0 = .02)
stopifnot(!exists(".Random.seed", .GlobalEnv))
set.seed(2); rng <- .Random.seed; hsim(h, 1, lambda_component0 = .02)
stopifnot(identical(rng, .Random.seed))

# Explicitly verify both dispatch choices; unsupported models start in R before
# any native RNG use, never by restarting a partly simulated native path.
.native_callbacks <- .native_numeric <- 0L
trace("hawkes_simulate_callbacks_cpp", where = asNamespace("emhawkes"), print = FALSE,
  tracer = quote(.GlobalEnv$.native_callbacks <- .GlobalEnv$.native_callbacks + 1L))
trace("hawkes_simulate_cpp", where = asNamespace("emhawkes"), print = FALSE,
  tracer = quote(.GlobalEnv$.native_numeric <- .GlobalEnv$.native_numeric + 1L))
trace("hawkes_simulate_equal_cpp", where = asNamespace("emhawkes"), print = FALSE,
  tracer = quote(.GlobalEnv$.native_numeric <- .GlobalEnv$.native_numeric + 1L))
compare(models[["2 4 combined"]])
h <- new("hspec", mu = function(param = c(level = .4)) matrix(param[1]), alpha = .03, beta = 1)
compare(h)
stopifnot(.native_callbacks == 1L, .native_numeric == 1L)
unsupported <- list()
h <- models[["2 4 combined"]]
x <- h; x@beta[,] <- matrix(1:8, 2); unsupported[[1]] <- x
x <- h; x@mu <- function(n) matrix(c(.3, .6)); unsupported[[2]] <- x
x <- h; x@impact <- function(lambda_component_n, ...) lambda_component_n * .02; unsupported[[3]] <- x
x <- h; x@rmark <- function(n) 1; unsupported[[4]] <- x
x <- h; x@mu <- matrix(c(0, .3)); unsupported[[5]] <- x
x <- h; x@impact <- function(N, lambda_component_n) lambda_component_n * .02; unsupported[[6]] <- x
Math.simulation_rate <- function(x, ...) NextMethod()
Ops.simulation_rate <- function(e1, e2) NextMethod()
x <- h; x@beta <- function(param = c(rate = 1)) {
  structure(matrix(param[1], 2, 4), class = c("simulation_rate", "matrix", "array"))
}; unsupported[[7]] <- x
for (x in unsupported) compare(x)
stopifnot(.native_callbacks == 1L, .native_numeric == 1L)
untrace("hawkes_simulate_callbacks_cpp", where = asNamespace("emhawkes"))
untrace("hawkes_simulate_cpp", where = asNamespace("emhawkes"))
untrace("hawkes_simulate_equal_cpp", where = asNamespace("emhawkes"))

# Deterministic callbacks permit independent R/C++ likelihood reconstruction.
for (id in c("1 1 current", "2 4 marked", "4 8 combined")) {
  h <- models[[id]]; path <- compare(h, 301); rng <- .Random.seed
  for (engine in c("r", "cpp")) {
    op <- options(emhawkes.engine = engine)
    replay <- logLik(h, path$inter_arrival, path$type, path$mark,
      lambda_component0 = .02, N0 = 2, Nc0 = 5, infer = TRUE)
    options(op)
    for (name in c("lambda", "rambda", "lambda_component", "rambda_component"))
      stopifnot(identical(as.numeric(path[[name]]), as.numeric(replay[[name]])))
    stopifnot(identical(rng, .Random.seed))
  }
}
