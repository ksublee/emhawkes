library(emhawkes)

# Retained callback promises must describe the call, including full histories.
# Verify both explicit arguments and the legacy ... context.
for (dots in c(FALSE, TRUE)) local({
  saved <- list()
  record <- function(label, n, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n) {
    saved[[length(saved) + 1L]] <<- function() list(label = label, n = n,
      mark = mark, type = type, inter_arrival = inter_arrival, N = N, Nc = Nc,
      lambda = lambda, lambda_component = lambda_component,
      lambda_component_n = lambda_component_n)
  }
  mu <- function(n, N, mark, type, inter_arrival, Nc, lambda,
                 lambda_component, lambda_component_n, unused = stop("unused default")) {
    record("mu", n, mark, type, inter_arrival, N, Nc, lambda,
           lambda_component, lambda_component_n)
    matrix(c(.4, .5) + .01 * N[max(1, n - 1), ], 2)
  }
  impact <- function(n, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n) {
    record("impact", n, mark, type, inter_arrival, N, Nc, lambda,
           lambda_component, lambda_component_n)
    .02 * lambda_component_n
  }
  rmark <- function(n, mark, type, inter_arrival, N, Nc, lambda,
                    lambda_component, lambda_component_n) {
    record("rmark", n, mark, type, inter_arrival, N, Nc, lambda,
           lambda_component, lambda_component_n)
    stats::rgeom(1, .7) + 1
  }
  if (dots) {
    narrow_mu <- mu; narrow_impact <- impact; narrow_rmark <- rmark
    mu <- function(...) narrow_mu(...)
    impact <- function(...) narrow_impact(...)
    rmark <- function(...) narrow_rmark(...)
    # The wrappers receive all supported context, including model parameters.
    for (name in c("narrow_mu", "narrow_impact", "narrow_rmark")) {
      fun <- get(name); formals(fun) <- c(formals(fun), alist(... = ))
      assign(name, fun)
    }
  }
  initial <- matrix(seq(.01, .06, length.out = 6), 2, 3)
  h <- new("hspec", mu = mu, alpha = matrix(.03, 2, 3),
    beta = matrix(c(.8, 1.1, 1.7, 2.2, 3.1, 4), 2, 3),
    impact = impact, rmark = rmark, type_col_map = list(c(3, 1), c(2, 3)))
  set.seed(119)
  path <- hsim(h, size = 12, lambda_component0 = initial, N0 = c(2, 4), Nc0 = c(5, 7))
  seed <- .Random.seed
  traces <- lapply(saved, function(f) f())
  stopifnot(identical(seed, .Random.seed), length(traces) == 2L + 3L * 11L)
  stopifnot(identical(vapply(traces, `[[`, "", "label"),
    c("mu", "mu", rep(c("rmark", "impact", "mu"), 11))))
  prefix <- function(value, through, mark = FALSE) {
    if (is.null(dim(value))) {
      result <- if (mark) c(0, rep(1, length(value) - 1)) else numeric(length(value))
      if (through) result[seq_len(through)] <- value[seq_len(through)]
    } else {
      result <- value; result[] <- 0
      if (through) result[seq_len(through), ] <- value[seq_len(through), , drop = FALSE]
    }
    result
  }
  for (i in seq_along(traces)) {
    z <- traces[[i]]; n <- z$n; event <- z$label != "mu"
    rows <- if (event) n else max(1, n - 1)
    marked_rows <- if (z$label == "rmark") n - 1 else rows
    intensity_rows <- if (z$label == "rmark") n - 1 else rows
    stopifnot(identical(z$type, prefix(path$type, if (i == 1) 0 else rows)),
      identical(z$inter_arrival, prefix(path$inter_arrival, if (i == 1) 0 else rows)),
      identical(z$mark, prefix(path$mark, if (i == 1) 0 else marked_rows, mark = TRUE)),
      identical(z$N, prefix(path$N, rows)), identical(z$Nc, prefix(path$Nc, marked_rows)))
    if (i == 1) {
      stopifnot(is.null(z$lambda), is.null(z$lambda_component), is.null(z$lambda_component_n))
    } else {
      stopifnot(identical(z$lambda, prefix(path$lambda, intensity_rows)),
        identical(z$lambda_component, prefix(path$lambda_component, intensity_rows)))
      current_row <- if (z$label == "rmark") n - 1 else rows
      current <- if (current_row == 1) initial else
        matrix(path$lambda_component[current_row, ], 2, byrow = TRUE)
      stopifnot(identical(z$lambda_component_n, current))
    }
  }
  # Replay the same model without recording callbacks that request full history.
  replay_h <- new("hspec", mu = function(n, N) {
    matrix(c(.4, .5) + .01 * N[max(1, n - 1), ], 2)
  }, alpha = h@alpha, beta = h@beta, type_col_map = h@type_col_map,
    impact = function(lambda_component_n) .02 * lambda_component_n)
  for (engine in c("r", "cpp")) {
    previous <- options(emhawkes.engine = engine)
    replay <- logLik(replay_h, path$inter_arrival, path$type, path$mark, N0 = c(2, 4),
      Nc0 = c(5, 7), lambda_component0 = initial, infer = TRUE)
    options(previous)
    for (field in c("lambda", "lambda_component", "rambda", "rambda_component")) {
      stopifnot(isTRUE(all.equal(as.numeric(replay[[field]]), as.numeric(path[[field]]),
        tolerance = 1e-12)))
    }
    stopifnot(identical(seed, .Random.seed))
  }
})

# Small-state arrival inputs preserve the history-input RNG contract.
for (helper in c("rarrival", "darrival")) for (equal in c(TRUE, FALSE)) {
  if (helper == "darrival" && !equal) next
  fun <- getFromNamespace(helper, "emhawkes")
  current <- matrix(c(.1, .2, .3, .4), 2)
  beta <- if (equal) matrix(c(1, 2), 2, 2) else matrix(c(1, 2, 3, 4), 2)
  args <- list(n = 3, mu = c(.3, .4), beta = beta, dimens = 2,
    all_rows_beta_equal = equal, rresidual = function(n) stats::rgamma(n, 2))
  history <- rbind(0, as.vector(t(current)), 0)
  set.seed(62); before <- do.call(fun, c(args, list(rambda_component = history))); seed <- .Random.seed
  set.seed(62); after <- do.call(fun, c(args, list(current_component = current)))
  stopifnot(identical(before, after), identical(seed, .Random.seed))
}

# No-event requests leave RNG untouched. Extinction still trims every output.
for (equal in c(TRUE, FALSE)) {
  h <- new("hspec", mu = 0, alpha = matrix(0, 1, 2),
    beta = matrix(if (equal) c(1, 1) else c(1, 2), 1), type_col_map = list(1:2))
  for (residual in c(FALSE, TRUE)) {
    if (residual && !equal) next
    if (residual) h <- set_residual(h, param = c(shape = 2),
      rresidual = function(n, shape) stats::rgamma(n, shape),
      dresidual = function(x, shape) stats::dgamma(x, shape),
      presidual = function(q, shape) stats::pgamma(q, shape),
      qresidual = function(p, shape) stats::qgamma(p, shape))
    set.seed(73); seed <- .Random.seed
    one <- hsim(h, size = 1, lambda_component0 = 0)
    stopifnot(identical(seed, .Random.seed), length(one$type) == 1)
    warnings <- character()
    extinct <- withCallingHandlers(hsim(h, size = 5, lambda_component0 = 0), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    stopifnot(length(warnings) == 1, grepl("no further event", warnings, fixed = TRUE),
      identical(one, extinct))
  }
}
