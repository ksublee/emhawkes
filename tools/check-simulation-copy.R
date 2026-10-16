# Rscript --vanilla tools/check-simulation-copy.R LIB OUT LABEL
# Save paths and RNG states for an exact cross-version comparison.
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); library(emhawkes)
source("tools/simulation-profile-cases.R"); results <- list()
record <- function(id, expr) {
  warnings <- character()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  if (inherits(value, "hreal")) value <- value[setdiff(names(value), "hspec")]
  results[[id]] <<- list(data = value, warnings = warnings,
    seed = if (exists(".Random.seed", .GlobalEnv)) .Random.seed else NULL)
}
cases <- c("univariate", "multivariate", "multi_equal", "multi_unequal", "large_equal",
  "large_unequal", "marked", "parameter_mu", "state_current", "state_dots", "state_observed", "signed")
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
  RNGkind(kind, "Inversion", "Rejection")
  for (id in cases) {
    x <- simulation_profile_case(id)
    set.seed(x$seed)
    record(paste(kind, id), hsim(x$object, size = 501, lambda_component0 = x$initial,
      N0 = seq_len(x$meta$d), Nc0 = seq_len(x$meta$d) + 2))
  }
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
h <- new("hspec", mu = c(.3, .4), alpha = matrix(.05, 2, 2), beta = matrix(1, 2, 2))
for (initial in c("automatic", "scalar", "matrix")) {
  args <- list(object = h, size = 81, N0 = 2, Nc0 = 3)
  if (initial != "automatic") args$lambda_component0 <- if (initial == "scalar") .07 else matrix(1:4 / 100, 2)
  set.seed(217); record(paste("initial", initial), do.call(hsim, args))
}
# Preserve parameter dimnames on the decayed callback state.
h@beta <- function(param = c(rate = 1)) matrix(param[1], 2, 2,
  dimnames = list(c("a", "b"), c("x", "y")))
h@impact <- function(lambda_component_n) { stopifnot(identical(dimnames(lambda_component_n), list(c("a", "b"), c("x", "y")))); matrix(0, 2, 2) }
set.seed(218); record("named_beta", hsim(h, size = 81, lambda_component0 = .03))
for (equal in c(TRUE, FALSE)) for (baseline in list(c(0, 0), c(0, .4))) {
  h <- new("hspec", mu = baseline, alpha = matrix(.02, 2, 3),
    beta = if (equal) matrix(c(1, 2), 2, 3) else matrix(1:6, 2, 3),
    type_col_map = list(c(3, 1), c(2, 3)))
  set.seed(219); record(paste("zero_mu", equal, paste(baseline, collapse = ":")),
    hsim(h, size = 81, lambda_component0 = .03))
  if (equal) {
    h <- set_residual(h, param = c(shape = 2),
      rresidual = function(n, shape) stats::rgamma(n, shape),
      dresidual = function(x, shape) stats::dgamma(x, shape),
      presidual = function(q, shape) stats::pgamma(q, shape),
      qresidual = function(p, shape) stats::qgamma(p, shape))
    set.seed(220); record(paste("gamma_residual", paste(baseline, collapse = ":")),
      hsim(h, size = 81, lambda_component0 = .03))
  }
}
# Non-exponential residuals with strictly positive baselines generate full paths.
h <- new("hspec", mu = c(.3, .4), alpha = matrix(.02, 2, 3), beta = matrix(c(1, 2), 2, 3),
  type_col_map = list(c(3, 1), c(2, 3)))
h <- set_residual(h, param = c(shape = 2),
  rresidual = function(n, shape) stats::rgamma(n, shape),
  dresidual = function(x, shape) stats::dgamma(x, shape),
  presidual = function(q, shape) stats::pgamma(q, shape),
  qresidual = function(p, shape) stats::qgamma(p, shape))
set.seed(224); record("gamma_residual_positive", hsim(h, size = 81, lambda_component0 = .03))
# Unequal signed kernels exercise validation and the thinning branch together.
h <- new("hspec", mu = c(.2, .3), alpha = matrix(c(4,-1,4,1,4,1,4,-1), 2, byrow = TRUE),
  beta = matrix(rep(c(10.56,13.44),4), 2, byrow = TRUE), type_col_map = list(1:2,3:4))
set.seed(221); record("signed_unequal", hsim(h, size = 81, lambda_component0 = 1))
# Callback closures retain context until after simulation, including ... snapshots.
for (dots in c(FALSE, TRUE)) local({
  saved <- list()
  capture <- function(label, n, mark, type, inter_arrival, N, Nc, lambda,
                      lambda_component, lambda_component_n, ...) {
    saved[[length(saved) + 1L]] <<- function() list(label, n, mark, type, inter_arrival,
      N, Nc, lambda, lambda_component, lambda_component_n)
  }
  mu <- function(n, N, ..., unused = stop("unused default")) {
    capture("mu", n = n, N = N, ...)
    matrix(c(.4, .5) + .01 * N[max(1, n - 1), ], 2)
  }
  impact <- function(n, lambda_component_n, ...) {
    capture("impact", n = n, lambda_component_n = lambda_component_n, ...)
    .02 * lambda_component_n
  }
  rmark <- function(...) { capture("rmark", ...); stats::rgeom(1, .7) + 1 }
  if (!dots) {
    for (name in c("mu", "impact", "rmark")) {
      fun <- get(name); body0 <- body(fun)
      # Keep the recorded arguments explicit while retaining the same calculation.
      formals(fun) <- alist(n = , mark = , type = , inter_arrival = , N = , Nc = ,
        lambda = , lambda_component = , lambda_component_n = , unused = stop("unused default"))
      label <- name
      body(fun) <- substitute({
        capture(LABEL, n, mark, type, inter_arrival, N, Nc, lambda, lambda_component, lambda_component_n)
        VALUE
      }, list(LABEL = label, VALUE = switch(name,
        mu = quote(matrix(c(.4, .5) + .01 * N[max(1, n - 1), ], 2)),
        impact = quote(.02 * lambda_component_n), rmark = quote(stats::rgeom(1, .7) + 1))))
      assign(name, fun)
    }
  }
  h <- new("hspec", mu = mu, alpha = matrix(.03, 2, 3), beta = matrix(1:6, 2),
    impact = impact, rmark = rmark, type_col_map = list(c(3,1),c(2,3)))
  set.seed(222); path <- hsim(h, size = 31, lambda_component0 = .03, N0 = c(2,4), Nc0 = c(5,7))
  record(paste("retained_callbacks", dots), list(path = path[setdiff(names(path), "hspec")],
    trace = lapply(saved, function(f) f())))
})
h <- new("hspec", mu = .4, alpha = .1, beta = 1)
rm(".Random.seed", envir = .GlobalEnv)
record("no_event_no_seed", hsim(h, size = 1, lambda_component0 = .03))
set.seed(223); record("no_event_existing_seed", hsim(h, size = 1, lambda_component0 = .03))
saveRDS(results, file.path(a[2], paste0(a[3], ".checks.rds")))
cat(length(results), "cross-version cases saved.\n")
