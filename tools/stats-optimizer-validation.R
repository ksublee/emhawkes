# Run against an installed development package. No package sources are changed.
library(emhawkes)
close <- function(a, b, tolerance) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance)))
rows <- list()
old <- options(emhawkes.engine = "auto")
for (seed in 1:4) for (mapped in c(FALSE, TRUE)) for (reduced in c(FALSE, TRUE)) {
  set.seed(8100 + seed)
  kernels <- if (mapped) 3 else 2
  mu <- function(param = c(log_mu1 = log(.4), log_mu2 = log(.6))) matrix(exp(param), 2)
  alpha <- matrix(rep(c(.04, .03), kernels), 2)
  beta <- matrix(rep(c(1.1, 1.4), kernels), 2)
  model <- new("hspec", mu = mu, alpha = alpha, beta = beta,
    eta = if (seed %% 2L == 0L) matrix(.01, 2, kernels) else NULL,
    rmark = function() 2,
    type_col_map = if (mapped) list(c(1,3), 2) else NULL)
  path <- hsim(model, size = 121 + seed * 10, lambda_component0 = .05)
  map <- emhawkes:::build_parameter_map(model, reduced)
  fixed <- !startsWith(names(map$start), "log_mu")
  jacobian <- list(mu = function(param) {
    answer <- diag(exp(param)); colnames(answer) <- names(param); answer
  })
  fits <- list()
  for (engine in c("r", "cpp")) for (derivative in c("numeric", "analytic")) {
    if (engine == "r" && derivative == "analytic") next
    options(emhawkes.engine = engine)
    args <- list(object = model, inter_arrival = path$inter_arrival, type = path$type,
      mark = path$mark,
      lambda_component0 = .05, reduced = reduced, fixed = fixed,
      gradient = derivative, jacobian = jacobian,
      control = list(maxit = 300, reltol = 1e-12), optimizer = "stats")
    fit <- suppressMessages(do.call(hfit, args))
    stopifnot(fit$code == 0L, fit_diagnostics(fit)$status == "ok")
    estimated <- fitted_hspec(fit)
    likelihoods <- vapply(c("r", "cpp"), function(replay_engine) {
      options(emhawkes.engine = replay_engine)
      as.numeric(logLik(estimated, path$inter_arrival, path$type, mark = path$mark, lambda_component0 = .05))
    }, 0)
    close(likelihoods, rep(as.numeric(logLik(fit)), 2), 1e-10)
    score <- fit_diagnostics(fit)$gradient
    stopifnot(max(abs(score)) < 1e-3)
    fits[[paste(engine, derivative)]] <- fit
    if (requireNamespace("maxLik", quietly = TRUE)) {
      options(emhawkes.engine = engine)
      args$optimizer <- "maxLik"
      args$control <- list(iterlim = 300, reltol = 1e-12)
      args$final_hessian <- "stable"
      reference <- suppressMessages(do.call(hfit, args))
      close(coef(fit), coef(reference), 2e-5)
      close(logLik(fit), logLik(reference), 1e-9)
      close(vcov(fit), vcov(reference), 3e-5)
    }
    rows[[length(rows) + 1L]] <- data.frame(seed = seed, mapped = mapped,
      reduced = reduced, engine = engine, gradient = derivative,
      code = fit$code, replay_error = max(abs(likelihoods - as.numeric(logLik(fit)))))
  }
  close(coef(fits[[1]]), coef(fits[[2]]), 2e-5)
  close(coef(fits[[2]]), coef(fits[[3]]), 2e-5)
}
options(old)
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
cat("Validated", nrow(result), "multivariate fits; max replay error:", max(result$replay_error), "\n")

# Poison maxLik entry points to catch accidental use by the stats provider.
if (requireNamespace("maxLik", quietly = TRUE)) {
  functions <- c("maxLik", "activePar", "objectiveFn", "hessian", "gradient", "returnCode", "returnMessage")
  for (name in functions) trace(name, quote(stop("Unexpected maxLik API call")),
    where = asNamespace("maxLik"), print = FALSE)
  h <- new("hspec", mu = .5, alpha = 0, beta = 1)
  fit <- suppressMessages(hfit(h, c(0,.2,.3,.5), optimizer = "stats",
    lambda_component0 = 0, fixed = c(FALSE,TRUE,TRUE)))
  stopifnot(fit_diagnostics(fit)$status == "ok")
  invisible(summary(fit)); invisible(fitted_hspec(fit))
  for (name in functions) untrace(name, where = asNamespace("maxLik"))
  cat("Stats path with poisoned maxLik entry points: OK\n")
}
