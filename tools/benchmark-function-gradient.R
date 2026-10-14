# Install the package, then run from the repository root:
# Rscript tools/benchmark-function-gradient.R output-directory [1001,5001] [3] [all|case] [compare|snapshot]
# snapshot works on the previous release and saves default-fit fields for comparison.
library(emhawkes)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else file.path("tools", "function-gradient-results")
sizes <- if (length(args) > 1) as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]]) else c(1001L, 5001L)
rounds <- if (length(args) > 2) as.integer(args[3]) else 3L
selection <- if (length(args) > 3) args[4] else "all"
snapshot <- length(args) > 4 && args[5] == "snapshot"
stopifnot(length(sizes) > 0, all(is.finite(sizes)), all(sizes >= 2), rounds >= 3)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "auto")

mu_log <- function(param = c(log_mu = log(.3))) matrix(exp(param))
alpha_log <- function(param = c(log_alpha = log(.4))) matrix(exp(param))
beta_log <- function(param = c(log_beta = log(1.2))) matrix(exp(param))
exp_jac <- function(param) matrix(exp(param), dimnames = list(NULL, names(param)))
alpha_zero <- function(param = c(log_a = log(.2), log_b = log(.08))) {
  matrix(c(exp(param["log_a"]), 0, exp(param["log_b"]), exp(param["log_a"])), 2)
}
alpha_zero_jac <- function(param) cbind(log_a = c(1, 0, 0, 1) * exp(param["log_a"]),
                                      log_b = c(0, 0, 1, 0) * exp(param["log_b"]))
shared_mu <- function(param = c(scale = log(.25))) matrix(c(1, 1.25) * exp(param), 2)
shared_beta <- function(param = c(scale = log(.25))) matrix(1 + exp(param), 2, 2)
shared_alpha <- function(param = c(a = log(.2), b = log(.06))) {
  matrix(exp(param[c("a", "b", "b", "a")]), 2)
}
shared_jac <- list(mu = function(param) matrix(c(1, 1.25) * exp(param), dimnames = list(NULL, "scale")),
  beta = function(param) matrix(rep(exp(param), 4), dimnames = list(NULL, "scale")),
  alpha = function(param) cbind(a = c(1, 0, 0, 1) * exp(param["a"]),
                              b = c(0, 1, 1, 0) * exp(param["b"])))
multi_alpha <- function(param = c(slow = log(.2), fast = log(.9))) matrix(exp(param), 1)
multi_beta <- function(param = numeric(0)) matrix(c(.8, 5), 1)
cases <- list(
  univariate_log = list(seed = 8741L, initial = NULL,
    truth = new("hspec", mu = .35, alpha = .5, beta = 1.5),
    start = new("hspec", mu = mu_log, alpha = alpha_log, beta = beta_log),
    jacobian = list(mu = exp_jac, alpha = exp_jac, beta = exp_jac)),
  bivariate_fixed_zero = list(seed = 8742L, initial = .1,
    truth = new("hspec", mu = c(.25, .35), alpha = matrix(c(.3, 0, .15, .3), 2), beta = matrix(1.2, 2, 2)),
    start = new("hspec", mu = c(.2, .3), alpha = alpha_zero,
      beta = function(param = c(log_decay = 0)) matrix(exp(param), 2, 2)),
    jacobian = list(alpha = alpha_zero_jac,
      beta = function(param) matrix(rep(exp(param), 4), dimnames = list(NULL, "log_decay")))),
  bivariate_shared_name = list(seed = 8743L, initial = NULL, log_bounds = TRUE,
    truth = new("hspec", mu = c(.3, .375), alpha = matrix(c(.25, .08, .08, .25), 2), beta = matrix(1.3, 2, 2)),
    start = new("hspec", mu = shared_mu, alpha = shared_alpha, beta = shared_beta),
    jacobian = shared_jac),
  multi_kernel_fixed_decay = list(seed = 8744L, initial = .1,
    truth = new("hspec", mu = .3, alpha = matrix(c(.25, 1.1), 1),
      beta = matrix(c(.8, 5), 1), type_col_map = list(1:2)),
    start = new("hspec", mu = .25, alpha = multi_alpha, beta = multi_beta, type_col_map = list(1:2)),
    jacobian = list(alpha = function(param) {
      value <- diag(exp(param)); colnames(value) <- names(param); value
    })))
if (selection != "all") {
  stopifnot(selection %in% names(cases))
  cases <- cases[selection]
}
timings <- comparisons <- estimates <- list()
snapshots <- list()
writeLines(trimws(capture.output(sessionInfo(),
  Sys.info()[c("sysname", "release", "version", "machine")]), which = "right"),
  file.path(output, "session-info.txt"))
for (size in sizes) for (case in names(cases)) {
  spec <- cases[[case]]
  set.seed(spec$seed)
  path <- hsim(spec$truth, size = size, lambda_component0 = .1)
  map <- emhawkes:::build_parameter_map(spec$start, TRUE)
  # Only the untransformed numeric baseline requires a positivity constraint.
  indices <- if (map$components$mu$kind == "numeric") unique(map$components$mu$index) else integer()
  constraints <- if (length(indices)) list(ineqA = diag(length(map$start))[indices, , drop = FALSE],
                                          ineqB = rep(-1e-5, length(indices))) else NULL
  # Broad common bounds prevent exp-scale proposals from overflowing the
  # existing R automatic initializer in the shared-name model.
  if (isTRUE(spec$log_bounds)) constraints <- list(
    ineqA = rbind(diag(length(map$start)), -diag(length(map$start))),
    ineqB = c(rep(8, length(map$start)), rep(2, length(map$start))))
  fit_once <- function(mode, stable = FALSE) {
    extra <- if (mode == "analytic") list(gradient = "analytic", jacobian = spec$jacobian) else list()
    suppressMessages(do.call(hfit, c(list(object = spec$start, inter_arrival = path$inter_arrival,
      type = path$type, lambda_component0 = spec$initial, method = "BFGS", constraints = constraints,
      final_hessian = if (stable) "stable" else "maxLik",
      control = list(iterlim = 500, reltol = 1e-10)), extra)))
  }
  fits <- list(numeric = fit_once("numeric"))
  fields <- c("estimate", "maximum", "gradient", "hessian", "code", "iterations")
  snapshots[[paste(case, size, sep = ":")]] <- fits$numeric[fields]
  saveRDS(snapshots, file.path(output, "default-fits.rds"))
  if (snapshot) next
  fits$analytic <- fit_once("analytic")
  delta_ll <- abs(as.numeric(logLik(fits$numeric)) - as.numeric(logLik(fits$analytic)))
  delta_coef <- max(abs(coef(fits$numeric) - coef(fits$analytic)))
  # Check the original R likelihood at the analytic estimate.
  original <- options(emhawkes.engine = "r")
  reference_ll <- suppressMessages(do.call(logLik, c(list(object = fitted_hspec(fits$analytic),
    inter_arrival = path$inter_arrival, type = path$type), fit_initial_state(fits$analytic))))
  options(original)
  reference_delta <- abs(as.numeric(reference_ll) - as.numeric(logLik(fits$analytic)))
  stopifnot(fits$numeric$code == 0, fits$analytic$code == 0,
            delta_ll < 1e-5, delta_coef < 1e-3, reference_delta < 1e-8)
  timed_batch <- function(mode, count) unname(system.time(for (i in seq_len(count)) fit_once(mode))[["elapsed"]])
  counts <- vapply(names(fits), function(mode) {
    n <- 1L
    repeat {
      if (timed_batch(mode, n) >= .2 || n >= 128L) return(n)
      n <- 2L * n
    }
  }, integer(1))
  seconds <- matrix(NA_real_, rounds, 2, dimnames = list(NULL, names(fits)))
  for (round in seq_len(rounds)) for (mode in if (round %% 2) names(fits) else rev(names(fits))) {
    gc()
    elapsed <- timed_batch(mode, counts[[mode]])
    seconds[round, mode] <- elapsed / counts[[mode]]
    timings[[length(timings) + 1L]] <- data.frame(case = case, rows = size, round = round,
      mode = mode, fits_in_batch = counts[[mode]], batch_seconds = elapsed, seconds = seconds[round, mode])
  }
  # The default Hessian differs by derivative path. Report it without assuming
  # identical uncertainty, and compare common optional stabilization separately.
  se <- lapply(fits, function(x) suppressWarnings(sqrt(diag(vcov(x)))))
  stable <- lapply(names(fits), function(mode) fit_once(mode, stable = TRUE))
  stable_se <- lapply(stable, function(x) suppressWarnings(sqrt(diag(vcov(x)))))
  se_delta <- if (all(is.finite(unlist(stable_se))))
    max(abs(stable_se[[1]] - stable_se[[2]]) / pmax(stable_se[[1]], 1e-12)) else NA_real_
  row <- data.frame(case = case, rows = size, parameters = length(map$start),
    numeric_seconds = median(seconds[, "numeric"]), analytic_seconds = median(seconds[, "analytic"]),
    speedup = median(seconds[, "numeric"]) / median(seconds[, "analytic"]),
    numeric_code = fits$numeric$code, analytic_code = fits$analytic$code,
    loglik_abs_difference = delta_ll, coefficient_max_abs_difference = delta_coef,
    r_loglik_abs_difference = reference_delta, stable_se_max_relative_difference = se_delta)
  comparisons[[length(comparisons) + 1L]] <- row
  estimates[[length(estimates) + 1L]] <- data.frame(case = case, rows = size,
    parameter = names(coef(fits$numeric)), numeric = coef(fits$numeric), analytic = coef(fits$analytic),
    numeric_se = se$numeric, analytic_se = se$analytic,
    stable_numeric_se = stable_se[[1]], stable_analytic_se = stable_se[[2]], row.names = NULL)
  print(row, row.names = FALSE)
  for (entry in c("timings", "comparisons", "estimates"))
    write.csv(do.call(rbind, get(entry)), file.path(output, paste0(entry, ".csv")), row.names = FALSE)
}
cat("Results:", normalizePath(output, winslash = "/"), "\n")
