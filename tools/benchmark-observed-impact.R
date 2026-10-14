# Rscript tools/benchmark-observed-impact.R output-directory [501,2001] [3] [all|case] [compare|snapshot]
library(emhawkes)
source("tools/observed-impact-cases.R")
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else "tools/observed-impact-results"
sizes <- if (length(args) > 1) as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]]) else c(501L, 2001L)
rounds <- if (length(args) > 2) as.integer(args[3]) else 3L
selection <- if (length(args) > 3) args[4] else "all"
snapshot <- length(args) > 4 && args[5] == "snapshot"
stopifnot(length(sizes) > 0, all(is.finite(sizes)), all(sizes >= 2), rounds >= 3)
cases <- observed_impact_cases()
if (selection != "all") { stopifnot(selection %in% names(cases)); cases <- cases[selection] }
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "auto")
writeLines(trimws(capture.output(sessionInfo(), Sys.info()[c("sysname", "release", "version", "machine")]), which = "right"), file.path(output, "session-info.txt"))
timings <- comparisons <- estimates <- list()
default_fits <- list()
for (size in sizes) for (case in names(cases)) {
  spec <- cases[[case]]
  set.seed(spec$seed)
  path <- hsim(spec$truth, size = size, lambda_component0 = .1, N0 = 2, Nc0 = 3)
  map <- emhawkes:::build_parameter_map(spec$start, TRUE)
  obs <- emhawkes:::normalize_observation(spec$start, path$inter_arrival, path$type, path$mark, N0 = 2, Nc0 = 3)
  constraints <- list(ineqA = diag(length(map$start)), ineqB = rep(-1e-5, length(map$start)))
  control <- list(iterlim = 500, reltol = 1e-10)
  public_fit <- function(mode = "numeric", stable = FALSE) suppressMessages(hfit(spec$start,
    path$inter_arrival, path$type, path$mark, lambda_component0 = spec$initial,
    N0 = 2, Nc0 = 3, constraints = constraints, control = control, gradient = mode,
    jacobian = if (mode == "analytic") spec$jacobian else NULL,
    final_hessian = if (stable) "stable" else "maxLik"))
  default <- public_fit()
  default_fits[[paste(case, size, sep = ":")]] <- default[c("estimate", "maximum", "gradient", "hessian", "code", "iterations")]
  saveRDS(default_fits, file.path(output, "default-fits.rds"))
  if (snapshot) next
  # Time the same optimizer separately from objective setup and metadata.
  # R uses the unmodified hfit objective; the other two use identical R callback
  # evaluation and the C++ kernel, with or without derivative calculation.
  r_objective <- maxLik::objectiveFn(default)
  stopifnot(is.function(r_objective))
  cpp_numeric <- emhawkes:::make_hawkes_objective(spec$start, obs, spec$initial, map, derivatives = FALSE)
  cpp_analytic <- emhawkes:::make_hawkes_objective(spec$start, obs, spec$initial, map, jacobian = spec$jacobian)
  objectives <- list(r_numeric = list(loglik = r_objective, gradient = NULL),
                     cpp_numeric = cpp_numeric, cpp_analytic = cpp_analytic)
  fit_once <- function(mode) {
    obj <- objectives[[mode]]
    maxLik::maxLik(logLik = obj$loglik, grad = obj$gradient, start = map$start,
      method = "BFGS", constraints = constraints, control = control)
  }
  fits <- lapply(names(objectives), fit_once)
  names(fits) <- names(objectives)
  public <- public_fit("analytic")
  stopifnot(max(abs(coef(public) - coef(fits$cpp_analytic))) < 1e-10,
            identical(coef(default), coef(fits$r_numeric)))
  reference <- suppressMessages(do.call(logLik, c(list(object = fitted_hspec(public),
    inter_arrival = path$inter_arrival, type = path$type, mark = path$mark), fit_initial_state(public))))
  ref_delta <- abs(as.numeric(reference) - as.numeric(logLik(public)))
  stopifnot(ref_delta < 1e-8)
  for (fit in fits) stopifnot(fit$code == 0,
    abs(as.numeric(logLik(default)) - as.numeric(logLik(fit))) < 1e-5,
    max(abs(coef(default) - coef(fit))) < 1e-3)
  timed_batch <- function(mode, count) unname(system.time(for (i in seq_len(count)) fit_once(mode))[["elapsed"]])
  counts <- vapply(names(fits), function(mode) {
    count <- 1L
    repeat {
      if (timed_batch(mode, count) >= .2 || count >= 128L) return(count)
      count <- count * 2L
    }
  }, integer(1))
  seconds <- matrix(NA_real_, rounds, length(fits), dimnames = list(NULL, names(fits)))
  for (round in seq_len(rounds)) for (mode in if (round %% 2) names(fits) else rev(names(fits))) {
    gc()
    elapsed <- timed_batch(mode, counts[[mode]])
    seconds[round, mode] <- elapsed / counts[[mode]]
    timings[[length(timings) + 1L]] <- data.frame(case = case, rows = size, round = round,
      mode = mode, fits_in_batch = counts[[mode]], batch_seconds = elapsed, seconds = seconds[round, mode])
  }
  # Identical final-curvature stabilization is checked outside the timed fits.
  stable <- lapply(c("numeric", "analytic"), function(mode) public_fit(mode, stable = TRUE))
  stable_se <- lapply(stable, function(x) suppressWarnings(sqrt(diag(vcov(x)))))
  stable_diagnostics <- lapply(stable, fit_diagnostics)
  stable_delta <- if (all(is.finite(unlist(stable_se)))) max(abs(stable_se[[1]] - stable_se[[2]]) / pmax(stable_se[[1]], 1e-12)) else NA_real_
  for (mode in names(fits)) {
    delta_ll <- abs(as.numeric(logLik(default)) - as.numeric(logLik(fits[[mode]])))
    delta_coef <- max(abs(coef(default) - coef(fits[[mode]])))
    stopifnot(default$code == 0, fits[[mode]]$code == 0, delta_ll < 1e-5, delta_coef < 1e-3)
    row <- data.frame(case = case, rows = size, parameters = length(map$start), mode = mode,
      seconds = median(seconds[, mode]), code = fits[[mode]]$code,
      loglik = as.numeric(logLik(fits[[mode]])), loglik_abs_difference = delta_ll,
      coefficient_max_abs_difference = delta_coef, r_loglik_abs_difference = ref_delta,
      stable_se_max_relative_difference = stable_delta,
      stable_numeric_status = stable_diagnostics[[1]]$status,
      stable_analytic_status = stable_diagnostics[[2]]$status,
      stable_numeric_reason = if (is.null(stable_diagnostics[[1]]$reason)) NA_character_ else stable_diagnostics[[1]]$reason,
      stable_analytic_reason = if (is.null(stable_diagnostics[[2]]$reason)) NA_character_ else stable_diagnostics[[2]]$reason,
      stable_numeric_nonfinite_se = sum(!is.finite(stable_se[[1]])),
      stable_analytic_nonfinite_se = sum(!is.finite(stable_se[[2]])))
    comparisons[[length(comparisons) + 1L]] <- row
    estimates[[length(estimates) + 1L]] <- data.frame(case = case, rows = size, mode = mode,
      parameter = names(coef(default)), estimate = coef(fits[[mode]]),
      default_se = suppressWarnings(sqrt(diag(vcov(fits[[mode]])))), row.names = NULL)
    print(row, row.names = FALSE)
  }
  for (entry in c("timings", "comparisons", "estimates"))
    write.csv(do.call(rbind, get(entry)), file.path(output, paste0(entry, ".csv")), row.names = FALSE)
}
cat("Results:", normalizePath(output, winslash = "/"), "\n")
