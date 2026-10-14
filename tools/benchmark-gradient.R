# Install the package first, then run:
# Rscript tools/benchmark-gradient.R [output-directory] [2001,10001] [3] [all|case-name] [all|curvature-only]
# Both fits use the C++ likelihood, BFGS, and the default final Hessian.
library(emhawkes)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else file.path("tools", "gradient-results")
sizes <- if (length(args) > 1L) as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]]) else c(2001L, 10001L)
rounds <- if (length(args) > 2L) as.integer(args[3]) else 3L
selection <- if (length(args) > 3L) args[4] else "all"
curvature_only <- length(args) > 4L && identical(args[5], "curvature-only")
stopifnot(length(sizes) > 0L, !anyNA(sizes), all(sizes >= 2L),
          length(rounds) == 1L, !is.na(rounds), rounds >= 3L)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "cpp")

unary_truth <- new("hspec", mu = .35, alpha = .5, beta = 1.5)
unary_start <- new("hspec", mu = .3, alpha = .4, beta = 1.2)
bivariate_truth <- new("hspec", mu = c(.25, .35),
  alpha = matrix(c(.25, .08, .08, .25), 2), beta = matrix(1.2, 2, 2))
bivariate_start <- new("hspec", mu = c(.2, .3),
  alpha = matrix(c(.2, .06, .06, .2), 2), beta = matrix(1, 2, 2))
multi_truth <- new("hspec", mu = c(.25, .35), alpha = matrix(.08, 2, 4),
  beta = matrix(rep(c(.6, 1.5, 3, 6), each = 2), 2),
  type_col_map = list(c(1, 3), c(2, 4)))
multi_start <- multi_truth
multi_start@mu <- matrix(c(.2, .3), 2)
multi_start@alpha[] <- .06
free_multi_truth <- new("hspec", mu = .3, alpha = matrix(c(.25, 1.1), 1),
  beta = matrix(c(.8, 5), 1), type_col_map = list(1:2))
free_multi_start <- new("hspec", mu = .25, alpha = matrix(c(.2, .9), 1),
  beta = matrix(c(.7, 4), 1), type_col_map = list(1:2))
cases <- list(
  univariate_fixed_initial = list(truth = unary_truth, start = unary_start, seed = 8741L, initial = .1),
  univariate_default_initial = list(truth = unary_truth, start = unary_start, seed = 8741L, initial = NULL),
  bivariate_shared = list(truth = bivariate_truth, start = bivariate_start, seed = 8742L, initial = .1),
  bivariate_default_initial = list(truth = bivariate_truth, start = bivariate_start, seed = 8742L, initial = NULL),
  bivariate_unshared_fixed_beta = list(truth = bivariate_truth, start = bivariate_start,
    seed = 8742L, initial = .1, reduced = FALSE, fix_beta = TRUE),
  multi_kernel_fixed_beta = list(truth = multi_truth, start = multi_start,
    seed = 8743L, initial = .1, fix_beta = TRUE),
  multi_kernel_free_beta = list(truth = free_multi_truth, start = free_multi_start,
    seed = 8744L, initial = .1, order_beta = TRUE))
if (selection != "all") {
  if (!selection %in% names(cases)) stop("Unknown case: ", selection)
  cases <- cases[selection]
}

timings <- comparisons <- estimates <- conditions <- list()
store_results <- function() {
  for (entry in c("timings", "comparisons", "estimates", "conditions")) {
    rows <- get(entry)
    if (length(rows)) write.csv(do.call(rbind, rows), file.path(output, paste0(entry, ".csv")), row.names = FALSE)
  }
}
capture.output(sessionInfo(), Sys.info()[c("sysname", "release", "version", "machine")],
  file = file.path(output, "session-info.txt"))
if (!curvature_only) for (size in sizes) for (case in names(cases)) {
  spec <- cases[[case]]
  reduced <- !identical(spec$reduced, FALSE)
  map <- emhawkes:::build_parameter_map(spec$start, reduced)
  count <- length(map$start)
  fixed <- rep(FALSE, count)
  if (isTRUE(spec$fix_beta)) fixed[map$components$beta$index] <- TRUE
  constraints <- list(ineqA = diag(count), ineqB = rep(-1e-5, count))
  if (isTRUE(spec$order_beta)) {
    order <- numeric(count)
    order[map$components$beta$index] <- c(-1, 1)
    constraints$ineqA <- rbind(constraints$ineqA, order)
    constraints$ineqB <- c(constraints$ineqB, -.2)
  }
  set.seed(spec$seed)
  # The observation design and simulation initial state are identical across modes.
  path <- hsim(spec$truth, size = size, lambda_component0 = .1)
  fit_once <- function(mode, final_hessian = "maxLik") {
    hfit(spec$start, path$inter_arrival, path$type,
      lambda_component0 = spec$initial, reduced = reduced, fixed = fixed,
      method = "BFGS", gradient = mode, final_hessian = final_hessian,
      constraints = constraints, control = list(iterlim = 500, reltol = 1e-10))
  }
  fit_quietly <- function(mode, record = FALSE, final_hessian = "maxLik") {
    withCallingHandlers(fit_once(mode, final_hessian),
      message = function(e) invokeRestart("muffleMessage"),
      warning = function(e) {
        if (record) conditions[[length(conditions) + 1L]] <<- data.frame(
          case = case, rows = size, mode = mode, message = conditionMessage(e))
        invokeRestart("muffleWarning")
      })
  }
  # Untimed warm-ups also retain the reproducibility and uncertainty comparison.
  fits <- lapply(c("numeric", "analytic"), fit_quietly, record = TRUE)
  names(fits) <- c("numeric", "analytic")
  timed_batch <- function(mode, count) {
    unname(system.time(for (i in seq_len(count)) fit_quietly(mode))[["elapsed"]])
  }
  # Batch short fits so the system clock's resolution cannot manufacture a speedup.
  batch_counts <- vapply(names(fits), function(mode) {
    count <- 1L
    repeat {
      if (timed_batch(mode, count) >= .2 || count >= 128L) return(count)
      count <- count * 2L
    }
  }, integer(1))
  seconds <- matrix(NA_real_, rounds, 2L, dimnames = list(NULL, names(fits)))
  for (round in seq_len(rounds)) {
    modes <- if (round %% 2L) names(fits) else rev(names(fits))
    for (mode in modes) {
      gc()
      elapsed <- timed_batch(mode, batch_counts[[mode]])
      seconds[round, mode] <- elapsed / batch_counts[[mode]]
      timings[[length(timings) + 1L]] <- data.frame(case = case, rows = size,
        round = round, mode = mode, fits_in_batch = batch_counts[[mode]],
        batch_seconds = elapsed, seconds = seconds[round, mode])
    }
  }
  numeric_est <- coef(fits$numeric)
  analytic_est <- coef(fits$analytic)
  standard_errors <- lapply(fits, function(fit) suppressWarnings(sqrt(diag(vcov(fit)))))
  comparable <- !fixed & is.finite(standard_errors$numeric) & is.finite(standard_errors$analytic)
  se_relative_difference <- if (any(comparable)) max(
    abs(standard_errors$analytic[comparable] - standard_errors$numeric[comparable]) /
      pmax(abs(standard_errors$numeric[comparable]), 1e-12)) else NA_real_
  # maxLik differentiates the supplied gradient for its final Hessian. Check
  # common likelihood-based stabilization separately and outside timed fits.
  stable_fits <- lapply(names(fits), fit_quietly, final_hessian = "stable")
  names(stable_fits) <- names(fits)
  stable_se <- lapply(stable_fits, function(fit) suppressWarnings(sqrt(diag(vcov(fit)))))
  stable_comparable <- !fixed & is.finite(stable_se$numeric) & is.finite(stable_se$analytic)
  stable_se_difference <- if (any(stable_comparable)) max(
    abs(stable_se$analytic[stable_comparable] - stable_se$numeric[stable_comparable]) /
      pmax(abs(stable_se$numeric[stable_comparable]), 1e-12)) else NA_real_
  row <- data.frame(case = case, rows = size, parameters = count, free_parameters = sum(!fixed),
    numeric_seconds = median(seconds[, "numeric"]), analytic_seconds = median(seconds[, "analytic"]),
    speedup = median(seconds[, "numeric"]) / median(seconds[, "analytic"]),
    numeric_code = fits$numeric$code, analytic_code = fits$analytic$code,
    numeric_iterations = fits$numeric$iterations, analytic_iterations = fits$analytic$iterations,
    numeric_loglik = as.numeric(logLik(fits$numeric)), analytic_loglik = as.numeric(logLik(fits$analytic)),
    loglik_abs_difference = abs(as.numeric(logLik(fits$numeric)) - as.numeric(logLik(fits$analytic))),
    coefficient_max_abs_difference = max(abs(numeric_est - analytic_est)),
    coefficient_max_scaled_difference = max(abs(numeric_est - analytic_est) / pmax(1, abs(numeric_est))),
    se_max_relative_difference = se_relative_difference,
    numeric_nonfinite_se = sum(!fixed & !is.finite(standard_errors$numeric)),
    analytic_nonfinite_se = sum(!fixed & !is.finite(standard_errors$analytic)),
    stable_se_max_relative_difference = stable_se_difference,
    stable_numeric_nonfinite_se = sum(!fixed & !is.finite(stable_se$numeric)),
    stable_analytic_nonfinite_se = sum(!fixed & !is.finite(stable_se$analytic)))
  comparisons[[length(comparisons) + 1L]] <- row
  estimates[[length(estimates) + 1L]] <- data.frame(case = case, rows = size,
    parameter = names(numeric_est), fixed = fixed, numeric = as.numeric(numeric_est),
    analytic = as.numeric(analytic_est), numeric_se = standard_errors$numeric,
    analytic_se = standard_errors$analytic, stable_numeric_se = stable_se$numeric,
    stable_analytic_se = stable_se$analytic, row.names = NULL)
  print(row, row.names = FALSE)
  store_results()
}
# A separate R central-difference implementation differentiates the analytic
# score at three steps. Its information matrix is compared with both default
# and stabilized SEs; this work is outside all timing measurements.
reference_results <- list()
saved_estimates <- read.csv(file.path(output, "estimates.csv"))
relative_error <- function(candidate, reference) {
  keep <- is.finite(candidate) & is.finite(reference)
  if (!any(keep)) return(NA_real_)
  max(abs(candidate[keep] - reference[keep]) / pmax(abs(reference[keep]), 1e-12))
}
for (size in sizes) for (case in names(cases)) {
  spec <- cases[[case]]
  map <- emhawkes:::build_parameter_map(spec$start, !identical(spec$reduced, FALSE))
  saved <- saved_estimates[saved_estimates$case == case & saved_estimates$rows == size, ]
  if (!nrow(saved)) next
  saved <- saved[match(names(map$start), saved$parameter), ]
  theta <- setNames(saved$analytic, saved$parameter)
  free <- which(!saved$fixed)
  set.seed(spec$seed)
  path <- hsim(spec$truth, size = size, lambda_component0 = .1)
  observation <- emhawkes:::normalize_observation(spec$start, path$inter_arrival, path$type)
  objective <- emhawkes:::make_hawkes_objective(spec$start, observation, spec$initial, map)
  previous_se <- NULL
  for (step in c(1e-3, 3e-4, 1e-4)) {
    hessian <- vapply(free, function(j) {
      delta <- pmax(abs(theta[j]), .1) * step
      up <- down <- theta
      up[j] <- up[j] + delta
      down[j] <- down[j] - delta
      (objective$gradient(up)[free] - objective$gradient(down)[free]) / (2 * delta)
    }, numeric(length(free)))
    information <- -(hessian + t(hessian)) / 2
    eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
    reference_se <- if (min(eigenvalues) > 0) sqrt(diag(solve(information))) else rep(NA_real_, length(free))
    reference_results[[length(reference_results) + 1L]] <- data.frame(case = case, rows = size,
      step = step, min_information_eigenvalue = min(eigenvalues),
      reference_se_relative_step_change = if (is.null(previous_se)) NA_real_ else relative_error(reference_se, previous_se),
      numeric_default_se_relative_error = relative_error(saved$numeric_se[free], reference_se),
      analytic_default_se_relative_error = relative_error(saved$analytic_se[free], reference_se),
      numeric_stable_se_relative_error = relative_error(saved$stable_numeric_se[free], reference_se),
      analytic_stable_se_relative_error = relative_error(saved$stable_analytic_se[free], reference_se))
    previous_se <- reference_se
  }
}
write.csv(do.call(rbind, reference_results), file.path(output, "hessian-reference.csv"), row.names = FALSE)
cat("Results:", normalizePath(output, winslash = "/"), "\n")
