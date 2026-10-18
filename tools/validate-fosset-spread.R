# Rscript --vanilla tools/validate-fosset-spread.R LIB OUT [equations]
# Run against the committed package plus inst/examples/fosset-spread.R.
args <- commandArgs(TRUE); task_out <- args[2]
.libPaths(c(args[1], .libPaths())); library(emhawkes)
dir.create(task_out, recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(args[1])))
source(system.file("examples", "fosset-spread.R", package = "emhawkes"))
checks <- new.env(parent = globalenv()); sys.source("tests/fosset-spread.R", envir = checks)
cat("Equation, boundary, replay, conditional MLE and curvature checks passed.\n")

# Direct interval integrals also verify the per-type residuals, not only total LL.
equation_rows <- list()
for (s0 in c(1, 3)) for (a in c(.2, .4, .7)) for (b in c(.5, 2, 3)) {
  dt <- c(0, .2, .1, .4, .5, .15, .3, .2)
  type <- c(0, 1, 2, 1, 1, 2, 2, 1); time <- cumsum(dt)
  spread <- s0 + cumsum(type == 1) - cumsum(type == 2)
  opening <- which(type == 1)
  theta <- c(log_nu_plus = 0, log_nu_minus = log(3), logit_a = qlogis(a), log_b = log(b))
  expected <- checks$fosset_reference(theta, dt, type, s0)
  integral <- matrix(0, length(dt), 2)
  for (n in 2:length(dt)) {
    past <- opening[opening < n]
    integral[n, 1] <- dt[n] + sum(a * (exp(-b * (time[n - 1] - time[past])) -
      exp(-b * (time[n] - time[past]))))
    integral[n, 2] <- 3 * (spread[n - 1] >= 2) * dt[n]
  }
  for (engine in c("r", "cpp")) {
    options(emhawkes.engine = engine)
    actual <- logLik(fosset_spec(a = a, b = b, spread0 = s0), dt, type,
      lambda_component0 = 0, infer = TRUE)
    residual_error <- 0
    for (k in 1:2) {
      index <- which(type == k)
      residual <- vapply(seq_len(length(index) - 1L), function(i)
        sum(integral[(index[i] + 1L):index[i + 1L], k]), 0)
      residual_error <- max(residual_error, abs(actual[[paste0("res_process", k)]] - residual))
    }
    equation_rows[[length(equation_rows) + 1L]] <- data.frame(spread0 = s0, a, b, engine,
      lambda_error = max(abs(actual$lambda - expected$lambda)),
      rambda_error = max(abs(actual$rambda - expected$rambda)),
      loglik_error = abs(actual$loglikelihood - expected$loglik), residual_error)
  }
}
equation_rows <- do.call(rbind, equation_rows)
stopifnot(max(as.matrix(equation_rows[, 5:8])) < 1e-10)
write.csv(equation_rows, file.path(task_out, "equations.csv"), row.names = FALSE)
if (length(args) >= 3L && args[3] == "equations") quit(status = 0)
options(emhawkes.engine = "auto")

set.seed(20261018)
path <- hsim(fosset_spec(), size = 4001, lambda_component0 = 0)
spread <- 1 + path$N[, 1] - path$N[, 2]
opening <- path$arrival[path$type == 1]
horizon <- tail(path$arrival, 1)
exposure <- sum(path$inter_arrival[-1] * (head(spread, -1) >= 2))
closing_count <- sum(path$type == 2)
true <- c(nu_plus = 1, nu_minus = 3, a = .4, b = 2)

fosset_curvature_check <- function(fit, loglik) {
  information <- stats::optimHess(coef(fit), function(theta) -loglik(theta),
    control = list(ndeps = rep(1e-3, 4)))
  package_se <- sqrt(diag(vcov(fit)))
  independent_se <- sqrt(diag(solve(information)))
  relative_error <- abs(package_se - independent_se) / package_se
  stopifnot(all(is.finite(independent_se)), max(relative_error) < 2e-4)
  data.frame(parameter = names(coef(fit)), package_se, independent_se, relative_error,
    row.names = NULL)
}

# Independent likelihood: only opening times enter the Hawkes recursion.
# The integral uses a closed-form sum over event times, not interval recursion.
reference_loglik <- function(theta) {
  p <- fosset_parameters(theta)
  memory <- numeric(length(opening))
  if (length(opening) > 1L) for (i in 2:length(opening))
    memory[i] <- exp(-p["b"] * (opening[i] - opening[i - 1])) * (1 + memory[i - 1])
  sum(log(p["nu_plus"] + p["a"] * p["b"] * memory)) - p["nu_plus"] * horizon -
    p["a"] * sum(-expm1(-p["b"] * (horizon - opening))) +
    closing_count * log(p["nu_minus"]) - p["nu_minus"] * exposure
}
starts <- list(c(.7, 2.5, .25, 1.2), c(1.4, 4, .55, 3), c(.5, 2, .2, .8))
rows <- list(); fit_records <- list(); fits <- list()
for (i in seq_along(starts)) {
  z <- starts[[i]]; h <- fosset_spec(z[1], z[2], z[3], z[4])
  for (engine in if (i == 1L) c("auto", "r") else "auto") {
    options(emhawkes.engine = engine)
    elapsed <- system.time(fit <- hfit(h, path$inter_arrival, path$type,
      lambda_component0 = 0, optimizer = "stats", method = "BFGS",
      control = list(maxit = 500, reltol = 1e-10)))["elapsed"]
    diagnostic <- fit_diagnostics(fit); theta <- coef(fit); p <- fosset_parameters(theta)
    stopifnot(fit$code == 0, diagnostic$status == "ok", all(is.finite(vcov(fit))),
      identical(names(theta), c("log_nu_plus", "log_nu_minus", "logit_a", "log_b")))
    se_theta <- sqrt(diag(vcov(fit)))
    se <- se_theta * c(p["nu_plus"], p["nu_minus"], p["a"] * (1 - p["a"]), p["b"])
    names(se) <- names(p)
    stopifnot(abs(p["nu_minus"] - closing_count / exposure) < 2e-4,
      abs(vcov(fit)[2, 2] - 1 / closing_count) < 1e-6,
      abs(as.numeric(logLik(fit)) - reference_loglik(theta)) < 1e-7)
    replay <- do.call(logLik, c(list(object = fitted_hspec(fit), inter_arrival = path$inter_arrival,
      type = path$type, infer = TRUE), fit_initial_state(fit)))
    stopifnot(abs(replay$loglikelihood - as.numeric(logLik(fit))) < 1e-7)
    key <- paste(engine, i, sep = "."); fits[[key]] <- fit
    rows[[length(rows) + 1L]] <- data.frame(engine, start = i, parameter = names(p),
      true = true, estimate = p, se, row.names = NULL)
    fit_records[[length(fit_records) + 1L]] <- data.frame(engine, start = i,
      actual_engine = paste(diagnostic$computation$engine$used, collapse = ","),
      elapsed = unname(elapsed), loglik = as.numeric(logLik(fit)),
      code = fit$code, curvature = diagnostic$status,
      narrowing_mle = unname(closing_count / exposure))
    cat(key, "fit completed:", paste(round(p, 5), collapse = ", "), "seconds", elapsed, "\n")
  }
}
fit_records <- do.call(rbind, fit_records)
stopifnot(diff(range(fit_records$loglik)) < 1e-5)
base_fit <- fits[[1]]
independent <- optim(coef(base_fit) + c(.03, -.03, .03, -.03),
  function(theta) -reference_loglik(theta), method = "BFGS",
  control = list(maxit = 500, reltol = 1e-12))
stopifnot(independent$convergence == 0, abs(-independent$value - as.numeric(logLik(base_fit))) < 1e-5,
  max(abs(fosset_parameters(independent$par) - fosset_parameters(coef(base_fit)))) < 1e-3)
write.csv(fosset_curvature_check(base_fit, reference_loglik),
  file.path(task_out, "curvature.csv"), row.names = FALSE)
write.csv(do.call(rbind, rows), file.path(task_out, "estimates.csv"), row.names = FALSE)
write.csv(fit_records, file.path(task_out, "fits.csv"), row.names = FALSE)
write.csv(data.frame(parameter = names(true), package = fosset_parameters(coef(base_fit)),
  independent = fosset_parameters(independent$par)), file.path(task_out, "independent-fit.csv"), row.names = FALSE)

# Exact global-hazard increments and conditional event probabilities.
# These are distribution checks, not exact comparisons between different samplers.
integrated <- path$inter_arrival[-1] * (1 + 3 * (head(spread, -1) >= 2)) +
  (head(path$rambda[, 1], -1) - 1) * (-expm1(-2 * path$inter_arrival[-1])) / 2
probability <- path$lambda[-1, 1] / rowSums(path$lambda[-1, , drop = FALSE])
martingale_z <- sum((path$type[-1] == 1) - probability) / sqrt(sum(probability * (1 - probability)))
stopifnot(abs(mean(integrated) - 1) < 5 / sqrt(length(integrated)), abs(martingale_z) < 5)
write.csv(data.frame(events = 4000, min_spread = min(spread),
  boundary_visits = sum(spread == 1), hazard_increment_mean = mean(integrated),
  event_type_martingale_z = martingale_z), file.path(task_out, "simulation.csv"), row.names = FALSE)

# Time-weighted occupation after burn-in, rather than the event-row proportion.
options(emhawkes.engine = "auto")
set.seed(8917); long <- hsim(fosset_spec(), size = 40001, lambda_component0 = 0)
long_spread <- 1 + long$N[, 1] - long$N[, 2]
keep <- 2002:length(long$type)
occupation <- sum(long$inter_arrival[keep] * (long_spread[keep - 1] >= 2)) /
  sum(long$inter_arrival[keep])
theoretical <- 1 / (3 * (1 - .4))
stopifnot(abs(occupation - theoretical) < .035, min(long_spread) == 1)
write.csv(data.frame(events = 40000, burn_in_events = 2000,
  open_probability = occupation, theoretical), file.path(task_out, "occupation.csv"), row.names = FALSE)

# Small fixed-seed recovery experiment: report all draws, including variation.
mc <- list()
for (seed in 4101:4106) {
  set.seed(seed); sample <- hsim(fosset_spec(), size = 3001, lambda_component0 = 0)
  fit <- hfit(fosset_spec(.7, 2.5, .25, 1.2), sample$inter_arrival, sample$type,
    lambda_component0 = 0, optimizer = "stats", method = "BFGS",
    control = list(maxit = 500, reltol = 1e-10))
  p <- fosset_parameters(coef(fit)); diagnostic <- fit_diagnostics(fit)
  stopifnot(fit$code == 0, diagnostic$status == "ok", all(is.finite(vcov(fit))))
  mc[[length(mc) + 1L]] <- data.frame(seed, parameter = names(p), true = true,
    estimate = p, row.names = NULL)
  cat("Recovery seed", seed, "completed.\n")
}
write.csv(do.call(rbind, mc), file.path(task_out, "recovery.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(task_out, "session-info.txt"))
saveRDS(list(path = path, fits = fits, independent = independent), file.path(task_out, "raw.rds"))
cat("Fosset spread validation passed.\n")
