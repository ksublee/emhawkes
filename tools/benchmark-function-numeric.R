# Run in fresh R processes from the repository root.
# Rscript tools/benchmark-function-numeric.R prepare OUTPUT CANDIDATE_LIBRARY
# Rscript tools/benchmark-function-numeric.R worker OUTPUT LIBRARY cran|candidate|previous [case]
# Rscript tools/benchmark-function-numeric.R compare OUTPUT CANDIDATE_LIBRARY CRAN_LIBRARY
args <- commandArgs(TRUE)
mode <- args[1]
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages(library(emhawkes))
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(args[3])))
fields <- c("estimate", "maximum", "gradient", "hessian", "code", "iterations")
if (mode == "prepare") {
  mu <- function(param = c(log_mu = log(.3))) matrix(exp(param))
  alpha <- function(param = c(log_alpha = log(.4))) matrix(exp(param))
  beta <- function(param = c(log_beta = log(1.2))) matrix(exp(param))
  fixed_zero <- function(param = c(log_a = log(.2), log_b = log(.08))) {
    matrix(c(exp(param["log_a"]), 0, exp(param["log_b"]), exp(param["log_a"])), 2)
  }
  decay <- function(param = c(log_decay = 0)) matrix(exp(param), 2, 2)
  shared_mu <- function(param = c(scale = log(.25))) matrix(c(1, 1.25) * exp(param), 2)
  shared_beta <- function(param = c(scale = log(.25))) matrix(1 + exp(param), 2, 2)
  shared_alpha <- function(param = c(a = log(.2), b = log(.06))) {
    matrix(exp(param[c("a", "b", "b", "a")]), 2)
  }
  multi_alpha <- function(param = c(log_slow = log(.1), log_fast = log(.25), log_cross = log(.01))) {
    slow <- exp(param["log_slow"]); fast <- exp(param["log_fast"]); cross <- exp(param["log_cross"])
    matrix(c(slow, cross, fast, cross, cross, slow, cross, fast), 2)
  }
  multi_beta <- function(param = numeric(0)) matrix(rep(c(.8, 4, .8, 4), each = 2), 2)
  cases <- list(
    univariate_log = list(seed = 8741L, initial = matrix(.1),
      truth = new("hspec", mu = .35, alpha = .5, beta = 1.5),
      start = new("hspec", mu = mu, alpha = alpha, beta = beta)),
    bivariate_fixed_zero = list(seed = 8742L, initial = matrix(.1, 2, 2),
      truth = new("hspec", mu = c(.25, .35), alpha = matrix(c(.3, 0, .15, .3), 2), beta = matrix(1.2, 2, 2)),
      start = new("hspec", mu = c(.2, .3), alpha = fixed_zero, beta = decay)),
    bivariate_shared_name = list(seed = 8743L, initial = NULL,
      truth = new("hspec", mu = c(.3, .375), alpha = matrix(c(.25, .08, .08, .25), 2), beta = matrix(1.3, 2, 2)),
      start = new("hspec", mu = shared_mu, alpha = shared_alpha, beta = shared_beta)),
    bivariate_multi_kernel = list(seed = 8744L, initial = matrix(.1, 2, 4),
      truth = new("hspec", mu = c(.25, .35),
        alpha = matrix(c(.12, .015, .35, .015, .015, .12, .015, .35), 2),
        beta = matrix(rep(c(.8, 4, .8, 4), each = 2), 2), type_col_map = list(1:2, 3:4)),
      start = new("hspec", mu = c(.2, .3), alpha = multi_alpha, beta = multi_beta,
        type_col_map = list(1:2, 3:4))))
  fixtures <- list()
  for (name in names(cases)) {
    spec <- cases[[name]]
    set.seed(spec$seed)
    path <- hsim(spec$truth, size = 2001L, lambda_component0 = .1)
    map <- emhawkes:::build_parameter_map(spec$start, TRUE)
    lower <- rep(-8, length(map$start))
    if (map$components$mu$kind == "matrix") lower[unique(map$components$mu$index)] <- 1e-5
    upper <- rep(2, length(map$start))
    p <- length(map$start)
    fixtures[[name]] <- list(object = spec$start, inter_arrival = path$inter_arrival,
      type = path$type, mark = path$mark, lambda_component0 = spec$initial, method = "BFGS",
      constraints = list(ineqA = rbind(diag(p), -diag(p)), ineqB = c(-lower, upper)),
      control = list(iterlim = 500, reltol = 1e-10))
  }
  saveRDS(fixtures, file.path(out, "fixtures.rds"))
  quit()
}
fixtures <- readRDS(file.path(out, "fixtures.rds"))
if (mode == "point_reference") {
  points <- readRDS(file.path(out, "point-models.rds"))
  values <- vapply(points, function(a) as.numeric(suppressMessages(do.call(logLik, a))), numeric(1))
  write.csv(data.frame(point = names(values), cran_loglik = values),
    file.path(out, "cran-point-values.csv"), row.names = FALSE)
  quit()
}
if (mode == "worker") {
  kind <- args[4]
  stopifnot(kind %in% c("cran", "candidate", "previous"))
  selection <- if (length(args) >= 5L) args[5] else names(fixtures)
  stopifnot(all(selection %in% names(fixtures)))
  result_file <- file.path(out, paste0(kind, "-fits.rds"))
  timing_file <- file.path(out, paste0(kind, "-timings.csv"))
  results <- if (file.exists(result_file)) readRDS(result_file) else list()
  timings <- list()
  if (file.exists(timing_file)) {
    existing <- read.csv(timing_file)
    existing <- existing[!existing$case %in% selection, ]
    if (!"fits_in_batch" %in% names(existing)) {
      existing$fits_in_batch <- 1L
      existing$batch_seconds <- existing$elapsed_seconds
    }
    if (nrow(existing)) timings <- split(existing, seq_len(nrow(existing)))
  }
  for (name in selection) {
    a <- fixtures[[name]]
    if (kind == "cran") {
      a$constraint <- a$constraints
      a$constraints <- NULL
    }
    fit_once <- function(engine) {
      previous <- options(emhawkes.engine = engine)
      on.exit(options(previous))
      suppressWarnings(suppressMessages(do.call(hfit, a)))
    }
    engines <- if (kind == "candidate") c(r_numeric = "r", cpp_numeric = "cpp") else
      setNames("auto", kind)
    fits <- list()
    counts <- stats::setNames(rep(1L, length(engines)), names(engines))
    for (engine in names(engines)) {
      warm_time <- unname(system.time(fits[[engine]] <- fit_once(engines[[engine]]))[["elapsed"]])
      if (warm_time < .2) counts[[engine]] <- min(64L, max(1L, ceiling(.5 / max(warm_time, .01))))
    }
    if (kind != "previous") for (round in 1:3) {
      order <- if (round %% 2) names(engines) else rev(names(engines))
      for (engine in order) {
        gc()
        elapsed <- unname(system.time(for (i in seq_len(counts[[engine]])) {
          fit <- fit_once(engines[[engine]])
        })[["elapsed"]])
        fits[[engine]] <- fit
        seconds <- elapsed / counts[[engine]]
        timings[[length(timings) + 1L]] <- data.frame(case = name, round = round, engine = engine,
          elapsed_seconds = seconds, fits_in_batch = counts[[engine]], batch_seconds = elapsed)
        cat(name, round, engine, seconds, "seconds; code", fit$code, "\n")
        flush.console()
      }
    }
    results[[name]] <- lapply(fits, function(fit) fit[fields])
    stopifnot(all(vapply(fits, function(fit) fit$code == 0L, logical(1))))
    saveRDS(results, file.path(out, paste0(kind, "-fits.rds")))
    if (length(timings)) write.csv(do.call(rbind, timings),
      file.path(out, paste0(kind, "-timings.csv")), row.names = FALSE)
  }
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
    file.path(out, paste0("session-", kind, ".txt")))
  quit()
}
if (mode == "compare") {
  cran <- readRDS(file.path(out, "cran-fits.rds"))
  candidate <- readRDS(file.path(out, "candidate-fits.rds"))
  previous <- readRDS(file.path(out, "previous-fits.rds"))
  timings <- rbind(read.csv(file.path(out, "cran-timings.csv")), read.csv(file.path(out, "candidate-timings.csv")))
  rows <- coefficients <- baseline <- engine_fields <- list()
  points <- native_values <- list()
  for (name in names(fixtures)) {
    a <- fixtures[[name]]
    map <- emhawkes:::build_parameter_map(a$object, TRUE)
    observation <- do.call(emhawkes:::normalize_observation,
      a[intersect(names(a), c("object", "inter_arrival", "type", "mark"))])
    native <- emhawkes:::make_hawkes_objective(a$object, observation,
      a$lambda_component0, map, derivatives = FALSE)
    parameters <- list(start = map$start, cran = cran[[name]]$cran$estimate,
      cpp = candidate[[name]]$cpp_numeric$estimate)
    for (label in names(parameters)) {
      point <- parameters[[label]][names(map$start)]
      key <- paste(name, label, sep = ":")
      values <- emhawkes:::restore_parameters(map, point)
      model <- a$object
      for (component in c("mu", "alpha", "beta", "eta")) slot(model, component) <- values[[component]]
      points[[key]] <- c(list(object = model),
        a[intersect(names(a), c("inter_arrival", "type", "mark", "lambda_component0"))])
      native_values[[key]] <- data.frame(point = key, cpp_loglik = as.numeric(native$loglik(point)))
    }
  }
  saveRDS(points, file.path(out, "point-models.rds"))
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(executable, c("--vanilla", "tools/benchmark-function-numeric.R", "point_reference",
    shQuote(out), shQuote(args[4])))
  stopifnot(status == 0L)
  point_values <- merge(read.csv(file.path(out, "cran-point-values.csv")), do.call(rbind, native_values), by = "point")
  point_values$absolute_difference <- abs(point_values$cran_loglik - point_values$cpp_loglik)
  stopifnot(nrow(point_values) == 12L, all(point_values$absolute_difference < 1e-8))
  write.csv(point_values, file.path(out, "pointwise-comparison.csv"), row.names = FALSE)
  se <- function(fit) suppressWarnings(sqrt(diag(-solve(fit$hessian))))
  for (name in names(fixtures)) {
    old <- previous[[name]]$previous
    r <- candidate[[name]]$r_numeric
    cpp <- candidate[[name]]$cpp_numeric
    cr <- cran[[name]]$cran
    same <- vapply(fields, function(field) identical(old[[field]], r[[field]]), logical(1))
    baseline[[name]] <- data.frame(case = name, field = fields, identical = same)
    engine_same <- vapply(fields, function(field) identical(r[[field]], cpp[[field]]), logical(1))
    engine_fields[[name]] <- data.frame(case = name, field = fields, identical = engine_same)
    stopifnot(all(same), all(engine_same))
    for (label in c("r_numeric", "cran")) {
      ref <- if (label == "r_numeric") r else cr
      stopifnot(setequal(names(ref$estimate), names(cpp$estimate)))
      indices <- names(cpp$estimate)
      delta <- max(abs(cpp$estimate - ref$estimate[indices]))
      ll_delta <- abs(cpp$maximum - ref$maximum)
      stopifnot(ll_delta < 1e-5, delta < 1e-3)
      rows[[paste(name, label)]] <- data.frame(case = name, reference = label,
        coefficient_max_abs_difference = delta, loglik_abs_difference = ll_delta,
        reference_code = ref$code, cpp_code = cpp$code,
        reference_iterations = ref$iterations, cpp_iterations = cpp$iterations,
        reference_nan_se = sum(!is.finite(se(ref))), cpp_nan_se = sum(!is.finite(se(cpp))),
        se_max_relative_difference = if (all(is.finite(c(se(ref), se(cpp)))))
          max(abs(se(cpp) - se(ref)[indices]) / pmax(se(ref)[indices], 1e-12)) else NA_real_)
    }
    coefficients[[name]] <- data.frame(case = name, parameter = names(cpp$estimate),
      cran = unname(cr$estimate[names(cpp$estimate)]), r_numeric = unname(r$estimate[names(cpp$estimate)]),
      cpp_numeric = unname(cpp$estimate), cran_se = unname(se(cr)[names(cpp$estimate)]),
      r_numeric_se = unname(se(r)[names(cpp$estimate)]), cpp_numeric_se = unname(se(cpp)))
  }
  write.csv(do.call(rbind, rows), file.path(out, "comparisons.csv"), row.names = FALSE)
  write.csv(do.call(rbind, coefficients), file.path(out, "estimates.csv"), row.names = FALSE)
  write.csv(do.call(rbind, baseline), file.path(out, "baseline-comparison.csv"), row.names = FALSE)
  write.csv(do.call(rbind, engine_fields), file.path(out, "engine-comparison.csv"), row.names = FALSE)
  write.csv(timings, file.path(out, "timings.csv"), row.names = FALSE)
  print(aggregate(elapsed_seconds ~ case + engine, timings, median), row.names = FALSE)
  print(do.call(rbind, rows), row.names = FALSE)
}
