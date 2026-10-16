# Rscript --vanilla tools/benchmark-score-curvature.R run|compare LIB FIXTURES PRIOR OUT
# FIXTURES: completed gradient-scaling output; PRIOR: stable-scaling output.
args <- commandArgs(TRUE)
action <- args[1]; lib <- args[2]; fixtures <- args[3]; prior <- args[4]; out <- args[5]
.libPaths(c(lib, .libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
options(emhawkes.engine = "auto")
dir.create(file.path(out, "samples"), recursive = TRUE, showWarnings = FALSE)
design <- read.csv(file.path(prior, "design.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(design) == 5L, all(design$events_used == 2500L))
path <- function(id) file.path(out, "samples", paste0(id, ".rds"))
fields <- c("estimate", "maximum", "gradient", "code", "iterations")
if (action == "run") {
  write.csv(design, file.path(out, "design.csv"), row.names = FALSE)
  writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "session-info.txt"))
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  for (id in design$id) {
    if (file.exists(path(id))) next
    cat("Start", id, "\n"); flush.console()
    status <- system2(rs, vapply(c("--vanilla", "tools/benchmark-score-curvature.R", "worker",
      lib, fixtures, prior, out, id), shQuote, ""))
    stopifnot(status == 0L)
  }
  quit()
}
if (action == "worker") {
  id <- args[6]
  x <- readRDS(file.path(fixtures, "fixtures", paste0(id, ".rds")))
  a <- x$arguments
  for (name in c("inter_arrival", "type", "mark")) a[[name]] <- a[[name]][seq_len(2501L)]
  a$jacobian <- x$jacobian
  affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
  if (affinity != "scheduler") {
    stopifnot(.Platform$OS.type == "windows", grepl("^[1-9][0-9]*$", affinity))
    actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf(
      "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), affinity))), stdout = TRUE)
    stopifnot(length(actual) == 1L, trimws(actual) == affinity)
  }
  warm <- a
  for (name in c("inter_arrival", "type", "mark")) warm[[name]] <- warm[[name]][seq_len(101L)]
  warm$control <- list(iterlim = 2)
  invisible(suppressWarnings(suppressMessages(do.call(hfit, warm))))
  notices <- character()
  invoke <- function(expr) withCallingHandlers(expr,
    warning = function(w) { notices <<- c(notices, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  invisible(gc())
  full_time <- system.time(fit <- invoke(do.call(hfit, c(a, list(final_hessian = "stable")))))
  d <- fit_diagnostics(fit)
  old <- readRDS(file.path(prior, "samples", paste0(id, ".stable.rds")))
  stopifnot(identical(fit[fields], old$value[fields]),
    identical(d$original_hessian, old$diagnostics$original_hessian),
    d$status == "ok", d$method == "scaled_score_Richardson",
    d$computation$gradient$used == "numeric", d$computation$engine$used == "cpp",
    is.null(d$score_fallback_reason))
  saveRDS(fit, file.path(out, "samples", paste0(id, ".fit.rds")))
  cat("Public fit", id, "seconds", full_time[["elapsed"]], "\n"); flush.console()

  # Matched final-curvature timings at exactly the same fitted vector. No refit,
  # tracing, or competing worker; include analytic-adapter construction.
  base <- fit; base$hessian <- d$original_hessian; base$varcovar <- NULL
  obs <- do.call(emhawkes:::normalize_observation, a[intersect(names(a),
    c("object", "inter_arrival", "type", "mark", "N", "Nc", "N0", "Nc0"))])
  map <- emhawkes:::build_parameter_map(a$object, TRUE)
  factory <- function() emhawkes:::make_hawkes_objective(a$object, obs,
    a$lambda_component0, map, jacobian = x$jacobian)
  finalize <- function(score) emhawkes:::finalize_hawkes_fit(base, a$object,
    a$constraints, NULL, list(), "BFGS", "stable", score_factory = score)
  # Alternate order across models to avoid a fixed ordering bias.
  modes <- if (match(id, design$id) %% 2L) c("likelihood", "score") else c("score", "likelihood")
  timing <- list(); curvature <- list()
  for (mode in modes) {
    invisible(gc())
    tm <- system.time(z <- invoke(finalize(if (mode == "score") factory else NULL)))
    timing[[mode]] <- unname(tm[c("elapsed", "user.self", "sys.self")])
    curvature[[mode]] <- list(diagnostics = z$diagnostics, covariance = vcov(z$fit),
      se = sqrt(diag(vcov(z$fit))))
    stopifnot(identical(z$fit[fields], base[fields]), z$diagnostics$status == "ok")
    cat("Curvature", id, mode, "seconds", tm[["elapsed"]], "\n"); flush.console()
  }
  stopifnot(identical(vcov(fit), curvature$score$covariance),
    identical(summary(fit)$estimate[, 2], sqrt(diag(vcov(fit)))))
  saveRDS(list(id = id, full_seconds = unname(full_time[["elapsed"]]),
    full_cpu_seconds = unname(full_time[["user.self"]] + full_time[["sys.self"]]),
    value = fit[fields], diagnostics = d, covariance = vcov(fit),
    se = sqrt(diag(vcov(fit))), timing = timing, curvature = curvature,
    notices = unique(notices), affinity = affinity), path(id))
  quit()
}
if (action == "compare") {
  rows <- coefficients <- steps <- notices <- list()
  reference <- read.csv(file.path(prior, "reference-standard-errors.csv"), stringsAsFactors = FALSE)
  for (id in design$id) {
    z <- readRDS(path(id)); old <- readRDS(file.path(prior, "samples", paste0(id, ".stable.rds")))
    same <- readRDS(file.path(prior, "samples", paste0(id, ".maxLik.rds")))
    stopifnot(identical(z$value, old$value[fields]), identical(z$value, same$value[fields]),
      z$diagnostics$status == "ok", all(is.finite(z$se)), all(z$se > 0))
    r <- reference[reference$id == id, ]; r <- r[match(names(z$se), r$parameter), ]
    error <- max(abs(z$se - old$se) / old$se)
    paired_error <- max(abs(z$curvature$score$se - z$curvature$likelihood$se) / z$curvature$likelihood$se)
    ref_error <- max(abs(z$se - r$se_score_reference) / r$se_score_reference)
    stopifnot(error < .001, paired_error < .001, ref_error < .001,
      identical(z$curvature$likelihood$covariance, old$covariance))
    row <- design[design$id == id, ]
    rows[[length(rows) + 1L]] <- data.frame(id, dimensions = row$d, kernel_columns = row$k,
      parameters = length(z$se), events_used = 2500, feedback = row$feedback,
      observed_baseline = row$observed, alpha_coupled = row$direct,
      full_seconds_previous_stable = old$seconds, full_seconds_score_stable = z$full_seconds,
      full_cpu_seconds_score_stable = z$full_cpu_seconds,
      full_time_ratio = old$seconds / z$full_seconds,
      curvature_seconds_likelihood = z$timing$likelihood[1], curvature_seconds_score = z$timing$score[1],
      curvature_cpu_likelihood = sum(z$timing$likelihood[2:3]),
      curvature_cpu_score = sum(z$timing$score[2:3]),
      curvature_speedup = z$timing$likelihood[1] / z$timing$score[1],
      max_relative_se_difference = error, max_relative_se_difference_paired = paired_error,
      max_relative_se_reference_difference = ref_error,
      relative_covariance_difference = norm(z$covariance - old$covariance, "F") / norm(old$covariance, "F"),
      max_parameter_difference = max(abs(z$value$estimate - old$value$estimate)),
      loglik_difference = abs(z$value$maximum - old$value$maximum),
      steps_score = nrow(z$diagnostics$checks), steps_likelihood = nrow(old$diagnostics$checks),
      code = z$value$code, nonfinite_se = sum(!is.finite(z$se)))
    coefficients[[length(coefficients) + 1L]] <- data.frame(id, parameter = names(z$se),
      estimate = as.numeric(z$value$estimate), se_likelihood = as.numeric(old$se), se_score = as.numeric(z$se),
      se_reference = r$se_score_reference)
    steps[[length(steps) + 1L]] <- data.frame(id, source = "analytic_score", z$diagnostics$checks)
    for (note in z$notices) notices[[length(notices) + 1L]] <- data.frame(id, text = note)
  }
  for (name in c("rows", "coefficients", "steps", "notices")) {
    value <- get(name); value <- if (length(value)) do.call(rbind, value) else data.frame(id = character(), text = character())
    write.csv(value, file.path(out, paste0(if (name == "rows") "comparisons" else name, ".csv")), row.names = FALSE)
  }
  print(do.call(rbind, rows), row.names = FALSE)
  quit()
}
stop("Unknown action: ", action)
