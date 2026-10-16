# Rscript --vanilla tools/benchmark-stable-final.R run|worker|repeat OLDLIB NEWLIB FIXTURES PRIOR OUT [ID MODE]
args <- commandArgs(TRUE)
action <- args[1]; oldlib <- args[2]; newlib <- args[3]
fixtures <- args[4]; prior <- args[5]; out <- args[6]
dir.create(file.path(out, "samples"), recursive = TRUE, showWarnings = FALSE)
design <- read.csv(file.path(prior, "design.csv"), stringsAsFactors = FALSE)
# Start with smaller cases, alternating version order across models.
design <- design[order(design$d, design$feedback == "dense", design$observed), ]
fields <- c("estimate", "maximum", "gradient", "code", "iterations")
path <- function(id, mode) file.path(out, "samples", paste0(id, ".", mode, ".rds"))
if (action == "run") {
  write.csv(design, file.path(out, "design.csv"), row.names = FALSE)
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  rows <- list(); coefficients <- list()
  for (i in seq_len(nrow(design))) {
    id <- design$id[i]
    modes <- if (i %% 2L) c("before", "after") else c("after", "before")
    for (mode in modes) {
      if (file.exists(path(id, mode))) next
      cat(format(Sys.time()), "Start", id, mode, "\n"); flush.console()
      status <- system2(rs, vapply(c("--vanilla", "tools/benchmark-stable-final.R", "worker",
        oldlib, newlib, fixtures, prior, out, id, mode), shQuote, ""))
      stopifnot(status == 0L)
    }
    old <- readRDS(path(id, "before")); new <- readRDS(path(id, "after"))
    stopifnot(identical(old$value, new$value), identical(old$covariance, new$covariance),
      identical(old$hessian, new$hessian), identical(old$summary, new$summary),
      identical(old$diagnostics$checks, new$diagnostics$checks),
      new$diagnostics$status == "ok", all(is.finite(new$se)), all(new$se > 0),
      isTRUE(new$diagnostics$computation$hessian$maxLik_skipped),
      is.null(new$diagnostics$original_hessian),
      !is.null(old$diagnostics$original_hessian))
    rows[[i]] <- data.frame(id, dimensions = design$d[i], kernel_columns = design$k[i],
      parameters = length(new$se), events = 2500, feedback = design$feedback[i],
      observed_baseline = design$observed[i], alpha_coupled = design$direct[i],
      before_seconds = old$elapsed, after_seconds = new$elapsed,
      before_cpu_seconds = old$cpu, after_cpu_seconds = new$cpu,
      speedup = old$elapsed / new$elapsed, reduction_percent = 100 * (1 - new$elapsed / old$elapsed),
      max_parameter_difference = max(abs(old$value$estimate - new$value$estimate)),
      loglik_difference = abs(old$value$maximum - new$value$maximum),
      max_se_difference = max(abs(old$se - new$se)), covariance_identical = identical(old$covariance, new$covariance),
      code = new$value$code, nonfinite_se = sum(!is.finite(new$se)))
    coefficients[[i]] <- data.frame(id, parameter = names(new$se), estimate_before = old$value$estimate,
      estimate_after = new$value$estimate, se_before = old$se, se_after = new$se)
    write.csv(do.call(rbind, rows), file.path(out, "comparisons.csv"), row.names = FALSE)
    write.csv(do.call(rbind, coefficients), file.path(out, "coefficients.csv"), row.names = FALSE)
    print(rows[[i]], row.names = FALSE); flush.console()
  }
  quit()
}
if (action == "worker") {
  id <- args[7]; mode <- args[8]
  lib <- if (mode == "before") oldlib else newlib
  .libPaths(c(lib, .libPaths())); library(emhawkes)
  stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
  options(emhawkes.engine = "auto")
  affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
  if (affinity != "scheduler") {
    stopifnot(.Platform$OS.type == "windows", grepl("^[1-9][0-9]*$", affinity))
    actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf(
      "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), affinity))), stdout = TRUE)
    stopifnot(length(actual) == 1L, trimws(actual) == affinity)
  }
  x <- readRDS(file.path(fixtures, "fixtures", paste0(id, ".rds")))
  a <- x$arguments
  for (name in c("inter_arrival", "type", "mark")) a[[name]] <- a[[name]][seq_len(2501L)]
  a$jacobian <- x$jacobian
  warm <- a
  for (name in c("inter_arrival", "type", "mark")) warm[[name]] <- warm[[name]][seq_len(101L)]
  warm$control <- list(iterlim = 2)
  invisible(suppressWarnings(suppressMessages(do.call(hfit, warm))))
  notices <- character()
  invoke <- function(expr) withCallingHandlers(expr,
    warning = function(w) { notices <<- c(notices, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  invisible(gc())
  tm <- system.time(fit <- invoke(do.call(hfit, c(a, list(final_hessian = "stable")))))
  d <- fit_diagnostics(fit)
  previous <- readRDS(file.path(prior, "samples", paste0(id, ".rds")))
  stopifnot(identical(fit[fields], previous$value), identical(vcov(fit), previous$covariance),
    d$status == "ok", d$method == "scaled_score_Richardson",
    d$computation$gradient$used == "numeric", identical(d$computation$engine$used, "cpp"),
    is.null(d$score_fallback_reason), length(notices) == 0L)
  z <- list(id = id, mode = mode, value = fit[fields], covariance = vcov(fit),
    hessian = maxLik::hessian(fit), summary = summary(fit)$estimate,
    se = sqrt(diag(vcov(fit))), diagnostics = d, elapsed = unname(tm[["elapsed"]]),
    cpu = unname(tm[["user.self"]] + tm[["sys.self"]]), notices = notices,
    affinity = affinity, arguments = list(method = a$method, control = a$control,
      events = length(a$inter_arrival) - 1L), session = capture.output(sessionInfo()))
  saveRDS(z, path(id, mode))
  writeLines(trimws(z$session, which = "right"), file.path(out, "session-info.txt"))
  cat(format(Sys.time()), "Finished", id, mode, "seconds", z$elapsed, "\n"); flush.console()
  quit()
}
if (action == "repeat") {
  id <- design$id[1]
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  rows <- list()
  for (trial in 1:3) {
    destination <- file.path(out, "repeats", as.character(trial))
    for (mode in if (trial %% 2L) c("before", "after") else c("after", "before")) {
      cat(format(Sys.time()), "Repeat", trial, id, mode, "\n"); flush.console()
      status <- system2(rs, vapply(c("--vanilla", "tools/benchmark-stable-final.R", "worker",
        oldlib, newlib, fixtures, prior, destination, id, mode), shQuote, ""))
      stopifnot(status == 0L)
    }
    old <- readRDS(file.path(destination, "samples", paste0(id, ".before.rds")))
    new <- readRDS(file.path(destination, "samples", paste0(id, ".after.rds")))
    stopifnot(identical(old$value, new$value), identical(old$covariance, new$covariance),
      identical(old$hessian, new$hessian), identical(old$summary, new$summary))
    rows[[trial]] <- data.frame(id, trial, before_seconds = old$elapsed, after_seconds = new$elapsed,
      before_cpu_seconds = old$cpu, after_cpu_seconds = new$cpu, parameters_identical = TRUE,
      covariance_identical = TRUE)
    write.csv(do.call(rbind, rows), file.path(out, "repeat-timings.csv"), row.names = FALSE)
  }
  print(do.call(rbind, rows), row.names = FALSE)
  quit()
}
stop("Unknown action")
