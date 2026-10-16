# Rscript --vanilla tools/compare-baseline-projection.R BEFORE AFTER OUT
args <- commandArgs(TRUE); before <- args[1]; after <- args[2]; out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- list.files(before, pattern = "[.]rds$")
stopifnot(length(files) > 0, setequal(files, list.files(after, pattern = "[.]rds$")))
relative <- function(x, y) if (!length(x)) 0 else max(abs(x - y) / pmax(1, abs(y)))
timing <- memory <- checks <- estimates <- fits <- list()
for (file in files) {
  b <- readRDS(file.path(before, file)); a <- readRDS(file.path(after, file))
  stopifnot(identical(b$id, a$id), identical(b$operation, a$operation), identical(b$round, a$round),
    identical(b$warnings, a$warnings))
  meta <- data.frame(id = b$id, operation = b$operation, events = b$events, round = b$round)
  for (version in c("before", "after")) {
    z <- if (version == "before") b else a
    timing[[length(timing) + 1L]] <- cbind(meta, version, seconds = z$seconds, warning_count = length(z$warnings))
    if (startsWith(z$operation, "fit_")) fits[[length(fits) + 1L]] <- cbind(meta, version,
      code = z$value$code, iterations = unname(z$value$iterations), status = z$value$curvature,
      finite_se = sum(is.finite(z$value$se)), warning_count = length(z$warnings),
      warnings = paste(z$warnings, collapse = " | "))
    if (!is.null(z$memory)) memory[[length(memory) + 1L]] <- cbind(meta, version, as.data.frame(as.list(z$memory)))
  }
  x <- a$value; y <- b$value
  if (b$operation == "path") {
    stopifnot(!y$compact, x$compact)
    checks[[length(checks) + 1L]] <- cbind(meta, quantity = "retained_bytes", before = y$bytes, after = x$bytes, error = NA_real_)
  } else if (startsWith(b$operation, "fit_")) {
    stopifnot(identical(x$code, y$code), identical(x$curvature, y$curvature),
      identical(x$iterations, y$iterations), identical(names(x$estimate), names(y$estimate)),
      all(is.finite(x$se)), all(is.finite(y$se)))
    for (quantity in c("estimate", "maximum", "gradient", "hessian", "covariance", "se")) {
      error <- relative(x[[quantity]], y[[quantity]])
      stopifnot(error < if (quantity %in% c("estimate", "maximum", "gradient")) 1e-7 else 1e-5)
      checks[[length(checks) + 1L]] <- cbind(meta, quantity, before = NA_real_, after = NA_real_, error)
    }
    se_relative <- max(abs(x$se - y$se) / abs(y$se))
    stopifnot(is.finite(se_relative), se_relative < 1e-5)
    checks[[length(checks) + 1L]] <- cbind(meta, quantity = "se_relative", before = NA_real_, after = NA_real_, error = se_relative)
    if (b$round == 1L) estimates[[length(estimates) + 1L]] <- cbind(meta,
      parameter = names(x$estimate), estimate_before = y$estimate, estimate_after = x$estimate,
      se_before = y$se, se_after = x$se)
  } else {
    stopifnot(identical(x$parameters, y$parameters), identical(x$loglikelihood, y$loglikelihood))
    error <- relative(x$gradient, y$gradient); stopifnot(error < 1e-10)
    checks[[length(checks) + 1L]] <- cbind(meta, quantity = "gradient", before = NA_real_, after = NA_real_, error)
  }
}
timing <- do.call(rbind, timing); memory <- do.call(rbind, memory); checks <- do.call(rbind, checks)
write.csv(timing, file.path(out, "timings.csv"), row.names = FALSE)
write.csv(memory, file.path(out, "memory.csv"), row.names = FALSE)
write.csv(checks, file.path(out, "checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, fits), file.path(out, "fits.csv"), row.names = FALSE)
write.csv(do.call(rbind, estimates), file.path(out, "estimates.csv"), row.names = FALSE)
medians <- aggregate(seconds ~ id + operation + events + version, subset(timing, round > 0), median)
summary <- reshape(medians, idvar = c("id", "operation", "events"), timevar = "version", direction = "wide")
summary$speedup <- summary$seconds.before / summary$seconds.after
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
print(aggregate(error ~ quantity, checks, max), row.names = FALSE)
print(memory, row.names = FALSE)
