# Rscript --vanilla tools/summarize-simulation-copy.R TASK OUT [NEW_LIB]
args <- commandArgs(TRUE); task <- args[1]; out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
old <- readRDS(file.path(task, "old.checks.rds")); new <- readRDS(file.path(task, "new.checks.rds"))
stopifnot(identical(old, new))
checks <- data.frame(case = names(old), path_warnings_rng_identical = TRUE)
write.csv(checks, file.path(out, "checks.csv"), row.names = FALSE)
files <- list.files(task, pattern = "^(old|new)\\..*\\.[0-9]+\\.rds$", full.names = TRUE)
stopifnot(length(files) == 40L)
values <- lapply(files, readRDS)
timing <- do.call(rbind, lapply(values, `[[`, "timings"))
write.csv(timing, file.path(out, "timings.csv"), row.names = FALSE)
summary <- list()
for (key in unique(paste(timing$case, timing$events))) {
  group <- values[vapply(values, function(x) paste(x$timings$case[1], x$timings$events[1]) == key, FALSE)]
  stopifnot(length(group) == 2L)
  labels <- vapply(group, function(x) x$timings$version[1], "")
  before <- group[[match("old", labels)]]; after <- group[[match("new", labels)]]
  stopifnot(identical(before$data, after$data), identical(before$seed, after$seed),
    identical(before$output_bytes, after$output_bytes))
  old_time <- median(before$timings$elapsed); new_time <- median(after$timings$elapsed)
  summary[[length(summary) + 1L]] <- data.frame(case = before$timings$case[1], events = before$timings$events[1],
    old_seconds = old_time, new_seconds = new_time, speedup = old_time / new_time,
    old_allocated_mib = before$allocated_bytes / 2^20, new_allocated_mib = after$allocated_bytes / 2^20,
    output_mib = before$output_bytes / 2^20, data_rng_identical = TRUE)
}
summary <- do.call(rbind, summary)
summary <- summary[order(summary$events, summary$case), ]
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(trimws(values[[1]]$session, which = "right"), file.path(out, "session-info.txt"))
# Independent inference replay of all twelve representative simulation fixtures.
.libPaths(c(if (length(args) >= 3L) args[3] else file.path(task, "new-lib"), .libPaths())); library(emhawkes)
source("tools/simulation-profile-cases.R")
cases <- unique(timing$case); replay_rows <- list()
for (id in cases) {
  x <- simulation_profile_case(id); set.seed(x$seed)
  path <- hsim(x$object, size = 1001, lambda_component0 = x$initial); seed <- .Random.seed
  for (j in seq_len(x$meta$d)) stopifnot(
    identical(as.numeric(path$N[, j]), cumsum(as.numeric(path$type == j))),
    identical(as.numeric(path$Nc[, j]), cumsum(path$mark * (path$type == j))))
  for (engine in c("r", "cpp")) {
    options(emhawkes.engine = engine)
    replay <- logLik(x$object, path$inter_arrival, path$type, path$mark,
      lambda_component0 = x$initial, infer = TRUE)
    stopifnot(identical(seed, .Random.seed))
    for (field in c("lambda", "lambda_component", "rambda", "rambda_component")) {
      error <- max(abs(replay[[field]] - path[[field]]) / pmax(1, abs(path[[field]])))
      stopifnot(is.finite(error), error < 1e-11)
      replay_rows[[length(replay_rows) + 1L]] <- data.frame(case = id, engine, field, scaled_error = error)
    }
  }
}
write.csv(do.call(rbind, replay_rows), file.path(out, "replay.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
cat(length(old), "edge cases and", nrow(summary), "timed case/size pairs exactly identical.\n")
