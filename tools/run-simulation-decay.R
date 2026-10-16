# Rscript --vanilla tools/run-simulation-decay.R BASE_LIB NEW_LIB RAW OUT
# Baseline: ec7768b. Run without other heavy jobs after compiling both versions.
a <- commandArgs(TRUE); stopifnot(length(a) == 4L)
libs <- setNames(normalizePath(a[1:2], winslash = "/", mustWork = TRUE), c("base", "native"))
raw <- a[3]; out <- a[4]
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
worker <- function(script, args, log) {
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(args)), stdout = log, stderr = log)
  if (status != 0) stop("Worker failed: ", log)
}
for (label in names(libs)) worker("tools/check-simulation-decay.R", c(libs[[label]], raw, label),
  file.path(raw, paste0(label, ".checks.log")))
checks <- lapply(names(libs), function(label) readRDS(file.path(raw, paste0(label, ".decay-checks.rds"))))
stopifnot(identical(checks[[1]], checks[[2]]))
write.csv(data.frame(combinations = length(checks[[1]]), output_identical = TRUE, rng_identical = TRUE),
  file.path(out, "validation.csv"), row.names = FALSE)
rm(checks); gc()
design <- rbind(
  data.frame(case = c("univariate", "multi_equal", "large_equal", "wide32", "wide128",
    "state_current", "mark_only", "marked", "heavy_impact", "impact32",
    "multi_unequal", "state_dots"), events = 10000L),
  data.frame(case = c("multi_equal", "large_equal", "state_current"), events = 100000L),
  data.frame(case = c("multi_equal", "state_current"), events = 100L),
  data.frame(case = "multi_equal", events = 0L))
write.csv(design, file.path(out, "design.csv"), row.names = FALSE)
metrics <- list()
pair <- function(id, events, round, mode, order) {
  values <- list()
  for (label in order) {
    key <- paste(label, id, events, round, mode, sep = ".")
    worker("tools/benchmark-native-simulation.R", c(libs[[label]], raw, label, id, events, round,
      mode, "tools/simulation-decay-cases.R"), file.path(raw, paste0(key, ".log")))
    values[[label]] <- readRDS(file.path(raw, paste0(key, ".rds")))
    metrics[[length(metrics) + 1L]] <<- values[[label]]$metrics
  }
  stopifnot(identical(values$base$data, values$native$data), identical(values$base$seed, values$native$seed))
  write.csv(do.call(rbind, metrics), file.path(out, "measurements.csv"), row.names = FALSE)
  cat(id, events, round, mode, "exact output/RNG match\n"); flush.console()
}
for (i in seq_len(nrow(design))) for (round in 1:3) {
  order <- if ((i + round) %% 2L) names(libs) else rev(names(libs))
  pair(design$case[i], design$events[i], round, "timing", order); gc()
}
for (mode in c("memory", "rss")) for (id in c("multi_equal", "wide128", "state_current")) {
  pair(id, 10000L, 0L, mode, names(libs)); gc()
}
metrics <- do.call(rbind, metrics)
summary <- aggregate(cbind(elapsed, cpu) ~ version + case + events,
  subset(metrics, mode == "timing"), median)
wide <- merge(subset(summary, version == "base"), subset(summary, version == "native"),
  by = c("case", "events"), suffixes = c("_base", "_native"))
wide$speedup <- wide$elapsed_base / wide$elapsed_native
write.csv(wide, file.path(out, "timing-summary.csv"), row.names = FALSE)
write.csv(subset(metrics, mode %in% c("memory", "rss")), file.path(out, "memory.csv"), row.names = FALSE)
writeLines(sub("[[:blank:]]+$", "", readRDS(file.path(raw, "native.multi_equal.10000.0.memory.rds"))$session),
  file.path(out, "session-info.txt"))
cat("All cross-version comparisons passed.\n")
