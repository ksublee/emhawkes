# Rscript --vanilla tools/run-simulation-callback.R BASE_LIB NEW_LIB RAW OUT
# Use baseline main 36e79e5. Run alone after compilation and tests have finished.
a <- commandArgs(TRUE); stopifnot(length(a) == 4L)
libs <- setNames(normalizePath(a[1:2], winslash = "/", mustWork = TRUE), c("base", "native"))
raw <- a[3]; out <- a[4]
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(script, args, log) {
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(args)), stdout = log, stderr = log)
  if (status != 0) stop("Worker failed: ", log)
}
for (label in names(libs)) run("tools/check-simulation-copy.R", c(libs[[label]], raw, label),
  file.path(raw, paste0(label, ".checks.log")))
stopifnot(identical(readRDS(file.path(raw, "base.checks.rds")),
                    readRDS(file.path(raw, "native.checks.rds"))))
design <- rbind(
  data.frame(case = "state_current", events = c(0L, 10L, 100L, 10000L, 100000L)),
  data.frame(case = "mark_only", events = c(100L, 10000L, 100000L)),
  data.frame(case = c("parameter_mu", "marked", "eta_only", "heavy_impact", "multi_equal",
                     "multi_unequal", "state_dots", "state_observed"), events = 10000L),
  data.frame(case = "multi_equal", events = c(0L, 10L, 100L)))
control <- design$case %in% c("multi_equal", "multi_unequal")
design <- rbind(design[control, ], design[!control, ])
write.csv(design, file.path(out, "design.csv"), row.names = FALSE)
metrics <- list()
pair <- function(id, events, round, mode, order) {
  values <- list()
  for (label in order) {
    key <- paste(label, id, events, round, mode, sep = ".")
    run("tools/benchmark-native-simulation.R", c(libs[[label]], raw, label, id, events, round,
      mode, "tools/simulation-callback-cases.R"), file.path(raw, paste0(key, ".log")))
    values[[label]] <- readRDS(file.path(raw, paste0(key, ".rds")))
    metrics[[length(metrics) + 1L]] <<- values[[label]]$metrics
  }
  stopifnot(identical(values$base$data, values$native$data), identical(values$base$seed, values$native$seed))
  write.csv(do.call(rbind, metrics), file.path(out, "measurements.csv"), row.names = FALSE)
  cat(id, events, round, mode, "exact output/RNG match\n"); flush.console()
}
for (i in seq_len(nrow(design))) for (round in 1:3) {
  order <- if ((i + round) %% 2L) c("base", "native") else c("native", "base")
  pair(design$case[i], design$events[i], round, "timing", order)
  gc()
}
for (mode in c("memory", "rss")) for (id in c("state_current", "mark_only", "marked",
  "heavy_impact", "multi_equal", "state_dots")) {
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
writeLines(sub("[[:blank:]]+$", "", readRDS(file.path(raw, "native.state_current.10000.0.memory.rds"))$session),
  file.path(out, "session-info.txt"))
cat("All cross-version comparisons passed.\n")
