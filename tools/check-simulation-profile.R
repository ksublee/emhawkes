# Rscript --vanilla tools/check-simulation-profile.R RAW_RESULTS PUBLISHED_RESULTS
# Validate path/RNG identities before comparing timings or diagnostic variants.
a <- commandArgs(TRUE); raw <- a[1]; out <- a[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- lapply(list.files(raw, pattern = "[.]rds$", full.names = TRUE), readRDS)
stopifnot(length(x) > 0)
rows <- lapply(x, function(z) {
  stopifnot(length(z$warnings) == 0)
  data.frame(id = z$id, dimensions = z$meta$d, kernels = z$meta$k, arrival = z$meta$arrival,
    events = z$events, variant = z$variant, mode = z$mode, round = z$round,
    seconds = z$seconds, cpu_seconds = z$cpu_seconds, returned_bytes = z$returned_bytes,
    heap_peak_mib = z$heap_peak_mib, rss_before = z$rss_before, rss_peak = z$rss_peak,
    fingerprint = z$fingerprint, rng = z$rng)
})
all <- do.call(rbind, rows)
checks <- list()
for (key in unique(paste(all$id, all$events))) {
  z <- all[paste(all$id, all$events) == key, ]
  stopifnot(length(unique(z$fingerprint)) == 1L, length(unique(z$rng)) == 1L)
  checks[[length(checks) + 1L]] <- data.frame(id = z$id[1], events = z$events[1],
    samples = nrow(z), reference_available = any(z$variant == "reference"),
    identical_path = TRUE, identical_rng = TRUE)
}
# The parameter-only baseline resolves once. Ignored dots change no model values.
for (ids in list(c("multi_equal", "parameter_mu"), c("state_current", "state_dots"))) {
  for (n in intersect(all$events[all$id == ids[1]], all$events[all$id == ids[2]])) {
    z <- all[all$id %in% ids & all$events == n, ]
    stopifnot(length(unique(z$fingerprint)) == 1L, length(unique(z$rng)) == 1L)
  }
}
write.csv(do.call(rbind, checks), file.path(out, "identity-checks.csv"), row.names = FALSE)
write.csv(all[all$mode == "time", ], file.path(out, "timings.csv"), row.names = FALSE)
timing <- all[all$mode == "time", ]; medians <- list()
for (key in unique(paste(timing$id, timing$events, timing$variant))) {
  z <- timing[paste(timing$id, timing$events, timing$variant) == key, ]
  medians[[length(medians) + 1L]] <- data.frame(id = z$id[1], dimensions = z$dimensions[1], kernels = z$kernels[1],
    arrival = z$arrival[1], events = z$events[1], variant = z$variant[1], rounds = nrow(z),
    median_seconds = median(z$seconds), min_seconds = min(z$seconds), max_seconds = max(z$seconds),
    median_cpu_seconds = median(z$cpu_seconds))
}
write.csv(do.call(rbind, medians), file.path(out, "timing-summary.csv"), row.names = FALSE)
memory <- shapes <- proposals <- phases <- attribution <- list()
for (z in x) {
  meta <- data.frame(id = z$id, dimensions = z$meta$d, kernels = z$meta$k, events = z$events, variant = z$variant)
  if (z$mode %in% c("memory", "peak")) {
    memory[[length(memory) + 1L]] <- cbind(meta, mode = z$mode,
      returned_bytes = z$returned_bytes, allocated_bytes = if (is.null(z$memory)) NA_real_ else z$memory$allocated_bytes,
      largest_block = if (is.null(z$memory)) NA_real_ else z$memory$largest_block,
      rss_before = z$rss_before, rss_peak = z$rss_peak, heap_peak_mib = z$heap_peak_mib)
    if (!is.null(z$memory)) {
      shapes[[length(shapes) + 1L]] <- cbind(meta, z$memory$shapes)
      allocations <- read.csv(file.path(raw, paste(z$id, z$events, z$variant, "allocations.csv", sep = ".")), stringsAsFactors = FALSE)
      stopifnot(sum(allocations$total_bytes) == z$memory$allocated_bytes,
        sum(allocations$count) == z$memory$allocation_count)
      attribution[[length(attribution) + 1L]] <- cbind(meta, head(allocations, 10L))
    }
  }
  if (z$mode == "probe") {
    count <- z$counts; n <- z$events; d <- z$meta$d
    expected_callbacks <- switch(z$id, marked = n, state_current = n, state_dots = n, state_observed = 3 * n + 2, 0)
    stopifnot(count$arrivals == n, count$callbacks == expected_callbacks,
      count$geom_calls == if (z$id %in% c("marked", "state_observed")) n else 0)
    if (z$meta$arrival == "direct") stopifnot(count$exp_calls == n, count$uniform_calls == n,
      count$exp_draws == n * d, count$uniform_draws == n * d)
    else stopifnot(count$exp_calls == count$uniform_calls, count$exp_calls == count$exp_draws,
      count$uniform_calls == count$uniform_draws, count$exp_calls >= n * d)
    proposals[[length(proposals) + 1L]] <- cbind(meta, as.data.frame(count),
      rejection_fraction = if (z$meta$arrival == "thinning") 1 - n * d / count$exp_calls else NA_real_)
  }
  if (z$mode == "profile") {
    prefix <- paste(z$id, z$events, z$variant, sep = ".")
    lines <- readLines(file.path(raw, paste(z$id, z$events, z$variant, "profile", z$round, "Rprof", sep = ".")), warn = FALSE)
    lines <- lines[!grepl("^#File|^GC profiling|^sample[.]interval|^memory", lines)]
    labels <- regmatches(lines, gregexpr('"[^"]+"', lines))
    phase <- vapply(labels, function(names) {
      names <- gsub('"', "", names, fixed = TRUE)
      if ("<GC>" %in% names) "garbage_collection"
      else if (any(names %in% c("handle_arrival", "rarrival"))) "arrival"
      else if ("call_component" %in% names) "callback_dispatch"
      else if ("valid_components" %in% names) "component_validation"
      else "event_loop_and_other"
    }, "")
    counts <- table(phase)
    phases[[length(phases) + 1L]] <- cbind(meta, phase = names(counts), samples = as.numeric(counts), percent = 100 * as.numeric(counts) / sum(counts))
    write.csv(z$profile$by.self, file.path(out, paste0(prefix, ".self.csv")))
    write.csv(z$profile$by.total, file.path(out, paste0(prefix, ".total.csv")))
  }
}
write.csv(do.call(rbind, memory), file.path(out, "memory.csv"), row.names = FALSE)
write.csv(do.call(rbind, shapes), file.path(out, "history-sized-allocations.csv"), row.names = FALSE)
write.csv(do.call(rbind, proposals), file.path(out, "rng-calls.csv"), row.names = FALSE)
write.csv(do.call(rbind, phases), file.path(out, "profile-phases.csv"), row.names = FALSE)
write.csv(do.call(rbind, attribution), file.path(out, "allocation-attribution.csv"), row.names = FALSE)
print(do.call(rbind, medians), row.names = FALSE)
cat("Verified", nrow(all), "samples; identical paths and RNG states within each model/size.\n")
