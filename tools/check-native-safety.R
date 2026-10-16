# Rscript --vanilla tools/check-native-safety.R basic|gc OUT [LIB]
a <- commandArgs(TRUE); stopifnot(length(a) %in% c(2L, 3L), a[1] %in% c("basic", "gc"))
if (length(a) == 3L) .libPaths(c(a[3], .libPaths()))
library(emhawkes)
source("tools/native-safety-cases.R")
dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "cpp")
fixtures <- native_safety_cases(); rows <- list()
run_gc <- function(fun, step) {
  on.exit(gctorture2(0), add = TRUE)
  if (step) gctorture2(step, wait = step, inhibit_release = TRUE)
  value <- fun()
  active <- gctorture2(0)
  stopifnot(active == step)
  list(value = value, active = active)
}
for (name in names(fixtures)) {
  set.seed(7361)
  expected <- fixtures[[name]](); expected_rng <- .Random.seed
  set.seed(7361)
  step <- if (a[1] == "gc") if (grepl("likelihood", name)) 10L else 1L else 0L
  actual <- run_gc(fixtures[[name]], step)
  stopifnot(identical(expected, actual$value), identical(expected_rng, .Random.seed))
  rows[[length(rows) + 1L]] <- data.frame(case = name, gc_step = step,
    active_gc_step = actual$active, output_identical = TRUE, rng_identical = TRUE)
  cat(name, "passed; GC step", step, "\n"); flush.console()
}
source("tests/native-unwind.R")
cat("Callback errors/interrupts and subsequent RNG scopes passed.\n")
write.csv(do.call(rbind, rows), file.path(a[2], "checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(a[2], "session-info.txt"))
