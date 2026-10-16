# Rscript --vanilla tools/run-simulation-copy.R OLD_LIB NEW_LIB WORK OUT
# WORK stores raw RDS paths and logs; OUT stores the small, shareable summaries.
a <- commandArgs(TRUE); stopifnot(length(a) == 4L)
libs <- setNames(normalizePath(a[1:2], winslash = "/", mustWork = TRUE), c("old", "new"))
work <- a[3]; out <- a[4]
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(script, arguments, log) {
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(arguments)),
    stdout = log, stderr = log)
  if (status != 0L) stop("Worker failed: ", log)
}
for (label in names(libs)) run("tools/check-simulation-copy.R", c(libs[[label]], work, label),
  file.path(work, paste0(label, ".checks.log")))
cases <- c("univariate", "multivariate", "multi_equal", "multi_unequal", "large_equal",
  "large_unequal", "marked", "parameter_mu", "state_current", "state_dots", "state_observed", "signed")
design <- rbind(data.frame(case = cases, events = 3000L),
  data.frame(case = c("univariate", "multivariate", "multi_equal", "multi_unequal", "marked",
    "state_current", "state_dots", "state_observed"), events = 10000L))
for (i in seq_len(nrow(design))) for (label in if (i %% 2L) c("old", "new") else c("new", "old")) {
  z <- design[i, ]; key <- paste(label, z$case, z$events, sep = ".")
  run("tools/benchmark-simulation-copy.R", c(libs[[label]], work, label, z$case, z$events,
    if (z$events == 3000L) "memory" else character()), file.path(work, paste0(key, ".log")))
  cat(key, "completed\n")
}
run("tools/summarize-simulation-copy.R", c(work, out, libs[["new"]]), file.path(work, "summary.log"))
installed <- do.call(rbind, lapply(names(libs), function(label) {
  pkg <- file.path(libs[[label]], "emhawkes")
  paths <- c(file.path(pkg, "DESCRIPTION"), list.files(file.path(pkg, "R"), full.names = TRUE),
    list.files(file.path(pkg, "libs"), recursive = TRUE, full.names = TRUE))
  paths <- paths[file.info(paths)$isdir == FALSE]
  data.frame(version = label, file = substring(paths, nchar(pkg) + 2L),
    md5 = unname(tools::md5sum(paths)))
}))
write.csv(installed, file.path(out, "installed-code.csv"), row.names = FALSE)
cat("All path/RNG comparisons and inference replay passed.\n")
