# Rscript --vanilla tools/benchmark-auto-engine.R LIB BID_ASK_FIXTURES IMPACT_FIXTURES OUT
# Fixture directories must include fixtures.rds and the previous cpp-CASE.rds fits.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
rows <- list()
for (group in c('bid_ask', 'impact')) {
  original <- args[if (group == 'bid_ask') 2L else 3L]
  out <- file.path(args[4], group)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(original, 'fixtures.rds'), out, overwrite = TRUE))
  script <- if (group == 'bid_ask') 'tools/benchmark-bid-ask-numeric.R' else 'tools/benchmark-impact-numeric.R'
  for (case in names(readRDS(file.path(out, 'fixtures.rds')))) {
    status <- system2(file.path(R.home('bin'), 'Rscript.exe'),
      c('--vanilla', shQuote(script), 'worker', shQuote(args[1]), 'auto', case, '1', shQuote(out)))
    stopifnot(status == 0L)
    value <- readRDS(file.path(out, paste0('auto-', case, '.rds')))
    previous <- readRDS(file.path(original, paste0('cpp-', case, '.rds')))
    stopifnot(identical(value$results[[1]], previous$results[[1]]))
    rows[[case]] <- data.frame(case, events = value$events,
      auto_seconds = value$times[1], previous_cpp_identical = TRUE)
  }
}
write.csv(do.call(rbind, rows), file.path(args[4], 'comparison.csv'), row.names = FALSE)
