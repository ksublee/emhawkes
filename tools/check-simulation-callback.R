# Rscript --vanilla tools/check-simulation-callback.R NEW_LIB OUT
a <- commandArgs(TRUE); stopifnot(length(a) == 2L)
.libPaths(c(a[1], .libPaths())); library(emhawkes)
source("tools/simulation-callback-cases.R")
reference <- getFromNamespace("simulate_hawkes_r", "emhawkes")
rows <- list()
for (id in c("parameter_mu", "state_current", "mark_only", "marked", "eta_only", "heavy_impact")) {
  x <- simulation_profile_case(id)
  set.seed(x$seed)
  path <- hsim(x$object, 1001, lambda_component0 = x$initial, N0 = 2, Nc0 = 5)
  rng <- .Random.seed
  set.seed(x$seed)
  old <- reference(x$object, 1001, lambda_component0 = x$initial, N0 = 2, Nc0 = 5)
  stopifnot(identical(old, path), identical(rng, .Random.seed))
  replays <- list()
  for (engine in c("r", "cpp")) {
    op <- options(emhawkes.engine = engine)
    replays[[engine]] <- logLik(x$object, path$inter_arrival, path$type, path$mark,
      lambda_component0 = x$initial, N0 = 2, Nc0 = 5, infer = TRUE)
    options(op)
    for (name in c("lambda", "rambda", "lambda_component", "rambda_component")) {
      error <- max(abs(path[[name]] - replays[[engine]][[name]]))
      stopifnot(error == 0)
      rows[[length(rows) + 1L]] <- data.frame(case = id, engine, field = name, absolute_error = error)
    }
    stopifnot(identical(rng, .Random.seed))
  }
  error <- abs(replays$r$loglikelihood - replays$cpp$loglikelihood) / 1000
  stopifnot(error <= 1e-12)
  rows[[length(rows) + 1L]] <- data.frame(case = id, engine = "r_cpp",
    field = "loglikelihood_per_event", absolute_error = error)
}
write.csv(do.call(rbind, rows), file.path(a[2], "replay.csv"), row.names = FALSE)
cat("Six independent replay fixtures passed without consuming RNG.\n")
