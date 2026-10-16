# Rscript --vanilla tools/check-simulation-replay.R LIB OUT HSIM_SOURCE
args <- commandArgs(TRUE); .libPaths(c(args[1], .libPaths())); library(emhawkes)
source("tools/simulation-profile-cases.R"); out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cases <- c("univariate", "multivariate", "multi_equal", "multi_unequal", "large_equal", "large_unequal",
  "marked", "parameter_mu", "state_current", "state_dots", "state_observed", "signed")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
variant <- simulation_profile_variant("both", args[3]); rows <- list()
fields <- c("lambda", "lambda_component", "rambda", "rambda_component")
for (id in cases) {
  x <- simulation_profile_case(id); set.seed(x$seed)
  path <- hsim(x$object, size = 1001, lambda_component0 = x$initial); seed <- .Random.seed
  set.seed(x$seed); other <- variant(x$object, size = 1001, lambda_component0 = x$initial)
  stopifnot(identical(seed, .Random.seed), identical(path[setdiff(names(path), "hspec")], other[setdiff(names(other), "hspec")]))
  for (i in seq_len(x$meta$d)) {
    stopifnot(identical(as.numeric(path$N[, i]), cumsum(as.numeric(path$type == i))),
      identical(as.numeric(path$Nc[, i]), cumsum(path$mark * (path$type == i))))
  }
  replay <- list()
  for (engine in c("r", "cpp")) {
    options(emhawkes.engine = engine)
    replay[[engine]] <- suppressMessages(logLik(x$object, path$inter_arrival, path$type, path$mark,
      lambda_component0 = x$initial, infer = TRUE))
    stopifnot(identical(seed, .Random.seed))
    for (field in fields) {
      reference <- path[[field]]; actual <- replay[[engine]][[field]]
      stopifnot(identical(dim(actual), dim(reference)))
      error <- max(abs(actual - reference) / pmax(1, abs(reference)))
      stopifnot(is.finite(error), error < 1e-11)
      rows[[length(rows) + 1L]] <- data.frame(id, events = 1000L, engine, field, scaled_error = error)
    }
  }
  difference <- abs(replay$r$loglikelihood - replay$cpp$loglikelihood)
  stopifnot(difference / 1000 < 1e-11)
  rows[[length(rows) + 1L]] <- data.frame(id, events = 1000L, engine = "r_cpp", field = "loglikelihood_per_event", scaled_error = difference / 1000)
}
write.csv(do.call(rbind, rows), file.path(out, "replay-checks.csv"), row.names = FALSE)
# Verify the 100k diagnostic prefix against the public simulator at 10k.
prefixes <- list()
for (id in c("large_equal", "state_current")) {
  x <- simulation_profile_case(id); set.seed(x$seed)
  reference <- hsim(x$object, size = if (id == "large_equal") 3001L else 10001L, lambda_component0 = x$initial)
  set.seed(x$seed); large <- variant(x$object, size = 100001L, lambda_component0 = x$initial)
  n <- length(reference$inter_arrival)
  for (field in setdiff(names(reference), "hspec")) {
    actual <- if (is.matrix(large[[field]])) large[[field]][seq_len(n), , drop = FALSE] else large[[field]][seq_len(n)]
    stopifnot(identical(actual, reference[[field]]))
  }
  prefixes[[length(prefixes) + 1L]] <- data.frame(id, events = 100000L, checked_prefix_events = n - 1L, identical = TRUE)
}
write.csv(do.call(rbind, prefixes), file.path(out, "large-prefix-checks.csv"), row.names = FALSE)
cat("Independent R/C++ intensity replay, counts and large-prefix checks passed.\n")
