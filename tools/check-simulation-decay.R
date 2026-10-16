# Rscript --vanilla tools/check-simulation-decay.R LIB OUT LABEL
a <- commandArgs(TRUE); stopifnot(length(a) == 3L)
.libPaths(c(a[1], .libPaths())); library(emhawkes)
source("tools/simulation-decay-cases.R")
ids <- c("univariate", "multi_equal", "large_equal", "signed", "state_current",
         "mark_only", "marked", "heavy_impact", "wide32", "wide128", "impact32")
results <- list()
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
  RNGkind(kind, "Inversion", "Rejection")
  for (id in ids) {
    x <- simulation_profile_case(id)
    for (size in c(1L, 2L, 31L, 251L)) for (seed in c(1L, 62L, 7361L)) {
      set.seed(seed)
      path <- hsim(x$object, size, lambda_component0 = x$initial, N0 = 2, Nc0 = 5)
      results[[length(results) + 1L]] <- list(case = id, kind = kind, size = size,
        seed = seed, data = path[setdiff(names(path), "hspec")], rng = .Random.seed)
    }
  }
}
saveRDS(results, file.path(a[2], paste0(a[3], ".decay-checks.rds")), compress = FALSE)
cat(length(results), "cross-version fixtures recorded\n")
