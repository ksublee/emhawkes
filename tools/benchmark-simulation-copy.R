# Rscript --vanilla tools/benchmark-simulation-copy.R LIB OUT LABEL CASE EVENTS [memory]
# Run versions in separate, sequential processes. The case factory is shared.
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(a[1])))
source("tools/simulation-profile-cases.R")
x <- simulation_profile_case(a[4]); events <- as.integer(a[5])
dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
run <- function(n = events) hsim(x$object, size = n + 1L, lambda_component0 = x$initial)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (i in 1:3) { set.seed(x$seed); invisible(run(100L)) }
rows <- list(); reference <- seed <- NULL
for (i in 1:3) {
  gc(); set.seed(x$seed)
  timing <- system.time(path <- run())
  data <- path[setdiff(names(path), "hspec")]
  if (i == 1) { reference <- data; seed <- .Random.seed } else
    stopifnot(identical(data, reference), identical(seed, .Random.seed))
  rows[[i]] <- data.frame(version = a[3], case = a[4], events, round = i,
    elapsed = unname(timing["elapsed"]), cpu = unname(timing["user.self"] + timing["sys.self"]))
}
allocated <- NA_real_
if (length(a) >= 6L && a[6] == "memory") {
  gc(); set.seed(x$seed); profile <- tempfile(fileext = ".Rprofmem")
  Rprofmem(profile); path <- tryCatch(run(), finally = Rprofmem(NULL))
  stopifnot(identical(reference, path[setdiff(names(path), "hspec")]), identical(seed, .Random.seed))
  con <- file(profile, "r"); allocated <- 0
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE); if (!length(lines)) break
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    allocated <- allocated + sum(bytes[is.finite(bytes)])
  }
  close(con); unlink(profile)
}
saveRDS(list(timings = do.call(rbind, rows), allocated_bytes = allocated,
  output_bytes = as.numeric(object.size(path)), data = reference, seed = seed,
  session = capture.output(sessionInfo())),
  file.path(a[2], paste(a[3], a[4], events, "rds", sep = ".")))
cat(a[3], a[4], events, "median_seconds", median(vapply(rows, function(r) r$elapsed, 0)), "\n")
