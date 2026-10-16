# POSIX SIGINT during the native thinning loop; bounded child process and replay.
# Rscript --vanilla tools/check-native-interrupt.R OUT [LIB]
a <- commandArgs(TRUE)
worker <- length(a) && a[1] == "--worker"
if (worker) a <- a[-1]
stopifnot(length(a) %in% c(1L, 2L))
if (length(a) == 2L) .libPaths(c(a[2], .libPaths()))
library(emhawkes)
if (.Platform$OS.type != "unix") stop("This OS-signal check requires POSIX signals.")

if (worker) {
  native <- getFromNamespace("hawkes_simulate_cpp", "emhawkes")
  mu <- 1e-30; initial <- matrix(5e7, 1, 2)
  beta <- matrix(c(1e30, 2e30), 1)
  # Rapid decay makes the first acceptance effectively impossible. Output is
  # only two rows; the long work is in the interruptible rejection loop.
  RNGkind("Wichmann-Hill", "Inversion", "Rejection"); set.seed(7361)
  before <- .Random.seed
  writeLines(as.character(Sys.getpid()), file.path(a[1], "ready.tmp"))
  stopifnot(file.rename(file.path(a[1], "ready.tmp"), file.path(a[1], "ready.pid")))
  condition <- tryCatch(native(mu, matrix(0, 1, 2), beta, initial, list(0:1),
    matrix(0, 2, 1), matrix(0, 2, 1), numeric(2), numeric(2), c(0, 1)),
    interrupt = function(e) e, error = function(e) e)
  after <- .Random.seed
  draws <- runif(3)
  recovery <- hsim(new("hspec", mu = .3, alpha = 0, beta = 1), 5, lambda_component0 = 0)
  saveRDS(list(interrupted = inherits(condition, "interrupt"), before = before, after = after,
    draws = draws, recovery = recovery, final_rng = .Random.seed), file.path(a[1], "result.tmp"))
  stopifnot(file.rename(file.path(a[1], "result.tmp"), file.path(a[1], "result.rds")))
  quit(save = "no", status = 0L)
}

dir.create(a[1], recursive = TRUE, showWarnings = FALSE)
out <- tempfile("signal-", tmpdir = normalizePath(a[1], winslash = "/")); dir.create(out)
script <- normalizePath("tools/check-native-interrupt.R", winslash = "/")
rscript <- file.path(R.home("bin"), "Rscript")
launcher <- normalizePath("tools/run-native-signal.py", winslash = "/")
status <- system2("python3", shQuote(c(launcher, rscript, script, out,
  if (length(a) == 2L) a[2])))
if (status != 0L) stop("Native signal worker failed; inspect ", file.path(out, "worker.log"))
result_file <- file.path(out, "result.rds")
result <- readRDS(result_file)
stopifnot(result$interrupted)

RNGkind("Wichmann-Hill", "Inversion", "Rejection"); set.seed(7361)
stopifnot(identical(.Random.seed, result$before))
matched <- NA_integer_
for (attempt in seq_len(4194303L)) {
  rexp(1, rate = 1e8); runif(1, 0, 1e8)
  # checkUserInterrupt occurs before draws at attempts 1024, 2048, ... .
  if (attempt %% 1024L == 1023L && identical(.Random.seed, result$after)) {
    matched <- attempt; break
  }
}
stopifnot(!is.na(matched), identical(runif(3), result$draws))
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")
expected <- r_sim(new("hspec", mu = .3, alpha = 0, beta = 1), 5, lambda_component0 = 0)
stopifnot(identical(expected, result$recovery), identical(.Random.seed, result$final_rng))
write.csv(data.frame(interrupted = TRUE, rejected_pairs = matched, rng_replayed = TRUE,
  recovery_identical = TRUE), file.path(a[1], "interrupt.csv"), row.names = FALSE)
cat("Native SIGINT unwound after", matched, "rejected pairs; RNG and recovery match R.\n")
