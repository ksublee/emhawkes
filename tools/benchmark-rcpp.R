# Install the compiled package before running this script.
# Rscript tools/benchmark-rcpp.R <output-directory> [1001,10001,100001] [all|likelihood|fit]
library(emhawkes)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else tempfile("emhawkes-benchmark-")
sizes <- if (length(args) > 1L) as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]]) else c(1001L, 10001L, 100001L)
selection <- if (length(args) > 2L) match.arg(args[3], c("all", "likelihood", "fit")) else "all"
if (!length(sizes) || anyNA(sizes) || any(sizes < 2L)) stop("Row counts must include an initial row and at least one event.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)

# Batch fast calls so the clock's resolution cannot manufacture a speedup.
timed_batch <- function(fun, count) {
  unname(system.time(for (i in seq_len(count)) fun())[["elapsed"]])
}
repetitions <- function(fun) {
  count <- 1L
  repeat {
    elapsed <- timed_batch(fun, count)
    if (elapsed >= .2 || count >= 65536L) return(count)
    count <- count * 2L
  }
}
measure_pair <- function(r, cpp) {
  r(); cpp()
  counts <- c(r = repetitions(r), cpp = repetitions(cpp))
  elapsed <- matrix(0, 3, 2, dimnames = list(NULL, c("r", "cpp")))
  for (round in 1:3) {
    for (engine in if (round %% 2) c("r", "cpp") else c("cpp", "r")) {
      gc()
      fun <- if (engine == "r") r else cpp
      elapsed[round, engine] <- timed_batch(fun, counts[engine]) / counts[engine]
    }
  }
  c(r_seconds = median(elapsed[, "r"]), cpp_seconds = median(elapsed[, "cpp"]),
    speedup = median(elapsed[, "r"]) / median(elapsed[, "cpp"]))
}

constant <- new("hspec", mu = c(.2, .3), alpha = matrix(c(.2, .1, .1, .2), 2), beta = matrix(1.2, 2, 2))
multi <- new("hspec", mu = c(.2, .3), alpha = matrix(.08, 2, 4),
             beta = matrix(rep(c(.8, 1.2, 1.6, 2), each = 2), 2),
             type_col_map = list(c(1, 3), c(2, 4)))
cases <- list(univariate = new("hspec", mu = .2, alpha = .2, beta = 1.2),
              constant = constant, multi_kernel = multi)
results <- list()
if (selection != "fit") for (size in sizes) {
  dt <- c(0, rep(c(.1, .2, .3), length.out = size - 1L))
  marks <- c(0, rep(c(1, 2, 1), length.out = size - 1L))
  for (name in names(cases)) {
    object <- cases[[name]]
    types <- c(0, rep(seq_len(object@dimens), length.out = size - 1L))
    observation <- emhawkes:::normalize_observation(object, dt, types, marks)
    callbacks <- emhawkes:::model_callback_plans(object)
    for (mode in c("prepared_loglik", "public_loglik", "inference")) {
      make_evaluation <- function(engine) {
        force(engine)
        if (mode == "public_loglik") {
          function() {
            previous <- options(emhawkes.engine = engine)
            on.exit(options(previous))
            logLik(object, dt, types, marks, lambda_component0 = .1)
          }
        } else {
          function() emhawkes:::loglik_hspec(object, observation, .1,
            infer = mode == "inference", callbacks = callbacks,
            show_initial_message = FALSE, engine = engine)
        }
      }
      r <- make_evaluation("r")
      cpp <- make_evaluation("cpp")
      reference <- r()
      candidate <- cpp()
      stopifnot(isTRUE(all.equal(reference, candidate, tolerance = 1e-10)))
      maximum_difference <- max(abs(unlist(reference, use.names = FALSE) -
                                    unlist(candidate, use.names = FALSE)))
      rm(reference, candidate)
      measured <- measure_pair(r, cpp)
      row <- data.frame(case = name, rows = size, dimensions = object@dimens,
                        kernels = ncol(object@beta), operation = mode,
                        t(measured), max_abs_difference = maximum_difference, row.names = NULL)
      results[[length(results) + 1L]] <- row
      print(row, row.names = FALSE)
    }
  }
}
if (selection != "fit") write.csv(do.call(rbind, results), file.path(output, "likelihood-timings.csv"), row.names = FALSE)
if (selection == "likelihood") {
  capture.output(sessionInfo(), file = file.path(output, "session-info.txt"))
  cat("Benchmark results:", output, "\n")
  quit(status = 0)
}

# A full fit includes constructor, shared-parameter reconstruction, optimizer,
# numerical derivatives, and the native boundary, not just the event loop.
set.seed(741)
truth <- new("hspec", mu = .35, alpha = .3, beta = 1.5)
fits <- list(univariate = list(path = hsim(truth, size = 2001, lambda_component0 = .1),
  start = new("hspec", mu = .25, alpha = .2, beta = 1.2), fixed = NULL))
set.seed(742)
truth <- new("hspec", mu = c(.25, .35), alpha = matrix(.06, 2, 4),
  beta = matrix(rep(c(.8, 1.2, 1.6, 2), each = 2), 2), type_col_map = list(c(1, 3), c(2, 4)))
path <- hsim(truth, size = 2001, lambda_component0 = .1)
start <- truth
start@mu <- matrix(c(.2, .3), 2)
start@alpha[] <- .04
fits$multi_kernel <- list(path = path, start = start, fixed = c(rep(FALSE, 3), rep(TRUE, 4)))
# Count engine evaluations in separate, untimed fits; tracing must not affect
# the speed measurements. The counter excludes proposals rejected before LL.
count_evaluations <- function(engine) {
  counter <- new.env(parent = emptyenv())
  counter$value <- 0L
  tracer <- substitute(assign("value", get("value", envir = COUNTER, inherits = FALSE) + 1L,
                              envir = COUNTER), list(COUNTER = counter))
  namespace <- asNamespace("emhawkes")
  trace("loglik_hspec", tracer = tracer, print = FALSE, where = namespace)
  on.exit(untrace("loglik_hspec", where = namespace))
  fit(engine)
  counter$value
}
fit_results <- list()
for (name in names(fits)) {
  specification <- fits[[name]]
  parameter_count <- length(emhawkes:::build_parameter_map(specification$start)$start)
  fit <- function(engine) {
    previous <- options(emhawkes.engine = engine)
    on.exit(options(previous))
    hfit(specification$start, specification$path$inter_arrival, specification$path$type,
      lambda_component0 = .1, fixed = specification$fixed,
      constraints = list(ineqA = diag(parameter_count), ineqB = rep(-1e-6, parameter_count)),
      control = list(iterlim = 200, reltol = 1e-9))
  }
  r_fit <- suppressMessages(fit("r"))
  cpp_fit <- suppressMessages(fit("cpp"))
  stopifnot(r_fit$code == 0, cpp_fit$code == 0,
            abs(as.numeric(logLik(r_fit)) - as.numeric(logLik(cpp_fit))) < 1e-7)
  row <- data.frame(case = name, rows = 2001, operation = "hfit",
    t(suppressMessages(measure_pair(function() fit("r"), function() fit("cpp")))), row.names = NULL)
  row$r_evaluations <- suppressMessages(count_evaluations("r"))
  row$cpp_evaluations <- suppressMessages(count_evaluations("cpp"))
  row$loglik_abs_difference <- abs(as.numeric(logLik(r_fit)) - as.numeric(logLik(cpp_fit)))
  row$coefficient_max_abs_difference <- max(abs(coef(r_fit) - coef(cpp_fit)))
  fit_results[[length(fit_results) + 1L]] <- row
}
fit_results <- do.call(rbind, fit_results)
print(fit_results, row.names = FALSE)
write.csv(fit_results, file.path(output, "fit-timings.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "session-info.txt"))
cat("Benchmark results:", output, "\n")
