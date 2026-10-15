# Run from the repository root, using an isolated installation of the tested commit.
# Rscript --vanilla tools/benchmark-auto-r.R prepare|run|compare LIB OUT
# Optional argument after OUT for run: repetitions (default 3).
# Internal worker: worker LIB OUT CASE ENGINE ROUND
args <- commandArgs(TRUE)
mode <- args[1]
lib <- args[2]
out <- args[3]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(emhawkes))
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
launch <- function(arguments) {
  status <- system2(rscript, c("--vanilla", vapply(arguments, shQuote, "")))
  stopifnot(status == 0L)
}
if (mode == "prepare") {
  stopifnot(!length(list.files(out, pattern = "[.](r|auto)[.][0-9]+[.]rds$")))
  options(emhawkes.engine = "auto")
  fixtures <- list()
  specifications <- list(
    univariate = list(seed = 9511L,
      truth = new("hspec", mu = .3, alpha = .2, beta = 1.2),
      start = new("hspec", mu = .25, alpha = .15, beta = 1)),
    multi_kernel = list(seed = 9512L,
      truth = new("hspec", mu = .3, alpha = matrix(c(.2, .8), 1),
        beta = matrix(c(1, 5), 1), type_col_map = list(1:2)),
      start = new("hspec", mu = .25, alpha = matrix(c(.15, .6), 1),
        beta = matrix(c(.8, 4), 1), type_col_map = list(1:2))))
  for (name in names(specifications)) {
    spec <- specifications[[name]]
    set.seed(spec$seed)
    data <- hsim(spec$truth, size = 2001L, lambda_component0 = .1)
    count <- length(emhawkes:::build_parameter_map(spec$start, TRUE)$start)
    fixtures[[name]] <- list(events = 2000L, arguments = list(object = spec$start,
      inter_arrival = data$inter_arrival, type = data$type, mark = data$mark,
      lambda_component0 = .1, method = "BFGS",
      constraints = list(ineqA = diag(count), ineqB = rep(-1e-5, count)),
      control = list(iterlim = 500, reltol = 1e-8)))
  }
  for (group in c("impact", "bid-ask")) {
    destination <- file.path(out, paste0("prepare-", group))
    launch(c(paste0("tools/benchmark-", group, "-numeric.R"), "prepare", lib, destination))
    values <- readRDS(file.path(destination, "fixtures.rds"))
    if (group == "impact") fixtures$flocking <- values$flocking else
      fixtures$bid_ask_current <- values$bid_ask_current
  }
  saveRDS(fixtures, file.path(out, "fixtures.rds"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
  quit()
}
fixtures <- readRDS(file.path(out, "fixtures.rds"))
if (mode == "run") {
  repetitions <- if (length(args) > 3L) as.integer(args[4]) else 3L
  for (case in names(fixtures)) for (round in seq_len(repetitions)) {
    engines <- if (round %% 2L) c("r", "auto") else c("auto", "r")
    for (engine in engines) {
      result_file <- file.path(out, paste(case, engine, round, "rds", sep = "."))
      if (file.exists(result_file)) {
        previous <- readRDS(result_file)
        stopifnot(identical(previous$affinity, Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")))
        next # Resume only with the same CPU policy.
      }
      cat("Starting", case, engine, "round", round, "\n"); flush.console()
      launch(c("tools/benchmark-auto-r.R", "worker", lib, out, case, engine, round))
    }
  }
  quit()
}
if (mode == "worker") {
  case <- args[4]; engine <- args[5]; round <- as.integer(args[6])
  options(emhawkes.engine = engine)
  affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
  if (affinity != "scheduler") {
    if (.Platform$OS.type != "windows") stop("Benchmark affinity is supported on Windows.")
    stopifnot(grepl("^[1-9][0-9]*$", affinity))
    command <- sprintf(paste0("$taskProcess = Get-Process -Id %d; ",
      "$taskProcess.ProcessorAffinity = [IntPtr]%s; ",
      "$taskProcess.ProcessorAffinity.ToInt64()"), Sys.getpid(), affinity)
    actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command",
      shQuote(command)), stdout = TRUE)
    stopifnot(length(actual) == 1L, trimws(actual) == affinity)
  }
  a <- fixtures[[case]]$arguments
  warnings <- character()
  fit_once <- function(arguments) withCallingHandlers(
    suppressMessages(do.call(hfit, arguments)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  warm <- a
  size <- min(101L, length(a$inter_arrival))
  for (name in intersect(names(a), c("inter_arrival", "type", "mark")))
    warm[[name]] <- a[[name]][seq_len(size)]
  for (name in intersect(names(a), c("N", "Nc")))
    warm[[name]] <- a[[name]][seq_len(size), , drop = FALSE]
  warm$control <- list(iterlim = 2L)
  invisible(fit_once(warm))
  warnings <- character()
  gc()
  time <- system.time(value <- fit_once(a))
  elapsed <- unname(time[["elapsed"]])
  stopifnot(attr(value, "emhawkes")$gradient == "numeric")
  fields <- value[c("estimate", "maximum", "gradient", "hessian", "code", "iterations")]
  saveRDS(list(case = case, engine = engine, round = round, elapsed = elapsed,
    result = fields, warnings = unique(warnings), affinity = affinity,
    cpu_user = unname(time[["user.self"]]), cpu_system = unname(time[["sys.self"]])),
    file.path(out, paste(case, engine, round, "rds", sep = ".")))
  cat(case, engine, round, elapsed, "seconds; code", value$code, "\n")
  quit()
}
stopifnot(mode == "compare")
files <- list.files(out, pattern = "[.](r|auto)[.][0-9]+[.]rds$", full.names = TRUE)
values <- lapply(files, readRDS)
metrics <- summary <- estimates <- timings <- list()
scaled <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
for (case in names(fixtures)) {
  fits <- Filter(function(x) x$case == case, values)
  rounds <- sort(unique(vapply(fits, `[[`, 0L, "round")))
  stopifnot(length(rounds) >= 3L)
  reference <- Filter(function(x) x$engine == "r", fits)[[1]]$result
  for (round in rounds) {
    pair <- Filter(function(x) x$round == round, fits)
    stopifnot(length(pair) == 2L)
    r <- Filter(function(x) x$engine == "r", pair)[[1]]
    auto <- Filter(function(x) x$engine == "auto", pair)[[1]]
    for (field in names(reference)) {
      a <- auto$result[[field]]; b <- r$result[[field]]
      metrics[[length(metrics) + 1L]] <- data.frame(case, round, field,
        identical = identical(a, b), max_abs_difference = max(abs(a - b)),
        max_scaled_difference = scaled(a, b))
    }
    stopifnot(identical(names(auto$result$estimate), names(r$result$estimate)),
      scaled(auto$result$estimate, r$result$estimate) < 1e-4,
      abs(auto$result$maximum - r$result$maximum) < 1e-7,
      identical(auto$result$code, r$result$code), auto$result$code == 0L)
    for (fit in pair) {
      stopifnot(identical(fit$result, Filter(function(x) x$engine == fit$engine, fits)[[1]]$result))
      timings[[length(timings) + 1L]] <- data.frame(case, engine = fit$engine,
        round, events = fixtures[[case]]$events, seconds = fit$elapsed,
        cpu_user = if (is.null(fit$cpu_user)) NA_real_ else fit$cpu_user,
        cpu_system = if (is.null(fit$cpu_system)) NA_real_ else fit$cpu_system,
        affinity = if (is.null(fit$affinity)) "scheduler" else fit$affinity,
        warnings = paste(fit$warnings, collapse = " | "))
      estimates[[length(estimates) + 1L]] <- data.frame(case, engine = fit$engine,
        round, parameter = names(fit$result$estimate), estimate = unname(fit$result$estimate))
    }
  }
  times <- do.call(rbind, timings)
  times <- times[times$case == case, ]
  r_time <- median(times$seconds[times$engine == "r"])
  auto_time <- median(times$seconds[times$engine == "auto"])
  summary[[case]] <- data.frame(case, events = fixtures[[case]]$events,
    repetitions = length(rounds), r_seconds = r_time, auto_seconds = auto_time,
    speedup = r_time / auto_time)
}
# Compare the two likelihood implementations at common parameter points.
points <- list()
for (case in names(fixtures)) {
  a <- fixtures[[case]]$arguments
  map <- emhawkes:::build_parameter_map(a$object, TRUE)
  observation <- do.call(emhawkes:::normalize_observation,
    a[intersect(names(a), c("object", "inter_arrival", "type", "mark", "N", "Nc", "N0", "Nc0"))])
  objective <- emhawkes:::make_hawkes_objective(a$object, observation,
    a$lambda_component0, map, derivatives = FALSE)
  fits <- Filter(function(x) x$case == case && x$round == 1L, values)
  parameters <- c(list(start = map$start), setNames(lapply(fits, function(x) x$result$estimate),
    vapply(fits, `[[`, "", "engine")))
  for (point in names(parameters)) {
    param <- parameters[[point]]
    model <- a$object
    components <- emhawkes:::restore_parameters(map, param)
    for (name in names(components)) slot(model, name) <- components[[name]]
    options(emhawkes.engine = "r")
    likelihood_args <- a[intersect(names(a), c("object", "inter_arrival", "type", "mark", "N", "Nc", "N0", "Nc0", "lambda_component0"))]
    likelihood_args$object <- model
    r_value <- suppressMessages(sum(do.call(logLik, likelihood_args)))
    cpp_value <- objective$loglik(param)
    stopifnot(abs(r_value - cpp_value) < 1e-9)
    points[[length(points) + 1L]] <- data.frame(case, point, r_loglik = r_value,
      cpp_loglik = cpp_value, absolute_difference = abs(r_value - cpp_value))
  }
}
write.csv(do.call(rbind, points), file.path(out, "common-points.csv"), row.names = FALSE)
for (name in c("metrics", "summary", "estimates", "timings"))
  write.csv(do.call(rbind, get(name)), file.path(out, paste0(name, ".csv")), row.names = FALSE)
print(do.call(rbind, summary), row.names = FALSE)
