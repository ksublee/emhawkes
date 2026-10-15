# From the repository root, with an isolated installation of the tested source:
# Rscript --vanilla tools/benchmark-scaling.R prepare|run|compare|snapshots LIB OUT
# run uses three timing rounds; memory is measured in separate processes.
# Optional: EMHAWKES_BENCH_AFFINITY=4 pins Windows workers to logical CPU 2.
args <- commandArgs(TRUE)
mode <- args[1]; lib <- args[2]; out <- args[3]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(emhawkes))
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
launch <- function(arguments) {
  status <- system2(rscript, c("--vanilla", vapply(arguments, shQuote, "")))
  stopifnot(status == 0L)
}
if (mode == "prepare") {
  stopifnot(!file.exists(file.path(out, "fixtures.rds")))
  fixtures <- list()
  # Independent event sequences exercise likelihood/inference scaling without
  # including simulation time. These are not used for estimation comparisons.
  for (n in c(1000L, 10000L, 100000L)) for (d in c(1L, 4L, 8L)) for (multiple in 1:2) {
    k <- d * multiple
    set.seed(61000L + d)
    dt <- c(0, rexp(n, rate = d * .4))
    types <- c(0L, sample.int(d, n, replace = TRUE))
    beta <- matrix(rep(seq(1.1, 2.1, length.out = k), each = d), d)
    h <- new("hspec", mu = rep(.3, d), alpha = matrix(.2 / k, d, k),
      beta = beta, type_col_map = lapply(seq_len(d), function(j) seq(j, k, by = d)))
    id <- sprintf("numeric-n%d-d%d-k%d", n, d, k)
    fixtures[[id]] <- list(kind = "numeric", events = n, dimensions = d, kernels = k,
      arguments = list(object = h, inter_arrival = dt, type = types,
                       lambda_component0 = matrix(.05 / k, d, k)))
  }
  # Keep few shared parameters to measure data scaling, rather than a growing
  # optimizer parameter count. These paths are simulated from their truth model.
  for (multiple in 1:2) {
    d <- 4L; k <- d * multiple
    h <- new("hspec", mu = rep(.3, d), alpha = matrix(.3 / k, d, k),
      beta = matrix(1.2, d, k),
      type_col_map = lapply(seq_len(d), function(j) seq(j, k, by = d)))
    set.seed(62000L + multiple)
    sim <- hsim(h, size = 10001L, lambda_component0 = .05)
    start <- h; start@mu[] <- .25; start@alpha[] <- .24 / k; start@beta[] <- 1
    for (n in c(1000L, 10000L)) {
      id <- sprintf("fit-n%d-d%d-k%d", n, d, k)
      fixtures[[id]] <- list(kind = "fit", events = n, dimensions = d, kernels = k,
        arguments = list(object = start, inter_arrival = sim$inter_arrival[seq_len(n + 1L)],
          type = sim$type[seq_len(n + 1L)], lambda_component0 = .05,
          method = "BFGS", constraints = list(ineqA = diag(3), ineqB = rep(-1e-5, 3)),
          control = list(iterlim = 500, reltol = 1e-8)))
    }
  }
  destination <- file.path(out, "bid-source")
  launch(c("tools/benchmark-bid-ask-numeric.R", "prepare", lib, destination))
  bid <- readRDS(file.path(destination, "fixtures.rds"))
  h <- bid$bid_ask_current$arguments$object
  set.seed(63001)
  sim <- hsim(h, size = 10001L, lambda_component0 = matrix(1, 4, 4))
  for (n in c(1000L, 5000L, 10000L)) for (variant in c(names(bid), "bid_ask_retained")) {
    a <- bid[[if (variant == "bid_ask_retained") "bid_ask" else variant]]$arguments
    if (variant == "bid_ask_retained") {
      fun <- a$object@impact
      capture <- new.env(parent = environment(fun))
      capture$.snapshots <- new.env(parent = emptyenv())
      environment(fun) <- capture
      body(fun) <- substitute({
        if (n %% 100L == 0L) assign(as.character(n), lambda_component, envir = .snapshots)
        ORIGINAL
      }, list(ORIGINAL = body(fun)))
      a$object@impact <- fun
    }
    a <- a[intersect(names(a), c("object", "lambda_component0"))]
    a$inter_arrival <- sim$inter_arrival[seq_len(n + 1L)]
    a$type <- sim$type[seq_len(n + 1L)]
    a$mark <- sim$mark[seq_len(n + 1L)]
    id <- paste0(variant, "-n", n)
    fixtures[[id]] <- list(kind = "bid", events = n, dimensions = 4L, kernels = 4L,
      arguments = a)
  }
  saveRDS(fixtures, file.path(out, "fixtures.rds"))
  writeLines(c(paste("source_commit", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("R", getRversion()), paste("package", packageVersion("emhawkes")),
    paste("maxLik", packageVersion("maxLik")), paste("Rcpp", packageVersion("Rcpp"))),
    file.path(out, "manifest.txt"))
  writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "session-info.txt"))
  quit()
}
fixtures <- readRDS(file.path(out, "fixtures.rds"))
if (mode == "run") {
  rounds <- if (length(args) >= 4L) as.integer(args[4]) else 3L
  stopifnot(rounds == 3L)
  for (id in names(fixtures)) {
    tasks <- if (fixtures[[id]]$kind == "numeric") c("likelihood", "inference") else
      if (fixtures[[id]]$kind == "fit") "fit" else c("prepare", "objective")
    for (task in tasks) for (round in c(seq_len(rounds), 0L)) {
      engines <- if (round %% 2L) c("r", "auto") else c("auto", "r")
      for (engine in engines) {
        file <- file.path(out, paste(id, task, engine, round, "rds", sep = "."))
        if (file.exists(file)) next
        cat(id, task, engine, round, "\n"); flush.console()
        launch(c("tools/benchmark-scaling.R", "worker", lib, out, id, task, engine, round))
      }
    }
  }
  quit()
}
if (mode == "worker") {
  id <- args[4]; task <- args[5]; engine <- args[6]; round <- as.integer(args[7])
  affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
  if (affinity != "scheduler") {
    stopifnot(.Platform$OS.type == "windows", grepl("^[1-9][0-9]*$", affinity))
    actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command",
      shQuote(sprintf("$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), affinity))), stdout = TRUE)
    stopifnot(length(actual) == 1L, trimws(actual) == affinity)
  }
  fixture <- fixtures[[id]]; a <- fixture$arguments
  options(emhawkes.engine = engine)
  observation_args <- a[intersect(names(a), c("object", "inter_arrival", "type", "mark"))]
  obs <- do.call(emhawkes:::normalize_observation, observation_args)
  map <- emhawkes:::build_parameter_map(a$object, TRUE)
  make <- function() {
    if (engine == "auto") emhawkes:::make_hawkes_objective(a$object, obs,
      a$lambda_component0, map, derivatives = FALSE) else {
      callbacks <- emhawkes:::model_callback_plans(a$object)
      list(loglik = function(param) {
        h <- a$object
        values <- emhawkes:::restore_parameters(map, param)
        for (name in names(values)) slot(h, name) <- values[[name]]
        emhawkes:::loglik_hspec_r(h, obs, a$lambda_component0,
          callbacks = callbacks, show_initial_message = FALSE)
      })
    }
  }
  if (task == "objective") objective <- make()
  operation <- switch(task,
    likelihood = function() do.call(logLik, a),
    inference = function() do.call(logLik, c(a, list(infer = TRUE))),
    fit = function() do.call(hfit, a),
    prepare = function() make(),
    objective = function() objective$loglik(map$start))
  # Warm the loaded code on a small prefix, without warming the measured
  # objective's one-point cache. Bid-ask preparation includes callback validation.
  warm <- a
  for (name in intersect(names(warm), c("inter_arrival", "type", "mark")))
    warm[[name]] <- warm[[name]][seq_len(min(51L, length(warm[[name]])))]
  if (task == "fit") { warm$control <- list(iterlim = 2L); suppressWarnings(suppressMessages(do.call(hfit, warm))) } else
    invisible(suppressMessages(do.call(logLik, warm)))
  warnings <- character()
  run <- function() withCallingHandlers(suppressMessages(operation()), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  gc()
  if (round == 0L) {
    baseline <- gc(reset = TRUE)
    profile <- tempfile()
    Rprofmem(profile)
    value <- run()
    Rprofmem(NULL)
    after <- gc()
    # Query before reading the allocation log so log parsing is excluded.
    peak_rss <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
      c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))), stdout = TRUE)) else NA_real_
    con <- file(profile, "r"); allocated <- allocations <- largest <- 0
    repeat {
      lines <- readLines(con, n = 10000L, warn = FALSE)
      if (!length(lines)) break
      bytes <- suppressWarnings(as.numeric(sub(" .*", "", substr(lines, 1L, 24L))))
      bytes <- bytes[is.finite(bytes)]
      allocated <- allocated + sum(bytes); allocations <- allocations + length(bytes)
      if (length(bytes)) largest <- max(largest, bytes)
    }
    close(con); unlink(profile)
    memory <- c(allocated_bytes = allocated, allocation_count = allocations,
      largest_block_bytes = largest, peak_r_heap_mb = sum(after[, 6]),
      incremental_r_heap_mb = sum(after[, 6]) - sum(baseline[, 2]), peak_rss_bytes = peak_rss)
    elapsed <- cpu <- NA_real_
  } else {
    batch <- if (fixture$kind == "numeric") {
      if (fixture$events == 1000L) { if (task == "likelihood") 200L else 50L } else
        if (fixture$events == 10000L) { if (task == "likelihood") 50L else 10L } else
          if (task == "likelihood") 5L else 2L
    } else 1L
    time <- system.time(for (iteration in seq_len(batch)) value <- run())
    elapsed <- unname(time[["elapsed"]]) / batch
    cpu <- unname(time[["user.self"]] + time[["sys.self"]]) / batch
    memory <- NULL
  }
  result <- switch(task,
    fit = value[c("estimate", "maximum", "gradient", "hessian", "code", "iterations")],
    prepare = NULL, value)
  computation <- if (task == "fit") fit_diagnostics(value)$computation else NULL
  saveRDS(list(id = id, task = task, engine = engine, round = round, seconds = elapsed,
    cpu_seconds = cpu, batch = if (round == 0L) 1L else batch, memory = memory, result = result, computation = computation,
    warnings = unique(warnings), affinity = affinity),
    file.path(out, paste(id, task, engine, round, "rds", sep = ".")))
  quit()
}
if (mode == "snapshots") {
  rows <- list()
  for (n in c(1000L, 5000L, 10000L)) {
    snapshots <- likelihoods <- list()
    for (engine in c("r", "auto")) {
      a <- readRDS(file.path(out, "fixtures.rds"))[[paste0("bid_ask_retained-n", n)]]$arguments
      options(emhawkes.engine = engine)
      if (engine == "r") likelihoods[[engine]] <- as.numeric(suppressMessages(do.call(logLik, a))) else {
        obs <- do.call(emhawkes:::normalize_observation, a[names(a) != "lambda_component0"])
        map <- emhawkes:::build_parameter_map(a$object, TRUE)
        objective <- emhawkes:::make_hawkes_objective(a$object, obs, a$lambda_component0,
          map, derivatives = FALSE)
        likelihoods[[engine]] <- as.numeric(objective$loglik(map$start))
      }
      kept <- as.list(environment(a$object@impact)$.snapshots)
      kept <- kept[order(as.integer(names(kept)))]
      stopifnot(length(kept) == n / 100L)
      for (name in names(kept)) {
        state <- kept[[name]]
        stopifnot(all(state[seq.int(as.integer(name) + 1L, nrow(state)), , drop = FALSE] == 0))
      }
      snapshots[[engine]] <- kept
    }
    stopifnot(isTRUE(all.equal(snapshots$r, snapshots$auto, tolerance = 1e-12)),
      isTRUE(all.equal(likelihoods$r, likelihoods$auto, tolerance = 1e-12)))
    rows[[length(rows) + 1L]] <- data.frame(events = n, snapshots = length(snapshots$r),
      matrices_exact = identical(snapshots$r, snapshots$auto), tolerance_1e_12 = TRUE,
      future_rows_zero = TRUE, loglik_exact = identical(likelihoods$r, likelihoods$auto),
      loglik_max_abs_difference = max(abs(likelihoods$r - likelihoods$auto)))
  }
  write.csv(do.call(rbind, rows), file.path(out, "bid-snapshots.csv"), row.names = FALSE)
  quit()
}
stopifnot(mode == "compare")
files <- list.files(out, pattern = "[.](r|auto)[.][0-9]+[.]rds$", full.names = TRUE)
rows <- checks <- memory <- cross_checks <- fit_fields <- estimates <- list()
for (id in names(fixtures)) {
  cat("Checking", id, "\n"); flush.console()
  expected <- if (fixtures[[id]]$kind == "numeric") c("likelihood", "inference") else
    if (fixtures[[id]]$kind == "fit") "fit" else c("prepare", "objective")
  for (task in expected) {
    paths <- file.path(out, as.vector(outer(c("r", "auto"), 0:3,
      function(engine, round) paste(id, task, engine, round, "rds", sep = "."))))
    stopifnot(all(file.exists(paths)))
    ref <- readRDS(file.path(out, paste(id, task, "r", 1L, "rds", sep = ".")))$result
    for (path in paths) {
      record <- readRDS(path)
      exact <- identical(record$result, ref)
      equal <- exact || isTRUE(all.equal(record$result, ref, tolerance = 1e-9))
      checks[[length(checks) + 1L]] <- data.frame(id, task, engine = record$engine,
        round = record$round, exact, tolerance_1e_9 = equal,
        max_abs_difference = if (length(ref)) max(abs(unlist(record$result, use.names = FALSE) - unlist(ref, use.names = FALSE))) else 0,
        warnings = paste(record$warnings, collapse = " | "))
      stopifnot(equal)
      if (task == "fit") stopifnot(record$result$code == 0L,
        identical(record$computation$engine$used, if (record$engine == "auto") "cpp" else "r"))
      if (task == "fit") {
        for (field in names(ref)) fit_fields[[length(fit_fields) + 1L]] <- data.frame(
          id, engine = record$engine, round = record$round, field,
          exact = identical(record$result[[field]], ref[[field]]),
          max_abs_difference = max(abs(record$result[[field]] - ref[[field]])))
        estimates[[length(estimates) + 1L]] <- data.frame(id, engine = record$engine,
          round = record$round, parameter = names(record$result$estimate),
          estimate = unname(record$result$estimate))
      }
      info <- fixtures[[id]]
      if (record$round == 0L) memory[[length(memory) + 1L]] <- cbind(
        data.frame(id, task, engine = record$engine), as.data.frame(as.list(record$memory))) else
        rows[[length(rows) + 1L]] <- data.frame(id, task, engine = record$engine,
          round = record$round, batch = record$batch, events = info$events, dimensions = info$dimensions,
          kernels = info$kernels, seconds = record$seconds, cpu_seconds = record$cpu_seconds,
          affinity = record$affinity)
    }
  }
}
for (n in c(1000L, 5000L, 10000L)) for (engine in c("r", "auto")) {
  get_value <- function(id) readRDS(file.path(out, paste(id, "objective", engine, 1L, "rds", sep = ".")))$result
  current <- get_value(paste0("bid_ask_current-n", n))
  history <- get_value(paste0("bid_ask-n", n))
  retained <- get_value(paste0("bid_ask_retained-n", n))
  stopifnot(identical(current, history), identical(current, retained))
  cross_checks[[length(cross_checks) + 1L]] <- data.frame(events = n, engine,
    current_history_exact = TRUE, retained_exact = TRUE, max_abs_difference = max(abs(unlist(current) - unlist(history))))
}
write.csv(do.call(rbind, cross_checks), file.path(out, "bid-variants.csv"), row.names = FALSE)
timings <- do.call(rbind, rows)
summary <- aggregate(cbind(seconds, cpu_seconds) ~ id + task + engine + events + dimensions + kernels,
  timings, median)
for (name in c("timings", "summary", "checks", "memory", "fit_fields", "estimates")) {
  value <- get(name)
  if (is.list(value) && !is.data.frame(value)) value <- do.call(rbind, value)
  write.csv(value, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}
print(summary, row.names = FALSE)
