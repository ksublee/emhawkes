# Fresh R process per package; run from the repository root.
# Rscript --vanilla tools/benchmark-cran-r-engine.R install|prepare|run|worker|compare [cran|current] [case]
args <- commandArgs(TRUE)
mode <- args[1]
out <- "tools/cran-r-engine-results"
dir.create(out, showWarnings = FALSE)
state_file <- file.path(out, "installation.rds")
if (mode == "install") {
  root <- tempfile("emhawkes-cran-r-engine-", tmpdir = dirname(tempdir()))
  dir.create(root)
  stage <- file.path(root, "current-source")
  dir.create(stage)
  files <- c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
    list.files("src", pattern = "[.](cpp|h)$|^Makevars", full.names = TRUE))
  for (file in files) {
    target <- file.path(stage, file)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file, target))
  }
  libs <- setNames(file.path(root, c("cran-library", "current-library")), c("cran", "current"))
  invisible(lapply(libs, dir.create))
  archive <- file.path(root, "emhawkes_0.9.8.tar.gz")
  download.file("https://cran.r-project.org/src/contrib/emhawkes_0.9.8.tar.gz", archive, mode = "wb")
  for (version in names(libs)) {
    package <- if (version == "cran") archive else stage
    status <- system2(file.path(R.home("bin"), "R.exe"), c("CMD", "INSTALL", "--no-multiarch",
      "--no-docs", "--no-help", paste0("--library=", shQuote(libs[[version]])), shQuote(package)))
    stopifnot(status == 0L)
  }
  saveRDS(list(libs = libs, dependencies = .libPaths()), state_file)
  write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
    file.path(out, "source-manifest.csv"), row.names = FALSE)
  quit()
}
state <- readRDS(state_file)
version <- if (mode == "worker") args[2] else "current"
stopifnot(dir.exists(state$libs[[version]]))
.libPaths(c(state$libs[[version]], state$dependencies, .libPaths()))
suppressPackageStartupMessages(library(emhawkes))
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))),
                    normalizePath(state$libs[[version]])))
options(emhawkes.engine = "auto")
if (mode == "prepare") {
  source("tools/observed-impact-cases.R")
  impact <- function(param = c(extra = .15), n, N, ...) {
    matrix(param["extra"] * (N[n, 1] %% 2))
  }
  truth <- new("hspec", mu = .3, alpha = .2, beta = 1.2, impact = impact)
  start <- new("hspec", mu = .25, alpha = .15, beta = 1, impact = impact)
  formals(start@impact)$param <- c(extra = .1)
  cases <- c(list(user_example = list(truth = truth, start = start, seed = 9351L,
    initial = .1)), observed_impact_cases()[c("bivariate_counts_marks", "multi_kernel_counts")])
  fixtures <- list()
  for (name in names(cases)) {
    spec <- cases[[name]]
    set.seed(spec$seed)
    counts <- if (name == "user_example") list() else list(N0 = 2, Nc0 = 3)
    data <- do.call(hsim, c(list(object = spec$truth, size = 2001L,
      lambda_component0 = .1), counts))
    p <- length(emhawkes:::build_parameter_map(spec$start, TRUE)$start)
    a <- list(object = spec$start, inter_arrival = data$inter_arrival,
      type = data$type, mark = data$mark, lambda_component0 = spec$initial,
      method = "BFGS", constraints = list(ineqA = diag(p), ineqB = rep(-1e-5, p)),
      control = list(iterlim = 500, reltol = 1e-8))
    if (length(counts)) {
      # CRAN hfit lacks Nc0; supply the same observed marked count path to both.
      a$N0 <- 2
      a$Nc <- data$Nc
    }
    fixtures[[name]] <- list(arguments = a, events = 2000L)
  }
  # Extract the two extended examples while intercepting their hfit calls.
  lines <- readLines("vignettes_ext/extended_example.Rmd", warn = FALSE)
  contains <- function(x, target) {
    if (!is.call(x)) return(FALSE)
    if (identical(x[[1]], as.name(target))) return(TRUE)
    any(vapply(as.list(x)[-1], contains, logical(1), target = target))
  }
  config <- list(flocking = c("### Hawkes flocking model", "fit_fl"),
                 bid_ask = c("### Bid-ask price model", "mle_ba"))
  for (name in names(config)) {
    start <- match(config[[name]][1], lines)
    headings <- grep("^#{2,3} ", lines)
    end <- min(c(headings[headings > start], length(lines) + 1L))
    env <- new.env(parent = globalenv())
    env$hfit <- function(...) list(arguments = list(...))
    for (opening in grep("^```\\{r", lines)) {
      if (opening <= start || opening >= end) next
      closing <- opening + which(lines[(opening + 1L):length(lines)] == "```")[1]
      for (expr in parse(text = lines[(opening + 1L):(closing - 1L)])) {
        if (!is.call(expr)) next
        assignment <- identical(expr[[1]], as.name("<-"))
        if (!assignment && !identical(expr[[1]], as.name("set.seed"))) next
        if (assignment && (contains(expr, "logLik") || contains(expr, "infer_lambda"))) next
        eval(expr, env)
        if (assignment && identical(expr[[2]], as.name(config[[name]][2])) && contains(expr, "hfit")) {
          a <- env[[config[[name]][2]]]$arguments
          if (is.null(names(a))) names(a) <- rep("", length(a))
          positional <- c("object", "inter_arrival", "type", "mark")
          for (i in seq_along(a)) if (names(a)[i] == "") names(a)[i] <- positional[i]
          a$final_hessian <- NULL
          fixtures[[name]] <- list(arguments = a, events = length(a$inter_arrival) - 1L)
          break
        }
      }
      if (name %in% names(fixtures)) break
    }
    stopifnot(name %in% names(fixtures))
  }
  saveRDS(fixtures, file.path(out, "fixtures.rds"))
  quit()
}
if (mode == "run") {
  cases <- names(readRDS(file.path(out, "fixtures.rds")))
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  for (name in cases) for (version in c("cran", "current")) {
    status <- system2(executable, c("--vanilla", "tools/benchmark-cran-r-engine.R", "worker", version, name))
    stopifnot(status == 0L)
  }
  quit()
}
if (mode == "worker") {
  fixture <- readRDS(file.path(out, "fixtures.rds"))[[args[3]]]
  a <- fixture$arguments
  if (version == "cran") {
    a$constraint <- a$constraints
    a$constraints <- NULL
  }
  ll_args <- a[intersect(names(a), c("object", "inter_arrival", "type", "mark",
    "N0", "Nc", "lambda_component0"))]
  initial_loglik <- tryCatch(suppressMessages(do.call(logLik, ll_args)), error = identity)
  if (inherits(initial_loglik, "error")) {
    message <- conditionMessage(initial_loglik)
    saveRDS(list(case = args[3], version = version, times = numeric(0),
      error = message, package_version = as.character(packageVersion("emhawkes")),
      maxLik_version = as.character(packageVersion("maxLik"))),
      file.path(out, paste0(version, "-", args[3], ".rds")))
    cat(version, args[3], "cannot evaluate the specified likelihood:", message, "\n")
    quit()
  }
  fit_once <- function() suppressWarnings(suppressMessages(do.call(hfit, a)))
  invisible(fit_once())
  times <- numeric(3)
  results <- list()
  for (round in 1:3) {
    gc()
    times[round] <- unname(system.time(fit <- fit_once())[["elapsed"]])
    results[[round]] <- fit[c("estimate", "maximum", "gradient", "hessian", "code", "iterations")]
    cat(version, args[3], round, times[round], "seconds; code", fit$code, "\n")
    flush.console()
    saveRDS(list(case = args[3], version = version, times = times[seq_len(round)],
      results = results, initial_loglik = initial_loglik,
      package_version = as.character(packageVersion("emhawkes")),
      maxLik_version = as.character(packageVersion("maxLik"))),
      file.path(out, paste0(version, "-", args[3], ".rds")))
  }
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
    file.path(out, paste0("session-", version, ".txt")))
  quit()
}
if (mode == "compare") {
  cases <- names(readRDS(file.path(out, "fixtures.rds")))
  rows <- timings <- estimates <- list()
  for (name in cases) {
    a <- readRDS(file.path(out, paste0("cran-", name, ".rds")))
    b <- readRDS(file.path(out, paste0("current-", name, ".rds")))
    stopifnot(a$maxLik_version == b$maxLik_version)
    for (value in list(a, b)) if (length(value$times)) {
      timings[[length(timings) + 1L]] <- data.frame(case = name, version = value$version,
        round = seq_along(value$times), elapsed_seconds = value$times)
    }
    events <- readRDS(file.path(out, "fixtures.rds"))[[name]]$events
    if (!is.null(a$error)) {
      rows[[name]] <- data.frame(case = name, events = events,
        cran_seconds = NA_real_, current_seconds = median(b$times),
        cran_loglik_length = NA_integer_, current_loglik_length = length(b$initial_loglik),
        cran_initial_loglik_sum = NA_real_, current_initial_loglik = sum(b$initial_loglik),
        estimate_max_abs_difference = NA_real_, maximum_abs_difference = NA_real_,
        cran_code = NA_integer_, current_code = tail(b$results, 1)[[1]]$code,
        cran_iterations = NA_integer_, current_iterations = tail(b$results, 1)[[1]]$iterations,
        cran_error = a$error)
      next
    }
    x <- tail(a$results, 1)[[1]]; y <- tail(b$results, 1)[[1]]
    stopifnot(identical(sort(names(x$estimate)), sort(names(y$estimate))), x$code == 0L, y$code == 0L)
    estimates[[name]] <- data.frame(case = name, parameter = names(x$estimate),
      cran = unname(x$estimate), current = unname(y$estimate[names(x$estimate)]))
    rows[[name]] <- data.frame(case = name, events = events,
      cran_seconds = median(a$times), current_seconds = median(b$times),
      cran_loglik_length = length(a$initial_loglik), current_loglik_length = length(b$initial_loglik),
      cran_initial_loglik_sum = sum(a$initial_loglik), current_initial_loglik = sum(b$initial_loglik),
      estimate_max_abs_difference = max(abs(x$estimate - y$estimate[names(x$estimate)])),
      maximum_abs_difference = abs(x$maximum - y$maximum), cran_code = x$code, current_code = y$code,
      cran_iterations = x$iterations, current_iterations = y$iterations, cran_error = "")
  }
  write.csv(do.call(rbind, rows), file.path(out, "comparison.csv"), row.names = FALSE)
  write.csv(do.call(rbind, timings), file.path(out, "timings.csv"), row.names = FALSE)
  write.csv(do.call(rbind, estimates), file.path(out, "estimates.csv"), row.names = FALSE)
  print(do.call(rbind, rows), row.names = FALSE)
}
