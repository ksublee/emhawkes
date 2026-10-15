# Complete hfit() timings in fresh R processes; numerical differentiation only.
# Prepare identical seeded fixtures with the pre-change package:
# Rscript tools/benchmark-impact-numeric.R prepare OLD_LIB OUT
# Optional fourth argument OUT copies fixtures from an existing third-argument RDS.
# Rscript tools/benchmark-impact-numeric.R worker LIB ENGINE CASE REPS OUT
# Rscript tools/benchmark-impact-numeric.R compare NEW_LIB OUT
# Rscript tools/benchmark-impact-numeric.R points LIB ENGINE OUT
# Run points with cpp first, then r and cran, for common-parameter checks.
# ENGINE: cran, old, r or cpp. Use the same maxLik installation in every worker.
args <- commandArgs(TRUE)
mode <- args[1]
lib <- args[2]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(emhawkes))
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
if (mode == "prepare") {
  out <- tail(args, 1)
  if (length(args) == 4L) {
    fixtures <- readRDS(args[3])[c("user_example", "bivariate_counts_marks", "flocking")]
  } else {
    source("tools/observed-impact-cases.R")
    impact <- function(param = c(extra = .15), n, N, ...) matrix(param["extra"] * (N[n, 1] %% 2))
    truth <- new("hspec", mu = .3, alpha = .2, beta = 1.2, impact = impact)
    start <- new("hspec", mu = .25, alpha = .15, beta = 1, impact = impact)
    formals(start@impact)$param <- c(extra = .1)
    cases <- c(list(user_example = list(truth = truth, start = start, seed = 9351L, initial = .1)),
      observed_impact_cases()["bivariate_counts_marks"])
    fixtures <- list()
    for (name in names(cases)) {
      spec <- cases[[name]]
      set.seed(spec$seed)
      counts <- if (name == "user_example") list() else list(N0 = 2, Nc0 = 3)
      data <- do.call(hsim, c(list(object = spec$truth, size = 2001L, lambda_component0 = .1), counts))
      p <- length(emhawkes:::build_parameter_map(spec$start, TRUE)$start)
      a <- list(object = spec$start, inter_arrival = data$inter_arrival, type = data$type,
        mark = data$mark, lambda_component0 = spec$initial, method = "BFGS",
        constraints = list(ineqA = diag(p), ineqB = rep(-1e-5, p)),
        control = list(iterlim = 500, reltol = 1e-8))
      if (length(counts)) { a$N0 <- 2; a$Nc <- data$Nc }
      fixtures[[name]] <- list(arguments = a, events = 2000L)
    }
    # Evaluate vignette assignments up to its first fit, intercepting that fit.
    lines <- readLines("vignettes_ext/extended_example.Rmd", warn = FALSE)
    begin <- match("### Hawkes flocking model", lines)
    headings <- grep("^#{2,3} ", lines)
    end <- min(c(headings[headings > begin], length(lines) + 1L))
    env <- new.env(parent = globalenv())
    env$hfit <- function(...) list(arguments = list(...))
    for (opening in grep("^```\\{r", lines)) {
      if (opening <= begin || opening >= end) next
      closing <- opening + which(lines[(opening + 1L):length(lines)] == "```")[1]
      for (expr in parse(text = lines[(opening + 1L):(closing - 1L)])) {
        if (!is.call(expr)) next
        assignment <- identical(expr[[1]], as.name("<-"))
        if (!assignment && !identical(expr[[1]], as.name("set.seed"))) next
        if (assignment && identical(expr[[2]], as.name("lambda_fl"))) next
        eval(expr, env)
        if (assignment && identical(expr[[2]], as.name("fit_fl"))) break
      }
      if (exists("fit_fl", envir = env, inherits = FALSE)) break
    }
    a <- env$fit_fl$arguments
    if (is.null(names(a))) names(a) <- rep("", length(a))
    for (i in seq_along(a)) if (names(a)[i] == "") names(a)[i] <- c("object", "inter_arrival", "type", "mark")[i]
    a$final_hessian <- NULL
    fixtures$flocking <- list(arguments = a, events = length(a$inter_arrival) - 1L)
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fixtures, file.path(out, "fixtures.rds"))
  quit()
}
if (mode == "points") {
  engine <- args[3]
  out <- args[4]
  fixtures <- readRDS(file.path(out, "fixtures.rds"))
  if (engine == "cpp") {
    points <- lapply(names(fixtures), function(case) {
      a <- fixtures[[case]]$arguments
      map <- emhawkes:::build_parameter_map(a$object, TRUE)
      fit <- readRDS(file.path(out, paste0("cpp-", case, ".rds")))$results[[1]]
      lapply(list(start = map$start, fitted = fit$estimate, nearby = fit$estimate * 1.001), function(param) {
        object <- a$object
        values <- emhawkes:::restore_parameters(map, param)
        for (name in names(values)) slot(object, name) <- values[[name]]
        object
      })
    })
    names(points) <- names(fixtures)
    saveRDS(points, file.path(out, "point-models.rds"))
  } else points <- readRDS(file.path(out, "point-models.rds"))
  rows <- list()
  for (case in names(fixtures)) for (point in names(points[[case]])) {
    a <- fixtures[[case]]$arguments
    a$object <- points[[case]][[point]]
    a <- a[intersect(names(a), c("object", "inter_arrival", "type", "mark", "N", "Nc", "N0", "Nc0", "lambda_component0"))]
    if (engine == "cpp") {
      map <- emhawkes:::build_parameter_map(a$object, TRUE)
      obs <- do.call(emhawkes:::normalize_observation, a[names(a) != "lambda_component0"])
      objective <- emhawkes:::make_hawkes_objective(a$object, obs, a$lambda_component0,
        map, derivatives = FALSE)
      value <- objective$loglik(map$start)
    } else {
      options(emhawkes.engine = "r")
      value <- suppressMessages(sum(do.call(logLik, a)))
    }
    rows[[paste(case, point)]] <- data.frame(case = case, point = point, engine = engine, loglik = value)
  }
  write.csv(do.call(rbind, rows), file.path(out, paste0("fixed-points-", engine, ".csv")), row.names = FALSE)
  print(do.call(rbind, rows), row.names = FALSE)
  quit()
}
if (mode == "compare") {
  out <- args[3]
  rows <- checks <- estimates <- timings <- list()
  for (case in c("user_example", "bivariate_counts_marks", "flocking")) {
    versions <- setNames(lapply(c("old", "cran", "r", "cpp"), function(engine) {
      readRDS(file.path(out, paste0(engine, "-", case, ".rds")))
    }), c("old", "cran", "r", "cpp"))
    stopifnot(length(unique(vapply(versions, `[[`, "", "maxLik_version"))) == 1L)
    ref <- versions$r$results[[1]]
    compiled <- versions$cpp$results[[1]]
    cran <- versions$cran$results[[1]]
    stopifnot(identical(ref, versions$old$results[[1]]), identical(ref, compiled))
    for (engine in names(versions)) {
      value <- versions[[engine]]
      stopifnot(all(vapply(value$results, function(x) identical(x, value$results[[1]]), logical(1))))
      timings[[paste(case, engine)]] <- data.frame(case = case, engine = engine,
        round = seq_along(value$times), elapsed_seconds = value$times)
    }
    rows[[case]] <- data.frame(case = case, events = versions$cpp$events,
      cran_seconds = median(versions$cran$times), r_seconds = median(versions$r$times),
      cpp_seconds = median(versions$cpp$times),
      cran_speedup = median(versions$cran$times) / median(versions$cpp$times))
    checks[[case]] <- data.frame(case = case, previous_default_identical = TRUE,
      cpp_six_fields_identical = TRUE,
      cran_estimate_max_abs_difference = max(abs(cran$estimate - compiled$estimate)),
      cran_maximum_abs_difference = abs(cran$maximum - compiled$maximum),
      cran_iterations = cran$iterations, cpp_iterations = compiled$iterations,
      cran_code = cran$code, cpp_code = compiled$code)
    estimates[[case]] <- data.frame(case = case, parameter = names(ref$estimate),
      cran = unname(cran$estimate), r = unname(ref$estimate), cpp = unname(compiled$estimate))
  }
  for (name in c("rows", "checks", "estimates", "timings")) {
    file <- switch(name, rows = "comparison", checks = "accuracy", name)
    write.csv(do.call(rbind, get(name)), file.path(out, paste0(file, ".csv")), row.names = FALSE)
  }
  print(do.call(rbind, rows), row.names = FALSE)
  print(do.call(rbind, checks), row.names = FALSE)
  quit()
}
stopifnot(mode == "worker")
engine <- args[3]
case <- args[4]
reps <- as.integer(args[5])
out <- args[6]
fixture <- readRDS(file.path(out, "fixtures.rds"))[[case]]
a <- fixture$arguments
options(emhawkes.engine = if (engine %in% c("cran", "old")) "auto" else engine)
if (engine == "cran") {
  a$constraint <- a$constraints
  a$constraints <- NULL
}
fit_once <- function(arguments = a) suppressWarnings(suppressMessages(do.call(hfit, arguments)))
# Exercise the public setup and optimizer before timing, on a shorter path.
warm <- a
for (name in intersect(names(a), c("inter_arrival", "type", "mark"))) warm[[name]] <- a[[name]][seq_len(min(101, length(a$inter_arrival))) ]
for (name in intersect(names(a), c("N", "Nc"))) warm[[name]] <- a[[name]][seq_len(min(101, length(a$inter_arrival))), , drop = FALSE]
warm$control <- list(iterlim = 2)
invisible(fit_once(warm))
ll_args <- a[intersect(names(a), c("object", "inter_arrival", "type", "mark", "N", "Nc", "N0", "Nc0", "lambda_component0"))]
options(emhawkes.engine = "r")
fixed_loglik <- suppressMessages(as.numeric(do.call(logLik, ll_args)))
options(emhawkes.engine = if (engine %in% c("cran", "old")) "auto" else engine)
results <- list()
times <- numeric()
for (i in seq_len(reps)) {
  gc()
  times[i] <- unname(system.time(fit <- fit_once())[["elapsed"]])
  results[[i]] <- fit[c("estimate", "maximum", "gradient", "hessian", "code", "iterations")]
  result <- list(case = case, engine = engine, events = fixture$events, times = times,
    results = results, fixed_loglik = fixed_loglik,
    package_version = as.character(packageVersion("emhawkes")),
    maxLik_version = as.character(packageVersion("maxLik")),
    initial_state = if (engine %in% c("r", "cpp")) fit_initial_state(fit) else NULL)
  saveRDS(result, file.path(out, paste0(engine, "-", case, ".rds")))
  cat(engine, case, i, times[i], "seconds; code", fit$code, "\n")
  flush.console()
}
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(out, paste0("session-", engine, ".txt")))
