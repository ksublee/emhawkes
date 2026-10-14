# Run from the repository root with the candidate package installed.
# Rscript tools/benchmark-r-engine.R [baseline_ref] [events] [repetitions] [output]
# Only the R likelihood function is exchanged in this disposable R session.
benchmark_r_engine <- function() {
  library(emhawkes)
  source("tools/observed-impact-cases.R")
  args <- commandArgs(trailingOnly = TRUE)
  ref <- if (length(args)) args[1] else "63dd422"
  events <- if (length(args) >= 2L) as.integer(args[2]) else 2000L
  repetitions <- if (length(args) >= 3L) as.integer(args[3]) else 3L
  output <- if (length(args) >= 4L) args[4] else tempfile("r-engine-results-")
  stopifnot(grepl("^[a-zA-Z0-9._/-]+$", ref), events >= 32L, repetitions >= 1L)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)

  old_source <- system2("git", c("show", paste0(ref, ":R/hllf.R")), stdout = TRUE)
  stopifnot(is.null(attr(old_source, "status")))
  definitions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("loglik_hspec_r")), as.list(parse(text = old_source)))
  stopifnot(length(definitions) == 1L)
  old_env <- new.env(parent = asNamespace("emhawkes"))
  eval(definitions[[1]], old_env)
  environment(old_env$loglik_hspec_r) <- asNamespace("emhawkes")
  current <- getFromNamespace("loglik_hspec_r", "emhawkes")
  engines <- list(before = old_env$loglik_hspec_r, after = current)
  on.exit(assignInNamespace("loglik_hspec_r", current, "emhawkes"), add = TRUE)
  options_old <- options(emhawkes.engine = "auto")
  on.exit(options(options_old), add = TRUE)

  # Same specification and seed as the example supplied for the user's timing.
  impact <- function(param = c(extra = .15), n, N, ...) {
    matrix(param["extra"] * (N[n, 1] %% 2))
  }
  truth <- new("hspec", mu = .3, alpha = .2, beta = 1.2, impact = impact)
  start <- new("hspec", mu = .25, alpha = .15, beta = 1, impact = impact)
  formals(start@impact)$param <- c(extra = .1)
  cases <- c(list(user_example = list(truth = truth, start = start, seed = 9351L,
    initial = .1)), observed_impact_cases()[c("bivariate_counts_marks", "multi_kernel_counts")])
  timings <- comparisons <- estimates <- list()
  fields <- c("estimate", "maximum", "gradient", "hessian", "code", "iterations")
  for (name in names(cases)) {
    spec <- cases[[name]]
    set.seed(spec$seed)
    counts <- if (name == "user_example") list() else list(N0 = 2, Nc0 = 3)
    data <- do.call(hsim, c(list(object = spec$truth, size = events + 1L,
      lambda_component0 = .1), counts))
    p <- length(emhawkes:::build_parameter_map(spec$start, TRUE)$start)
    fit_args <- c(list(object = spec$start, inter_arrival = data$inter_arrival,
      type = data$type, mark = data$mark, lambda_component0 = spec$initial,
      method = "BFGS", constraints = list(ineqA = diag(p), ineqB = rep(-1e-5, p)),
      control = list(iterlim = 500, reltol = 1e-8)), counts)
    fit_once <- function() suppressMessages(do.call(hfit, fit_args))
    reference <- NULL
    for (engine in names(engines)) {
      assignInNamespace("loglik_hspec_r", engines[[engine]], "emhawkes")
      fit <- fit_once()
      if (is.null(reference)) reference <- fit[fields]
      stopifnot(identical(reference, fit[fields]))
    }
    for (round in seq_len(repetitions)) {
      order <- if (round %% 2L) names(engines) else rev(names(engines))
      for (engine in order) {
        assignInNamespace("loglik_hspec_r", engines[[engine]], "emhawkes")
        gc()
        elapsed <- unname(system.time(fit <- fit_once())[["elapsed"]])
        exact <- vapply(fields, function(field) identical(reference[[field]], fit[[field]]), logical(1))
        stopifnot(all(exact), fit$code == 0L)
        timings[[length(timings) + 1L]] <- data.frame(case = name, round = round,
          engine = engine, elapsed_seconds = elapsed)
        comparisons[[length(comparisons) + 1L]] <- data.frame(case = name,
          round = round, engine = engine, field = fields, identical = exact)
        cat(name, round, engine, elapsed, "seconds; all fit fields identical\n")
        flush.console()
        write.csv(do.call(rbind, timings), file.path(output, "timings.csv"), row.names = FALSE)
      }
    }
    estimates[[name]] <- data.frame(case = name, parameter = names(fit$estimate),
      estimate = unname(fit$estimate), standard_error = sqrt(diag(-solve(fit$hessian))))
    saveRDS(list(reference = reference, result = fit[fields]), file.path(output, paste0(name, ".rds")))
  }
  write.csv(do.call(rbind, comparisons), file.path(output, "comparisons.csv"), row.names = FALSE)
  write.csv(do.call(rbind, estimates), file.path(output, "estimates.csv"), row.names = FALSE)
  writeLines(c(paste("Baseline R likelihood:", ref), paste("Events:", events),
    paste("Repetitions:", repetitions), capture.output(sessionInfo())), file.path(output, "session-info.txt"))
  cat("Results:", normalizePath(output, winslash = "/"), "\n")
}
benchmark_r_engine()
