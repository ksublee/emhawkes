library(emhawkes)

# Copy only the package and its compiled dependency into an isolated library.
# The new R process must not find maxLik anywhere in its library search path.
isolated_stats_session <- function() {
  root <- tempfile("stats-without-maxLik-")
  dir.create(root)
  isolated <- file.path(root, "library")
  dir.create(isolated)
  for (package in c("emhawkes", "Rcpp")) {
    stopifnot(file.copy(find.package(package), isolated, recursive = TRUE))
  }
  script <- file.path(root, "check.R")
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(c(args[1], .Library))",
    "stopifnot(!requireNamespace('maxLik', quietly = TRUE))",
    "library(emhawkes)",
    "stopifnot(!'maxLik' %in% loadedNamespaces())",
    "h <- new('hspec', mu = .5, alpha = 0, beta = 1)",
    "dt <- c(0,.2,.3,.5)",
    "fit <- hfit(h, dt, optimizer = 'stats', lambda_component0 = 0,",
    "  fixed = c(FALSE,TRUE,TRUE), gradient = 'analytic',",
    "  control = list(reltol = 1e-12))",
    "stopifnot(abs(coef(fit)[1]-3) < 1e-6, abs(vcov(fit)[1,1]-3) < 1e-6)",
    "saveRDS(fit, args[2])",
    "fit <- readRDS(args[2])",
    "stopifnot(inherits(summary(fit), 'summary.emhawkes_stats_fit'))",
    "invisible(capture.output(print(fit)))",
    "stopifnot(fit_diagnostics(fit)$status == 'ok')",
    "restored <- fitted_hspec(fit)",
    "replayed <- do.call(logLik, c(list(object = restored, inter_arrival = dt), fit_initial_state(fit)))",
    "stopifnot(abs(as.numeric(replayed)-as.numeric(logLik(fit))) < 1e-12)",
    "error <- tryCatch(hfit(h, dt), error = identity)",
    "stopifnot(inherits(error,'error'), grepl(\"optimizer = 'stats'\", conditionMessage(error), fixed=TRUE))",
    "stopifnot(!'maxLik' %in% loadedNamespaces())",
    "cat('Isolated stats fit and saved-result round trip: OK\\n')"
  ), script)
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  output <- system2(executable, c("--vanilla", shQuote(script), shQuote(isolated),
    shQuote(file.path(root, "fit.rds"))), stdout = TRUE, stderr = TRUE)
  cat(output, sep = "\n")
  stopifnot(is.null(attr(output, "status")))
  # Every removed file is contained in this test's own temporary directory.
  stopifnot(startsWith(normalizePath(root, winslash = "/"),
                       paste0(normalizePath(tempdir(), winslash = "/"), "/")))
  unlink(root, recursive = TRUE)
}
isolated_stats_session()
