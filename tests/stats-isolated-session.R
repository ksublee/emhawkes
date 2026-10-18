library(emhawkes)

# Give the child its own emhawkes/Rcpp library, retaining R's standard/site
# libraries. maxLik may legitimately be installed there. Always test the stats
# path without loading maxLik; test its absence only when it is actually absent.
isolated_stats_session <- function() {
  root <- tempfile("stats-without-maxLik-")
  dir.create(root)
  stopifnot(startsWith(normalizePath(root, winslash = "/"),
                       paste0(normalizePath(tempdir(), winslash = "/"), "/")))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  isolated <- file.path(root, "library")
  dir.create(isolated)
  for (package in c("emhawkes", "Rcpp")) {
    stopifnot(file.copy(find.package(package), isolated, recursive = TRUE))
  }
  script <- file.path(root, "check.R")
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(c(args[1], .Library))",
    "stopifnot(!('package:methods' %in% search()))",
    "maxlik_available <- length(find.package('maxLik', quiet = TRUE)) > 0L",
    "if (identical(args[3], 'true')) stopifnot(!maxlik_available)",
    "library(emhawkes)",
    "stopifnot(!('package:methods' %in% search()))",
    "stopifnot(!'maxLik' %in% loadedNamespaces())",
    "h <- methods::new('hspec', mu = .5, alpha = 0, beta = 1)",
    "stopifnot(length(h@dimens) == 1L, h@dimens == 1, identical(dim(h@mu), c(1L, 1L)))",
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
    "stopifnot(length(restored@dimens) == 1L, restored@dimens == 1, identical(dim(restored@mu), c(1L, 1L)))",
    "replayed <- do.call(logLik, c(list(object = restored, inter_arrival = dt), fit_initial_state(fit)))",
    "stopifnot(abs(as.numeric(replayed)-as.numeric(logLik(fit))) < 1e-12)",
    "if (!maxlik_available) {",
    "  error <- tryCatch(hfit(h, dt), error = identity)",
    "  stopifnot(inherits(error,'error'), grepl(\"optimizer = 'stats'\", conditionMessage(error), fixed=TRUE))",
    "}",
    "stopifnot(!'maxLik' %in% loadedNamespaces())",
    "cat('Stats fit and saved-result round trip: OK; maxLik visible:', maxlik_available, '\\n')"
  ), script)
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  require_missing <- identical(Sys.getenv("EMHAWKES_TEST_REQUIRE_NO_MAXLIK"), "true")
  output <- suppressWarnings(system2(executable, c("--vanilla",
    "--default-packages=datasets,utils,grDevices,graphics,stats",
    shQuote(script), shQuote(isolated),
    shQuote(file.path(root, "fit.rds")), if (require_missing) "true" else "false"),
    stdout = TRUE, stderr = TRUE))
  cat(output, sep = "\n")
  if (!is.null(attr(output, "status")))
    stop("Fresh-session stats validation failed; see the subprocess output above.", call. = FALSE)
}
isolated_stats_session()
