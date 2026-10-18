library(emhawkes)

# S4 dispatch must not depend on methods being attached by R's startup defaults.
check_minimal_s4_session <- function() {
  libraries <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(libraries, script)), add = TRUE)
  saveRDS(list(paths = .libPaths(), package = find.package("emhawkes"),
    version = as.character(utils::packageVersion("emhawkes"))), libraries)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "libraries <- readRDS(args[1])",
    ".libPaths(libraries$paths)",
    "stopifnot(!('package:methods' %in% search()))",
    "if (args[2] == 'namespace') loadNamespace('emhawkes') else library(emhawkes)",
    "stopifnot(identical(normalizePath(find.package('emhawkes')),",
    "                    normalizePath(libraries$package)))",
    "stopifnot(as.character(utils::packageVersion('emhawkes')) == libraries$version)",
    "stopifnot(!('package:methods' %in% search()))",
    "h <- methods::new('hspec', mu = .5, alpha = 0, beta = 1)",
    "stopifnot(length(h@dimens) == 1L, h@dimens == 1,",
    "          identical(dim(h@mu), c(1L, 1L)),",
    "          identical(dim(h@alpha), c(1L, 1L)),",
    "          identical(dim(h@beta), c(1L, 1L)))",
    "text <- capture.output(methods::show(h))",
    "stopifnot(any(grepl('1-dimensional Hawkes process', text, fixed = TRUE)))",
    "for (engine in c('r', 'cpp')) {",
    "  options(emhawkes.engine = engine)",
    "  value <- logLik(h, c(0, .2, .3, .5), lambda_component0 = 0)",
    "  stopifnot(abs(as.numeric(value) - (3 * log(.5) - .5)) < 1e-12)",
    "  stopifnot(identical(value, stats::logLik(h, c(0, .2, .3, .5), lambda_component0 = 0)))",
    "}",
    "mapping <- list(c(1L, 3L), 2L)",
    "mapped <- methods::new('hspec', mu = c(.2, .3),",
    "  alpha = matrix(.05, 2, 3), beta = matrix(1, 2, 3), type_col_map = mapping)",
    "stopifnot(mapped@dimens == 2, identical(dim(mapped@alpha), c(2L, 3L)),",
    "          identical(mapped@type_col_map, mapping))",
    "invalid <- tryCatch(methods::new('hspec', mu = -.1, alpha = 0, beta = 1), error = identity)",
    "stopifnot(inherits(invalid, 'error'))",
    "stopifnot(!('package:methods' %in% search()))",
    "cat('S4 initialization, display, validation and likelihood: OK\\n')"
  ), script)
  executable <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  for (mode in c("namespace", "attached")) {
    output <- suppressWarnings(system2(executable,
      c("--vanilla", "--default-packages=datasets,utils,grDevices,graphics,stats",
        shQuote(script), shQuote(libraries), mode), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(output, "status"))) {
      stop(paste(c(paste("Minimal-session S4 dispatch failed with", mode), output), collapse = "\n"),
           call. = FALSE)
    }
  }
}
check_minimal_s4_session()
