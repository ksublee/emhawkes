library(emhawkes)

# S4 dispatch must not depend on methods being attached by R's startup defaults.
check_minimal_s4_session <- function() {
  libraries <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(libraries, script)), add = TRUE)
  saveRDS(.libPaths(), libraries)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(readRDS(args[1]))",
    "stopifnot(!('package:methods' %in% search()))",
    "library(emhawkes)",
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
  output <- suppressWarnings(system2(executable,
    c("--vanilla", "--default-packages=datasets,utils,grDevices,graphics,stats",
      shQuote(script), shQuote(libraries)), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(output, "status"))) {
    stop(paste(c("Minimal-session S4 dispatch failed:", output), collapse = "\n"),
         call. = FALSE)
  }
}
check_minimal_s4_session()
