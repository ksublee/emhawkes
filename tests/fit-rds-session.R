# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)

# Reading a fit in a new session must register the maxLik S3 methods before
# coefficient extraction or model reconstruction, without attaching maxLik.
check_rds_session <- function() {
  fit_path <- tempfile(fileext = ".rds")
  libraries_path <- tempfile(fileext = ".rds")
  script_path <- tempfile(fileext = ".R")
  on.exit(unlink(c(fit_path, libraries_path, script_path)), add = TRUE)

  previous_options <- options(emhawkes.engine = "cpp")
  on.exit(options(previous_options), add = TRUE)
  inter_arrival <- c(0, .2, .3, .5, 1)
  specification <- methods::new("hspec", mu = .3, alpha = 0, beta = 1)
  fit <- suppressMessages(hfit(specification, inter_arrival,
    lambda_component0 = 0, N0 = 3, Nc0 = 7,
    fixed = c(FALSE, TRUE, TRUE)))
  stopifnot(fit$code == 0)
  saveRDS(fit, fit_path)
  # Use the library containing the package under test, including R CMD check.
  saveRDS(.libPaths(), libraries_path)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(readRDS(args[1]))",
    "stopifnot(!('package:methods' %in% search()))",
    "library(emhawkes)",
    "stopifnot(!('package:methods' %in% search()))",
    "fit <- readRDS(args[2])",
    "options(emhawkes.engine = 'r')",
    "calculation <- fit_diagnostics(fit)$computation",
    "stopifnot(calculation$recorded, identical(calculation$engine$used, 'cpp'),",
    "          calculation$engine$requested == 'cpp',",
    "          calculation$gradient$used == 'numeric',",
    "          calculation$hessian$used == 'maxLik_numeric')",
    "stopifnot(!('package:maxLik' %in% search()))",
    "coefficients <- stats::coef(fit)",
    "expected_names <- names(attr(fit, 'emhawkes')$parameter_map$start)",
    "stopifnot(is.numeric(coefficients), length(coefficients) == 3L,",
    "          identical(names(coefficients), expected_names),",
    "          all(is.finite(coefficients)),",
    "          abs(coefficients[['mu1']] - 2) < 1e-5)",
    "stopifnot(inherits(summary(fit), 'summary.maxLik'))",
    "covariance <- stats::vcov(fit)",
    "stopifnot(is.matrix(covariance), identical(dim(covariance), c(3L, 3L)),",
    "          all(is.finite(covariance)), covariance[1, 1] > 0)",
    "estimated <- fitted_hspec(fit)",
    "initial <- fit_initial_state(fit)",
    "stopifnot(methods::is(estimated, 'hspec'),",
    "          length(estimated@dimens) == 1L, estimated@dimens == 1,",
    "          identical(dim(estimated@mu), c(1L, 1L)),",
    "          isTRUE(all.equal(initial,",
    "            list(lambda_component0 = matrix(0), N0 = 3, Nc0 = 7))))",
    "replayed <- do.call(logLik, c(list(object = estimated,",
    "  inter_arrival = c(0, .2, .3, .5, 1)), initial))",
    "stopifnot(isTRUE(all.equal(as.numeric(stats::logLik(fit)),",
    "                          as.numeric(replayed), tolerance = 1e-8)))"
  ), script_path)

  executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  rscript <- file.path(R.home("bin"), executable)
  output <- suppressWarnings(system2(rscript,
    args = c("--vanilla", "--default-packages=datasets,utils,grDevices,graphics,stats",
             shQuote(script_path), shQuote(libraries_path),
             shQuote(fit_path)), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(c("Fresh-session fit restoration failed:", output), collapse = "\n"),
         call. = FALSE)
  }
}

check_rds_session()

}
