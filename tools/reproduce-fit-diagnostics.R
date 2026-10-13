# Run only the assignments needed to reproduce the four affected example fits.
# Rscript tools/reproduce-fit-diagnostics.R <output-directory> [cases] [auto|maxLik]
suppressPackageStartupMessages(library(emhawkes))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 3L) stop("Supply an output directory, optional cases, and optional final-Hessian choice.")
cases <- list(
  multi_kernel = list(heading = "### Multi-kernel model", fit = "fit_mk"),
  synchronized = list(heading = "### Synchronized intensity model", fit = "fit_sy"),
  bid_ask = list(heading = "### Bid-ask price model", fit = "mle_ba"),
  frpp = list(heading = "### Two-dimensional model", fit = "mle_dh"))
selected <- if (length(args) >= 2L) strsplit(args[2L], ",", fixed = TRUE)[[1L]] else names(cases)
if (!length(selected) || any(!selected %in% names(cases))) stop("Unknown diagnostic case.")
final_hessian <- if (length(args) == 3L) match.arg(args[3L], c("maxLik", "auto")) else "maxLik"
source("tools/diagnose-fit.R")
lines <- readLines("vignettes_ext/extended_example.Rmd", warn = FALSE)
headings <- grep("^#{2,3} ", lines)
chunks <- grep("^```\\{r", lines)
dir.create(args[1L], recursive = TRUE, showWarnings = FALSE)
for (name in selected) {
  configuration <- cases[[name]]
  first <- match(configuration$heading, lines)
  if (is.na(first)) stop("Example heading is unavailable: ", configuration$heading)
  following <- headings[headings > first]
  last <- if (length(following)) following[1L] else length(lines) + 1L
  environment <- new.env(parent = globalenv())
  environment$hfit <- function(...) {
    arguments <- list(...)
    arguments$final_hessian <- final_hessian
    do.call(emhawkes::hfit, arguments)
  }
  cat("Reproducing", name, "\n")
  for (opening in chunks[chunks > first & chunks < last]) {
    closing <- opening + which(lines[(opening + 1L):length(lines)] == "```")[1L]
    expressions <- parse(text = lines[(opening + 1L):(closing - 1L)])
    for (expression in expressions) {
      if (is.call(expression) && as.character(expression[[1L]])[1L] %in%
          c("<-", "=", "set.seed", "library")) eval(expression, envir = environment)
      if (exists(configuration$fit, envir = environment, inherits = FALSE)) break
    }
    if (exists(configuration$fit, envir = environment, inherits = FALSE)) break
  }
  fit <- environment[[configuration$fit]]
  if (is.null(fit)) stop("The example did not produce its expected fit.")
  destination <- file.path(args[1L], name)
  dir.create(destination, showWarnings = FALSE)
  saveRDS(fit, file.path(destination, "fit.rds"))
  result <- diagnose_fit(fit)
  write.csv(result$curvature, file.path(destination, "curvature.csv"), row.names = FALSE)
  write.csv(result$standard_errors, file.path(destination, "standard-errors.csv"), row.names = FALSE)
  saveRDS(result, file.path(destination, "diagnostics.rds"))
  print(result$curvature, row.names = FALSE)
}
capture.output(sessionInfo(), file = file.path(args[1L], "session-info.txt"))
