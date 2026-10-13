# Compare saved runs from reproduce-fit-diagnostics.R in legacy/ and stable/.
# Rscript tools/check-hawkes-standard-errors.R <output-directory>
suppressPackageStartupMessages(library(emhawkes))
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply the directory containing legacy/ and stable/ runs.")
source("tools/diagnose-fit.R")
rows <- list()
for (name in c("multi_kernel", "synchronized", "bid_ask")) {
  legacy <- readRDS(file.path(args[1L], "legacy", name, "fit.rds"))
  fit <- readRDS(file.path(args[1L], "stable", name, "fit.rds"))
  for (field in c("estimate", "maximum", "gradient", "code", "iterations")) {
    if (!identical(fit[[field]], legacy[[field]])) stop(name, ": changed optimizer field ", field)
  }
  diagnostics <- fit_diagnostics(fit)
  stopifnot(identical(diagnostics$status, "ok"))
  p <- coef(fit)
  active <- maxLik::activePar(fit)
  fn <- function(free) {
    proposal <- p
    proposal[active] <- free
    sum(maxLik::objectiveFn(fit)(proposal))
  }
  independent <- central_fit_derivatives(fn, p[active], pmax(abs(p[active]), .1), .0003)
  information <- -independent$hessian
  stopifnot(min(eigen(information, symmetric = TRUE)$values) > 0)
  reference <- sqrt(diag(solve(information)))
  se <- sqrt(diag(vcov(fit)))[active]
  relative <- max(abs(se / reference - 1))
  stopifnot(all(is.finite(se)), relative < .002)
  comparison <- data.frame(parameter = names(p)[active], estimate = unname(p[active]),
    legacy_se = suppressWarnings(sqrt(diag(vcov(legacy))))[active],
    new_se = unname(se), independent_central_se = unname(reference))
  write.csv(comparison, file.path(args[1L], paste0(name, "-standard-errors.csv")), row.names = FALSE)
  rows[[name]] <- data.frame(case = name, coefficient_difference = 0,
    loglikelihood_difference = 0, code = fit$code, iterations = fit$iterations,
    maximum_relative_se_difference = relative)
}
write.csv(do.call(rbind, rows), file.path(args[1L], "unchanged-estimates.csv"), row.names = FALSE)
print(do.call(rbind, rows), row.names = FALSE)
