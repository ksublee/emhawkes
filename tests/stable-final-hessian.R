# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)
h <- new("hspec", mu = 1, alpha = 0, beta = 1)
dt <- c(0, .5, 1.5, .25)
fixed <- c(FALSE, TRUE, TRUE)
fields <- c("estimate", "maximum", "gradient", "code", "iterations")
quiet <- function(args) suppressWarnings(suppressMessages(do.call(hfit, args)))
arguments <- list(object = h, inter_arrival = dt, lambda_component0 = 0,
                  fixed = fixed, final_hessian = "stable")
# Explicit TRUE retains the previous two-pass stable calculation.
compare <- function(args) {
  set.seed(240)
  old <- quiet(c(args, list(finalHessian = TRUE)))
  set.seed(240)
  new <- quiet(args)
  for (name in fields) stopifnot(identical(old[[name]], new[[name]]))
  stopifnot(identical(vcov(old), vcov(new)),
    identical(maxLik::hessian(old), maxLik::hessian(new)),
    identical(summary(old)$estimate, summary(new)$estimate),
    isTRUE(fit_diagnostics(new)$computation$hessian$maxLik_skipped),
    !fit_diagnostics(old)$computation$hessian$maxLik_skipped,
    is.null(fit_diagnostics(new)$original_hessian),
    !is.null(fit_diagnostics(old)$original_hessian),
    is.null(fit_diagnostics(new)$computation$hessian$finalHessian),
    isTRUE(fit_diagnostics(old)$computation$hessian$finalHessian))
  list(old = old, new = new)
}
previous <- options(emhawkes.engine = "auto")
for (engine in c("auto", "r")) {
  options(emhawkes.engine = engine)
  for (method in c("BFGS", "CG", "conjugate-gradient", "NM", "Nelder-Mead", "SANN")) {
    a <- c(arguments, list(method = method, control = list(iterlim = 50)))
    z <- compare(a)
    stopifnot(fit_diagnostics(z$new)$status == "ok")
  }
}
options(emhawkes.engine = "auto")
invisible(compare(c(arguments, list(gradient = "analytic"))))
invisible(compare(c(arguments, list(grad = function(param) c(3 / param[1] - sum(dt), 0, 0)))))

# Count actual likelihood work, including missing-Jacobian fallback.
evaluations <- 0L
impact <- function(n, ...) {
  evaluations <<- evaluations + 1L
  matrix(0)
}
counted <- arguments
counted$object <- new("hspec", mu = 1, alpha = 0, beta = 1, impact = impact)
run <- function(a) {
  evaluations <<- 0L
  fit <- quiet(a)
  list(fit = fit, evaluations = evaluations)
}
a <- run(c(counted, list(finalHessian = TRUE)))
b <- run(counted)
stopifnot(a$evaluations > b$evaluations,
  identical(a$fit[fields], b$fit[fields]), identical(vcov(a$fit), vcov(b$fit)),
  fit_diagnostics(b$fit)$method == "scaled_central_Richardson",
  nzchar(fit_diagnostics(b$fit)$score_fallback_reason))
plain <- counted; plain$final_hessian <- NULL
a <- run(plain)
plain$final_hessian <- "maxLik"
b <- run(plain)
stopifnot(identical(a$fit[fields], b$fit[fields]), a$evaluations == b$evaluations,
  !fit_diagnostics(a$fit)$computation$hessian$maxLik_skipped)

# Preserve explicit overrides, supplied Hessians, constraints, Newton methods,
# custom objectives and residual likelihoods.
not_skipped <- function(a) {
  x <- quiet(a)
  stopifnot(!fit_diagnostics(x)$computation$hessian$maxLik_skipped)
  invisible(x)
}
hess <- function(param) diag(c(-3 / param[1]^2, 0, 0))
not_skipped(c(arguments, list(hess = hess)))
for (method in c("NR", "BFGSR", "bfgs-r"))
  not_skipped(c(arguments, list(method = method)))
for (flag in list(TRUE, 1, FALSE, 0, "BHHH"))
  not_skipped(c(arguments, list(finalHessian = flag)))
# maxLik also accepts abbreviated names through ...; keep those overrides.
for (name in c("finalH", "finalHess"))
  not_skipped(c(arguments, setNames(list(TRUE), name)))
for (constraints in list(
    list(ineqA = matrix(c(1, 0, 0), 1), ineqB = 0),
    list(eqA = matrix(c(0, 1, 0), 1), eqB = 0)))
  not_skipped(c(arguments, list(constraints = constraints)))
custom <- arguments
custom$mylogLik <- function(object) -(object@mu - 1)^2
not_skipped(custom)
residual <- arguments
residual$object <- new("hspec", mu = 1, alpha = 0, beta = 1,
  dresidual = function(x) dexp(x), presidual = function(q) pexp(q),
  qresidual = function(p) qexp(p))
not_skipped(residual)

# Unidentified curvature stays unavailable, without using an optimizer matrix.
singular <- arguments; singular$fixed <- c(FALSE, TRUE, FALSE)
a <- quiet(c(singular, list(finalHessian = TRUE)))
b <- quiet(singular)
stopifnot(identical(a[fields], b[fields]), identical(vcov(a), vcov(b)),
  fit_diagnostics(b)$status == "unavailable",
  fit_diagnostics(b)$computation$hessian$used == "unavailable",
  all(is.na(summary(b)$estimate[c(1, 3), "Std. error"])))
file <- tempfile(fileext = ".rds")
saveRDS(b, file)
restored <- readRDS(file); unlink(file)
stopifnot(identical(fit_diagnostics(b), fit_diagnostics(restored)),
  is.null(fit_diagnostics(restored)$original_hessian))
options(previous)

}
