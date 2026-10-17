library(emhawkes)
quiet <- function(...) suppressMessages(hfit(..., optimizer = "stats"))
close <- function(x, y, tolerance = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
}
fails <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(fit = value, warnings = warnings)
}
old <- options(emhawkes.engine = "auto")
h <- new("hspec", mu = .5, alpha = 0, beta = 1)
dt <- c(0, .2, .3, .5)
fixed <- c(FALSE, TRUE, TRUE)
control <- list(maxit = 500, reltol = 1e-12)

# Independent Poisson solution and information, with either likelihood engine.
for (engine in c("r", "cpp")) for (mode in c("numeric", "analytic")) {
  if (engine == "r" && mode == "analytic") next
  options(emhawkes.engine = engine)
  fit <- quiet(h, dt, lambda_component0 = 0, fixed = fixed,
    gradient = mode, control = control)
  stopifnot(inherits(fit, "emhawkes_stats_fit"), !inherits(fit, "maxLik"), fit$code == 0,
    fit_diagnostics(fit)$status == "ok",
    fit_diagnostics(fit)$computation$optimizer$requested == "stats",
    fit_diagnostics(fit)$computation$engine$used == engine)
  close(coef(fit), c(3, 0, 1)); close(vcov(fit)[1, 1], 3)
  close(fit$hessian[1, 1], -1/3)
  close(logLik(fit), 3 * log(3) - 3)
  stopifnot(identical(attr(logLik(fit), "df"), 1L), all(diag(vcov(fit))[-1] == 0))
  restored <- fitted_hspec(fit)
  close(logLik(restored, dt, lambda_component0 = 0), logLik(fit), 1e-12)
  info <- summary(fit)
  stopifnot(inherits(info, "summary.emhawkes_stats_fit"), is.matrix(info$coefficients),
    all(is.na(info$coefficients[-1, "Pr(>|z|)"])))
  invisible(capture.output(print(fit)))
  invisible(capture.output(print(info)))
}

# Fixed coordinates, derivative sources, disabled and raw final curvature.
for (value in list(fixed, 2:3, c("alpha1", "beta1"))) {
  fit <- quiet(h, dt, lambda_component0 = 0, fixed = value, control = control)
  close(coef(fit), c(3, 0, 1))
}
raw <- quiet(h, dt, lambda_component0 = 0, fixed = fixed,
  final_hessian = "optimizer", control = control)
close(vcov(raw)[1, 1], 3, 1e-5)
stopifnot(fit_diagnostics(raw)$computation$hessian$used == "stats_numeric")
disabled <- quiet(h, dt, lambda_component0 = 0, fixed = fixed, finalHessian = FALSE)
stopifnot(is.null(disabled$hessian), is.na(vcov(disabled)[1, 1]),
  fit_diagnostics(disabled)$computation$hessian$used == "disabled")
supplied <- quiet(h, dt, lambda_component0 = 0, fixed = fixed,
  grad = function(p) c(3/p[1]-1, 0, 0),
  hess = function(p) diag(c(-3/p[1]^2, 0, 0)), control = control)
close(vcov(supplied)[1, 1], 3)
stopifnot(fit_diagnostics(supplied)$computation$gradient$used == "user",
          fit_diagnostics(supplied)$computation$hessian$used == "user")

# Alias normalization and strict interior start; active boundary has no Wald SE.
restriction <- list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -.1)
fit <- quiet(h, dt, lambda_component0 = 0, fixed = fixed,
  constraints = restriction, control = control, outer.eps = 1e-12)
alias <- quiet(h, dt, lambda_component0 = 0, fixed = fixed,
  constraint = restriction, control = control, outer.eps = 1e-12)
stopifnot(identical(coef(fit), coef(alias)), identical(vcov(fit), vcov(alias)),
  fit$optimizer == "constrOptim", fit_diagnostics(fit)$status == "ok")
close(coef(fit)[1], 3)
boundary_h <- new("hspec", mu = 5, alpha = 0, beta = 1)
boundary <- capture(quiet(boundary_h, dt, lambda_component0 = 0, fixed = fixed,
  constraints = list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -4), control = control))
close(coef(boundary$fit)[1], 4, 1e-4)
stopifnot(is.na(vcov(boundary$fit)[1, 1]), length(boundary$warnings) > 0,
  fit_diagnostics(boundary$fit)$reason == "constraint_near_boundary")

# Bound-constrained maximum and raw-Hessian boundary protection.
for (curvature in c("stable", "optimizer")) {
  box <- capture(quiet(h, dt, lambda_component0 = 0, fixed = fixed, method = "L-BFGS-B",
    lower = c(.1, 0, 1), upper = c(2, 0, 1), final_hessian = curvature))
  close(coef(box$fit)[1], 2, 1e-12)
  stopifnot(is.na(vcov(box$fit)[1, 1]),
    fit_diagnostics(box$fit)$reason == "constraint_near_boundary")
}

# Fixed contributions in A theta + B, with an independently known quadratic optimum.
q <- new("hspec", mu = 2, alpha = .2, beta = 1)
quad <- function(object, param) -.5 * sum((param - c(1, .7, 2))^2)
fix <- c(FALSE, TRUE, FALSE)
linear <- list(ineqA = matrix(c(1, 2, 0), 1), ineqB = -1.6)
answer <- quiet(q, mylogLik = quad, fixed = fix, constraints = linear,
  control = control, outer.eps = 1e-12)
close(coef(answer), c(1.2, .2, 2), 1e-5)
fixed_only <- list(ineqA = matrix(c(0, 1, 0), 1), ineqB = -.2)
answer <- quiet(q, mylogLik = quad, fixed = fix, constraints = fixed_only, control = control)
close(coef(answer), c(1, .2, 2))
stopifnot(answer$optimizer == "optim")
all_fixed <- quiet(q, mylogLik = quad, fixed = rep(TRUE, 3))
stopifnot(all_fixed$code == 0L, all(vcov(all_fixed) == 0), identical(coef(all_fixed), c(mu1=2, alpha1=.2, beta1=1)))

# Likelihood derivative attributes and explicit derivative precedence.
tagged <- function(object, param) structure(quad(object, param), gradient = -(param-c(1,.7,2)), hessian = -diag(3))
tag <- quiet(q, mylogLik = tagged, fixed = fix, control = control)
close(coef(tag), c(1, .2, 2)); close(diag(vcov(tag)), c(1, 0, 1))
stopifnot(fit_diagnostics(tag)$computation$gradient$used == "attribute",
  fit_diagnostics(tag)$computation$hessian$used == "attribute")
precedence <- quiet(q, mylogLik = tagged, fixed = fix, control = control,
  grad = function(p) -(p-c(1,.7,2)), hess = function(p) -diag(3))
stopifnot(fit_diagnostics(precedence)$computation$gradient$used == "user",
  fit_diagnostics(precedence)$computation$hessian$used == "user")

# Correlated curvature: conditional optimum and inverse information are known.
Q <- matrix(c(4,1,2, 1,3,.5, 2,.5,5), 3)
correlated <- function(object, param) {
  delta <- param-c(1,.7,2)
  -.5 * drop(crossprod(delta, Q %*% delta))
}
for (derivatives in c(FALSE, TRUE)) {
  extra <- if (derivatives) list(grad = function(p) -drop(Q %*% (p-c(1,.7,2))),
                                  hess = function(p) -Q) else list()
  correlated_fit <- do.call(quiet, c(list(object = q, mylogLik = correlated,
    fixed = fix, control = control), extra))
  close(coef(correlated_fit), c(1.125,.2,2), 1e-6)
  close(vcov(correlated_fit)[c(1,3),c(1,3)], solve(Q[c(1,3),c(1,3)]), 1e-6)
  stopifnot(all(vcov(correlated_fit)[2,] == 0), all(vcov(correlated_fit)[,2] == 0))
}
matrix_score <- function(object, param) {
  gradient <- -(param-c(1,.7,2))
  structure(quad(object,param), gradient = rbind(gradient/2,gradient/2))
}
matrix_fit <- quiet(q, mylogLik = matrix_score, fixed = fix, control = control)
close(coef(matrix_fit), c(1,.2,2))
stopifnot(fit_diagnostics(matrix_fit)$computation$gradient$used == "attribute")

# Sharing and parameter-only functions survive the new result type.
functional <- new("hspec", mu = function(param = c(log_rate = 0)) matrix(exp(param), 2),
  alpha = matrix(0, 2, 3), beta = matrix(1, 2, 3), type_col_map = list(c(1,3), 2))
types <- c(0, 1, 2, 1)
map <- emhawkes:::build_parameter_map(functional, TRUE)
fixed_functional <- names(map$start) != "log_rate"
f <- quiet(functional, dt, types, lambda_component0 = 0, fixed = fixed_functional,
  control = control, final_hessian = "optimizer")
close(coef(f)["log_rate"], log(1.5))
restored <- fitted_hspec(f)
close(restored@mu(), rep(1.5, 2))
stopifnot(identical(restored@type_col_map, functional@type_col_map))
close(logLik(restored, dt, types, lambda_component0 = 0), logLik(f), 1e-12)

# Observed-event baseline: independent nonhomogeneous Poisson MLE and variance.
observed <- new("hspec", mu = function(param = c(m = 0), n) matrix(exp(param) * (1+n/10)),
  alpha = function(param = numeric(0)) matrix(0), beta = 1)
observed_jacobian <- function(param, n) {
  J <- matrix(0, 1, length(param), dimnames = list(NULL, names(param)))
  J[, "m"] <- exp(param["m"]) * (1+n/10); J
}
expected_m <- log(3 / sum(dt[-1] * (1+(2:4)/10)))
for (engine in c("r", "cpp")) {
  options(emhawkes.engine = engine)
  modes <- if (engine == "r") "numeric" else c("numeric", "analytic")
  for (mode in modes) {
    observed_fit <- quiet(observed, dt, lambda_component0 = 0, fixed = "beta1",
      gradient = mode, jacobian = list(mu = observed_jacobian), control = control)
    close(coef(observed_fit)["m"], expected_m)
    close(vcov(observed_fit)["m", "m"], 1/3)
    stopifnot(fit_diagnostics(observed_fit)$status == "ok")
  }
}

# Other advertised optim methods still solve a known objective.
for (method in c("CG", "Nelder-Mead")) {
  other <- capture(quiet(q, mylogLik = quad, fixed = fix, method = method,
                          control = list(maxit = 1000, reltol = 1e-12)))
  close(coef(other$fit), c(1, .2, 2), 1e-5)
  stopifnot(other$fit$code == 0L)
}

# Invalid requests fail explicitly rather than changing providers or dropping constraints.
base <- list(object = h, inter_arrival = dt, lambda_component0 = 0, fixed = fixed)
bad <- function(extra, pattern) fails(do.call(quiet, c(base, extra)), pattern)
bad(list(method = "NR"), "arg")
bad(list(final_hessian = "maxLik"), "Use final_hessian")
bad(list(iterlim = 2), "accepts only")
bad(list(control = list(iterlim = 2)), "maxit")
bad(list(control = list(fnscale = -1)), "must be positive")
bad(list(constraints = list(eqA = matrix(c(1,0,0),1), eqB = -1)), "equality constraints")
bad(list(constraints = list(ineqA = diag(2), ineqB = c(0,0))), "one column")
bad(list(constraints = list(ineqA = matrix(c(1,0,0),1), ineqB = -.5)), "strictly inside")
bad(list(lower = .1), "L-BFGS-B")
bad(list(grad = function(p) 1), "full-parameter score")
bad(list(hess = function(p) matrix(1)), "full-parameter Hessian")
bad(list(grad = function(p) setNames(rep(1,3), c("wrong","names","here"))), "full-parameter score")
bad(list(outer.eps = 1e-10), "movable linear")
bad(list(outer.iterations = 0), "positive finite integer")
fails(quiet(q, mylogLik = quad, fixed = fix,
  constraints = list(ineqA = matrix(c(0,1,0),1), ineqB = -.3)), "Fixed parameters violate")
fails(quiet(h, mylogLik = function(object, param) -Inf), "starting log likelihood")
limited <- capture(quiet(h, dt, lambda_component0 = 0, fixed = fixed, control = list(maxit = 1)))
stopifnot(limited$fit$code != 0L, !fit_diagnostics(limited$fit)$optimizer_success,
  any(grepl("normal convergence", limited$warnings, fixed = TRUE)))

# An existing maxLik installation is not needed to consume a saved stats fit.
saved <- tempfile(fileext = ".rds")
saveRDS(f, saved)
copy <- readRDS(saved)
stopifnot(identical(coef(copy), coef(f)), identical(vcov(copy), vcov(f)),
  identical(fit_diagnostics(copy), fit_diagnostics(f)))
unlink(saved)
options(old)
