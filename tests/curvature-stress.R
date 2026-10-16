library(emhawkes)
previous <- options(emhawkes.engine = "auto")
fit_case <- function(h, dt, engine = "auto", initial = 0, ...) {
  old <- options(emhawkes.engine = engine)
  on.exit(options(old))
  notices <- character()
  fit <- withCallingHandlers(suppressMessages(hfit(h, dt, lambda_component0 = initial,
    reduced = FALSE, final_hessian = "stable", control = list(iterlim = 300, reltol = 1e-11), ...)),
    warning = function(w) { notices <<- c(notices, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(fit = fit, diagnostics = fit_diagnostics(fit), notices = notices)
}
model <- function(beta, alpha = c(.18, .22)) new("hspec", mu = .5,
  alpha = matrix(alpha, 1), beta = matrix(beta, 1), type_col_map = list(1:2))
set.seed(1710)
dt <- suppressMessages(hsim(model(c(1, 1.01)), size = 401, lambda_component0 = 0))$inter_arrival
fixed <- c(FALSE, FALSE, FALSE, TRUE, TRUE)

# Independent fixed-beta information: the likelihood is linear inside log(lambda).
# This pair-delay design does not use the package's likelihood or score recursion.
linear_reference <- function(dt, beta, parameters) {
  t <- cumsum(dt)[-1L]; lag <- outer(t, t, "-")
  Z <- matrix(1, length(t), 1 + length(beta))
  for (j in seq_along(beta)) {
    W <- matrix(0, length(t), length(t)); earlier <- lower.tri(W)
    W[earlier] <- exp(-beta[j] * lag[earlier])
    Z[, j + 1L] <- rowSums(W)
  }
  lambda <- drop(Z %*% parameters)
  scale <- pmax(abs(parameters), .1)
  M <- sweep(sweep(Z, 1L, lambda, "/"), 2L, scale, "*")
  sv <- svd(M)
  C <- (sv$v %*% diag(1 / sv$d^2, length(sv$d)) %*% t(sv$v)) * outer(scale, scale)
  list(covariance = C, ratio = (min(sv$d) / max(sv$d))^2)
}
close <- fit_case(model(c(1, 1.0002)), dt, fixed = fixed)
r <- linear_reference(dt, c(1, 1.0002), coef(close$fit)[1:3])
stopifnot(close$diagnostics$status == "ok", close$diagnostics$scaled_condition > 1e8,
  "weak_curvature" %in% close$diagnostics$warnings,
  any(grepl("weak curvature", close$notices, fixed = TRUE)),
  max(abs(sqrt(diag(vcov(close$fit)))[1:3] / sqrt(diag(r$covariance)) - 1)) < .001)

# Equal rates identify only the sum of two separate kernel amplitudes.
equal <- fit_case(model(c(1, 1)), dt, fixed = fixed)
p <- coef(equal$fit); direction <- c(0, .05, -.05, 0, 0)
fn <- maxLik::objectiveFn(equal$fit)
stopifnot(length(p) == 5L, equal$diagnostics$status == "unavailable",
  all(is.na(diag(vcov(equal$fit))[1:3])), all(diag(vcov(equal$fit))[4:5] == 0),
  abs(fn(p + direction) - fn(p - direction)) < 1e-9,
  linear_reference(dt, c(1, 1), p[1:3])$ratio < 1e-20)
near <- fit_case(model(c(1, 1.00001)), dt, fixed = fixed)
stopifnot(near$diagnostics$status == "unavailable", all(is.na(diag(vcov(near$fit))[1:3])))

# Explicit R likelihood differences also match the independent covariance.
r_fit <- fit_case(model(c(1, 1.5)), dt, engine = "r", fixed = fixed)
r <- linear_reference(dt, c(1, 1.5), coef(r_fit$fit)[1:3])
stopifnot(r_fit$diagnostics$status == "ok", r_fit$diagnostics$curvature_source == "likelihood",
  max(abs(sqrt(diag(vcov(r_fit$fit)))[1:3] / sqrt(diag(r$covariance)) - 1)) < .003)

# Genuine Hawkes boundaries, with nonzero excitation in the baseline-rate case.
regular <- c(0, rep(1, 40))
a <- fit_case(new("hspec", mu = 1, alpha = .1, beta = 1), regular,
  fixed = c(FALSE, FALSE, TRUE), constraints = list(ineqA = matrix(c(0, 1, 0), 1), ineqB = 0))
b <- fit_case(new("hspec", mu = 2.5, alpha = .2, beta = 1), regular,
  fixed = c(FALSE, TRUE, TRUE), constraints = list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -2))
for (z in list(a, b)) stopifnot(z$diagnostics$status == "unavailable",
  z$diagnostics$reason == "constraint_near_boundary",
  all(is.na(diag(vcov(z$fit))[maxLik::activePar(z$fit)])),
  any(grepl("standard errors are unavailable", z$notices, fixed = TRUE)))

# At zero excitation, a known nonzero initial transient can still identify beta.
for (initial in c(0, 30)) {
  set.seed(411)
  sample <- suppressMessages(hsim(new("hspec", mu = .5, alpha = 0, beta = 1),
    size = 151, lambda_component0 = initial))
  z <- fit_case(new("hspec", mu = .5, alpha = 0, beta = .8), sample$inter_arrival,
    initial = initial, fixed = c(FALSE, TRUE, FALSE),
    constraints = list(ineqA = diag(3)[c(1, 3), ], ineqB = c(-.001, -.05)))
  if (initial == 0) stopifnot(z$diagnostics$status == "unavailable") else
    stopifnot(z$diagnostics$status == "ok", all(is.finite(diag(vcov(z$fit)))) )
}

# Nonconvergence must not be hidden by a failed curvature calculation.
early <- fit_case(new("hspec", mu = 1, alpha = 0, beta = 1), c(0, .5, 1.5, .25),
  fixed = c(FALSE, TRUE, FALSE), iterlim = 1)
stopifnot(early$diagnostics$status == "unavailable", !early$diagnostics$optimizer_success,
  "optimizer_not_converged" %in% early$diagnostics$warnings,
  any(grepl("normal convergence", early$notices, fixed = TRUE)),
  any(grepl("standard errors are unavailable", early$notices, fixed = TRUE)))
options(previous)
