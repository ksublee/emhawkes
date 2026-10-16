# Rscript --vanilla tools/check-stable-stress.R LIB OUT [BASELINE_OUT]
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(2L, 3L))
.libPaths(c(args[1], .libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(args[1])))
options(emhawkes.engine = "auto")
out <- args[2]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
quiet <- function(expr) suppressMessages(expr)
simulate <- function(h, events, seed, initial = 0) {
  set.seed(seed)
  quiet(hsim(h, size = events + 1L, lambda_component0 = initial))$inter_arrival
}
model <- function(beta = c(1, 1.01), alpha = c(.18, .22), mu = .5) {
  new("hspec", mu = mu, alpha = matrix(alpha, 1), beta = matrix(beta, 1),
      type_col_map = if (length(beta) > 1L) list(seq_along(beta)) else NULL)
}
cases <- list()
add <- function(id, family, h, dt, initial = 0, ...) {
  cases[[id]] <<- list(id = id, family = family, arguments = c(list(object = h,
    inter_arrival = dt, lambda_component0 = initial, reduced = FALSE,
    control = list(iterlim = 300, reltol = 1e-11)), list(...)))
}
dt <- simulate(model(), 400, 1710)
for (gap in c(0, 1e-5, 1e-4, .001, .01, .1, .5))
  add(paste0("known-beta-gap-", gap), "known_beta", model(c(1, 1 + gap)), dt,
      fixed = c(FALSE, FALSE, FALSE, TRUE, TRUE))
for (beta in list(c(1, 1.01), c(.5, 2), c(2, .5)))
  add(paste0("free-beta-start-", paste(beta, collapse = "-")), "free_beta", model(beta), dt,
      constraints = list(ineqA = diag(5), ineqB = c(-.001, 0, 0, -.05, -.05)))
for (n in c(1, 2, 5))
  add(paste0("few-", n), "few_events", model(1, .05, .7), c(0, rep(1, n)),
      constraints = list(ineqA = diag(3), ineqB = c(-.001, 0, -.05)))
for (n in c(12, 40, 150))
  add(paste0("weak-events-", n), "weak_excitation", model(1, .05, .7),
      simulate(model(1, .01, .7), n, 712 + n),
      constraints = list(ineqA = diag(3), ineqB = c(-.001, 0, -.05)))
for (seed in 1:12)
  add(paste0("weak-seed-", seed), "weak_excitation", model(1, .05, .7),
      simulate(model(1, .01, .7), 40, seed),
      constraints = list(ineqA = diag(3), ineqB = c(-.001, 0, -.05)))
for (initial in c(0, 30))
  add(paste0("zero-alpha-initial-", initial), "initial_state", model(.8, 0),
      simulate(model(1, 0), 150, 411, initial), initial = initial,
      fixed = c(FALSE, TRUE, FALSE),
      constraints = list(ineqA = diag(3)[c(1, 3), ], ineqB = c(-.001, -.05)))
add("boundary-alpha", "boundary", model(1, .1, 1), c(0, rep(1, 40)),
    fixed = c(FALSE, FALSE, TRUE), constraints = list(ineqA = matrix(c(0, 1, 0), 1), ineqB = 0))
add("boundary-mu", "boundary", model(1, .2, 2.5), c(0, rep(1, 40)),
    fixed = c(FALSE, TRUE, TRUE), constraints = list(ineqA = matrix(c(1, 0, 0), 1), ineqB = -2))
add("interior", "interior", model(1, .2), simulate(model(1, .4), 400, 1710),
    fixed = c(FALSE, FALSE, TRUE), constraints = list(ineqA = diag(3)[1:2, ], ineqB = c(-.001, 0)))
add("early-close", "early_stop", model(c(2, .5)), dt, iterlim = 1)
add("early-unidentified", "early_stop", model(1, 0, 1), c(0, .5, 1.5, .25),
    fixed = c(FALSE, TRUE, FALSE), iterlim = 1)

# Independent direct-pair R calculation, without the package event recursion.
# Parameters are (mu, alpha[1:K], beta[1:K]); initial components are constants.
reference <- function(dt, p, initial, active) {
  K <- (length(p) - 1L) / 2L; ia <- 1L + seq_len(K); ib <- 1L + K + seq_len(K)
  t <- cumsum(dt)[-1L]; n <- length(t); T <- tail(t, 1L)
  lag <- outer(t, t, "-"); before <- lower.tri(lag)
  G <- matrix(1, n, length(p)); lambda <- rep(p[1], n)
  second <- exposure <- first_exposure <- second_exposure <- vector("list", K)
  initial <- rep(initial, length.out = K)
  FT <- F2T <- numeric(K)
  for (j in seq_len(K)) {
    b <- p[ib[j]]; a <- p[ia[j]]
    W <- matrix(0, n, n); W[before] <- exp(-b * lag[before])
    S <- rowSums(W); S1 <- -rowSums(lag * W); S2 <- rowSums(lag^2 * W)
    transient <- initial[j] * exp(-b * t)
    lambda <- lambda + a * S + transient
    G[, ia[j]] <- S; G[, ib[j]] <- a * S1 - t * transient
    u <- b * (T - t)
    exposure[[j]] <- sum(-expm1(-u) / b)
    first_exposure[[j]] <- sum((exp(-u) * (1 + u) - 1) / b^2)
    second_exposure[[j]] <- sum((2 - exp(-u) * (u^2 + 2 * u + 2)) / b^3)
    v <- b * T; FT[j] <- -expm1(-v) / b
    F2T[j] <- (2 - exp(-v) * (v^2 + 2 * v + 2)) / b^3
    second[[j]] <- list(cross = S1, beta = a * S2 + t^2 * transient)
  }
  stopifnot(all(lambda > 0))
  H <- -crossprod(sweep(G, 1L, lambda, "/"))
  for (j in seq_len(K)) {
    z <- sum(second[[j]]$cross / lambda) - first_exposure[[j]]
    H[ia[j], ib[j]] <- H[ia[j], ib[j]] + z
    H[ib[j], ia[j]] <- H[ib[j], ia[j]] + z
    H[ib[j], ib[j]] <- H[ib[j], ib[j]] + sum(second[[j]]$beta / lambda) -
      p[ia[j]] * second_exposure[[j]] - initial[j] * F2T[j]
  }
  loglik <- sum(log(lambda)) - p[1] * T - sum(p[ia] * unlist(exposure)) - sum(initial * FT)
  scale <- pmax(abs(p[active]), .1)
  info <- -H[active, active, drop = FALSE] * outer(scale, scale)
  # For known beta the Hessian is a Gram matrix. SVD avoids squaring its
  # condition number in the reference inversion of nearly coincident kernels.
  if (!any(active[ib])) {
    M <- sweep(sweep(G[, active, drop = FALSE], 1L, lambda, "/"), 2L, scale, "*")
    sv <- svd(M); condition <- (max(sv$d) / min(sv$d))^2
    positive <- min(sv$d)^2 > 1e-10 * max(sv$d)^2
    C <- if (positive) (sv$v %*% diag(1 / sv$d^2, length(sv$d)) %*% t(sv$v)) * outer(scale, scale) else NULL
  } else {
    ev <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
    positive <- min(ev) > 1e-10 * max(abs(ev))
    condition <- if (positive) max(ev) / min(ev) else Inf
    C <- if (positive) solve(info) * outer(scale, scale) else NULL
  }
  list(loglik = loglik, hessian = H[active, active, drop = FALSE],
       covariance = C, condition = condition, positive = positive)
}
fields <- c("estimate", "maximum", "gradient", "code", "iterations")
rows <- coefficients <- warnings <- slices <- comparisons <- list()
for (id in names(cases)) {
  case <- cases[[id]]; a <- case$arguments; notices <- character()
  tm <- system.time(fit <- withCallingHandlers(quiet(do.call(hfit, c(a, list(final_hessian = "stable")))),
    warning = function(w) { notices <<- c(notices, conditionMessage(w)); invokeRestart("muffleWarning") }))
  d <- fit_diagnostics(fit); p <- coef(fit); active <- maxLik::activePar(fit)
  r <- reference(a$inter_arrival, p, a$lambda_component0, active)
  stopifnot(abs(r$loglik - as.numeric(logLik(fit))) < 1e-8 * max(1, abs(r$loglik)))
  se <- sqrt(diag(vcov(fit)))
  se_error <- h_error <- NA_real_
  if (d$status == "ok") {
    stopifnot(r$positive)
    se_error <- max(abs(se[active] / sqrt(diag(r$covariance)) - 1))
    h_error <- norm(maxLik::hessian(fit)[active, active, drop = FALSE] - r$hessian, "F") /
      max(norm(r$hessian, "F"), .Machine$double.eps)
    stopifnot(se_error < .01, h_error < .001)
  } else stopifnot(all(is.na(se[active])))
  z <- list(value = fit[fields], covariance = vcov(fit), hessian = maxLik::hessian(fit),
    diagnostics = d, notices = notices, reference = r)
  if (length(args) > 2L) {
    old <- readRDS(file.path(args[3], paste0(id, ".rds")))
    stopifnot(identical(old$value, z$value), identical(old$covariance, z$covariance),
      identical(old$hessian, z$hessian), identical(old$diagnostics$status, d$status),
      identical(old$diagnostics$reason, d$reason), identical(old$diagnostics$checks, d$checks))
    if (isFALSE(d$optimizer_success)) {
      stopifnot("optimizer_not_converged" %in% d$warnings,
        any(grepl("normal convergence", notices, fixed = TRUE)))
    } else stopifnot(identical(old$notices, notices))
    comparisons[[id]] <- data.frame(id, optimizer_fields_identical = TRUE,
      hessian_identical = TRUE, covariance_identical = TRUE, curvature_checks_identical = TRUE,
      before_warning_count = length(old$notices), after_warning_count = length(notices))
  }
  # Conditional likelihood slices hold all other fitted parameters fixed.
  # They are not profile likelihoods and do not establish interval coverage.
  offsets <- if (case$family == "known_beta") c(-.1, -.05, 0, .05, .1) else
    if (case$family == "weak_excitation" && d$status == "ok")
      unique(c(0, p[2], p[2] + se[2], p[2] + 2 * se[2])) else numeric(0)
  for (offset in offsets) {
    proposal <- p
    if (case$family == "known_beta") {
      proposal[2:3] <- p[2:3] + c(offset, -offset)
      variable <- "amplitude_transfer"
    } else {
      proposal[2] <- offset; variable <- "alpha"
    }
    slices[[length(slices) + 1L]] <- data.frame(id, variable, value = offset,
      loglik_difference = as.numeric(maxLik::objectiveFn(fit)(proposal)) - as.numeric(logLik(fit)))
  }
  saveRDS(z, file.path(out, paste0(id, ".rds")))
  rows[[id]] <- data.frame(id, family = case$family, events = length(a$inter_arrival) - 1L,
    free_parameters = sum(active), loglik = as.numeric(logLik(fit)), code = fit$code, status = d$status,
    reason = if (is.null(d$reason)) "" else d$reason,
    condition = if (is.null(d$scaled_condition)) NA_real_ else d$scaled_condition,
    reference_condition = r$condition, reference_positive = r$positive,
    relative_se_error = se_error, relative_hessian_error = h_error,
    free_se_available = all(is.finite(se[active])), warnings = paste(notices, collapse = " | "),
    optimizer_warning = any(grepl("normal convergence", notices, fixed = TRUE)),
    elapsed = unname(tm[["elapsed"]]))
  coefficients[[id]] <- data.frame(id, parameter = names(p), estimate = unname(p),
    active = active, se = unname(se))
  for (notice in notices) warnings[[length(warnings) + 1L]] <- data.frame(id, warning = notice)
  cat(id, d$status, d$reason, "relative SE error", se_error, "\n"); flush.console()
}
write.csv(do.call(rbind, rows), file.path(out, "cases.csv"), row.names = FALSE)
write.csv(do.call(rbind, coefficients), file.path(out, "coefficients.csv"), row.names = FALSE)
write.csv(do.call(rbind, warnings), file.path(out, "warnings.csv"), row.names = FALSE)
write.csv(do.call(rbind, slices), file.path(out, "conditional-slices.csv"), row.names = FALSE)
if (length(comparisons)) write.csv(do.call(rbind, comparisons), file.path(out, "unchanged-results.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "session-info.txt"))
cat("Cases:", length(cases), "\n")
