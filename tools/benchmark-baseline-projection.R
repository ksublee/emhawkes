# Rscript --vanilla tools/benchmark-baseline-projection.R LIB FIXTURES OUT ID OP ROUND
# FIXTURES is the output of benchmark-gradient-scaling.R prepare.
# OP: path, score, cold, fit_numeric, fit_analytic. ROUND=0 profiles memory.
a <- commandArgs(TRUE)
lib <- a[1]; fixtures <- a[2]; out <- a[3]; id <- a[4]; op <- a[5]; round <- as.integer(a[6])
.libPaths(c(lib, .libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
options(emhawkes.engine = "auto")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows" && nzchar(Sys.getenv("EMHAWKES_BENCH_AFFINITY"))) {
  mask <- Sys.getenv("EMHAWKES_BENCH_AFFINITY")
  stopifnot(grepl("^[1-9][0-9]*$", mask))
  actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf(
    "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), mask))), stdout = TRUE)
  stopifnot(trimws(actual) == mask)
}
x <- readRDS(file.path(fixtures, "fixtures", paste0(sub("-mixedmu$", "", id), ".rds")))
if (grepl("-mixedmu$", id)) {
  # A full-rank, row-stochastic reparameterization. At the common generating
  # parameters it preserves the simulated model but makes mu derivatives dense.
  transformed <- local({
    old <- x$arguments$object@mu; env <- environment(old)
    weight <- env$mu_weight; group <- env$group; mu_names <- env$mu_names
    d <- length(weight); stopifnot(length(mu_names) == d)
    mix <- diag(.8, d) + matrix(.2 / d, d, d)
    mu <- function(param, n, N, ...) matrix(weight * as.vector(mix %*% exp(param[group])) *
      (1 + .08 * (N[max(1L, n - 1L), ] %% 3)), d)
    formals(mu)$param <- formals(old)$param
    jac <- function(param, n, N, ...) {
      J <- matrix(0, d, length(param), dimnames = list(NULL, names(param)))
      part <- sweep(mix, 2, exp(param[mu_names[group]]), "*")
      J[, mu_names] <- part * (weight * (1 + .08 * (N[max(1L, n - 1L), ] %% 3)))
      J
    }
    list(mu = mu, jac = jac)
  })
  x$arguments$object@mu <- transformed$mu; x$jacobian$mu <- transformed$jac
}
if (op == "path") {
  # Recycle observed data to study storage at 100k events, not parameter recovery.
  for (name in c("inter_arrival", "type", "mark"))
    x$arguments[[name]] <- c(x$arguments[[name]][1], rep(x$arguments[[name]][-1], length.out = 100000L))
}
fit <- startsWith(op, "fit_")
if (fit) for (name in c("inter_arrival", "type", "mark"))
  x$arguments[[name]] <- head(x$arguments[[name]], if (op == "fit_numeric") 751L else 1001L)
setup <- function(x) {
  b <- x$arguments
  obs <- emhawkes:::normalize_observation(b$object, b$inter_arrival, b$type, b$mark)
  map <- emhawkes:::build_parameter_map(b$object, TRUE)
  list(map = map, obs = obs, objective = emhawkes:::make_hawkes_objective(b$object, obs,
    b$lambda_component0, map, jacobian = x$jacobian))
}
warm <- x
for (name in c("inter_arrival", "type", "mark")) warm$arguments[[name]] <- head(warm$arguments[[name]], 51L)
w <- setup(warm); invisible(w$objective$score(w$map$start))
if (fit) {
  warm$arguments$control <- list(iterlim = 2)
  invisible(suppressWarnings(suppressMessages(do.call(hfit, c(warm$arguments,
    list(gradient = if (op == "fit_numeric") "numeric" else "analytic", jacobian = x$jacobian, final_hessian = "stable"))))))
}
s <- setup(x); counter <- 0L
if (op == "path") {
  plan <- emhawkes:::hawkes_component_plan(x$arguments$object, s$obs, s$map, x$jacobian, TRUE)
  p <- s$map$start; values <- plan$values(p)
}
notices <- character()
operation <- function() {
  if (fit) return(suppressMessages(do.call(hfit, c(x$arguments,
    list(gradient = if (op == "fit_numeric") "numeric" else "analytic", jacobian = x$jacobian, final_hessian = "stable")))))
  if (op == "path") return(plan$baseline_jacobian(p, values))
  counter <<- counter + 1L; q <- s$map$start; q[1] <- q[1] + counter * 1e-7
  if (op == "score") ll <- s$objective$loglik(q)
  g <- if (op == "cold") s$objective$score(q) else s$objective$gradient(q)
  list(parameters = q, loglikelihood = s$objective$loglik(q), gradient = g)
}
invoke <- function() withCallingHandlers(operation(), warning = function(w) {
  notices <<- c(notices, conditionMessage(w)); invokeRestart("muffleWarning")
})
rss <- function() if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
  c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))), stdout = TRUE)) else NA_real_
gc(); before <- gc(reset = TRUE); rss_before <- rss(); memory <- NULL
if (round == 0L) {
  profile <- tempfile(); Rprofmem(profile)
  value <- tryCatch(invoke(), finally = Rprofmem(NULL))
  after <- gc(); rss_after <- rss()
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(profile, warn = FALSE))))
  bytes <- bytes[is.finite(bytes)]; unlink(profile)
  memory <- c(allocated_bytes = sum(bytes), largest_bytes = max(c(0, bytes)),
    heap_peak_mib = sum(after[, 6]), rss_before = rss_before, rss_peak = rss_after)
  seconds <- NA_real_
} else seconds <- unname(system.time(value <- invoke())[["elapsed"]])
if (fit) {
  diagnostics <- fit_diagnostics(value)
  value <- list(estimate = coef(value), maximum = as.numeric(logLik(value)), gradient = value$gradient,
    hessian = value$hessian, covariance = vcov(value), se = sqrt(diag(vcov(value))),
    code = value$code, iterations = value$iterations, curvature = diagnostics$status,
    computation = diagnostics$computation)
}
if (op == "path") value <- list(bytes = as.numeric(object.size(value)),
  compact = is.list(value), dimensions = lapply(if (is.list(value)) value else list(path = value), dim))
saveRDS(list(id = id, operation = op, round = round, seconds = seconds, memory = memory,
  value = value, warnings = unique(notices), events = length(x$arguments$inter_arrival) - 1L),
  file.path(out, paste(id, op, round, "rds", sep = ".")))
cat(id, op, round, seconds, "\n")
