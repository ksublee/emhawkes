# Rscript --vanilla tools/benchmark-function-replay.R LIB OUT [parameters|observed]
# Same installed source, public logLik calls, numerical likelihood and inference.
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(args[1])))
out <- args[2]
group <- if (length(args) >= 3L) match.arg(args[3], c("parameters", "observed")) else "parameters"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
if (affinity != "scheduler") {
  stopifnot(.Platform$OS.type == "windows", grepl("^[1-9][0-9]*$", affinity))
  actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command",
    shQuote(sprintf("$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), affinity))), stdout = TRUE)
  stopifnot(length(actual) == 1L, trimws(actual) == affinity)
}
wrap <- function(value) {
  nr <- nrow(value); nc <- ncol(value)
  p <- setNames(as.vector(value), paste0("v", seq_along(value)))
  force(nr); force(nc); force(p)
  function(param = p) matrix(param, nr, nc)
}
# Both flocking interfaces implement the same rule. Rewriting the indexing in
# the example avoids extracting a whole count column for every single event.
flocking_interface <- function(object, scalar) {
  fun <- object@impact
  changed <- 0L
  rewrite <- function(expr) {
    if (!is.call(expr)) return(expr)
    column <- NULL
    if (identical(expr[[1]], as.name("["))) {
      if (length(expr) == 3L && is.call(expr[[2]]) &&
          identical(expr[[2]][[1]], as.name("[")) && identical(expr[[2]][[2]], as.name("N")) &&
          length(expr[[2]]) == 4L && is.character(expr[[2]][[4]])) column <- expr[[2]][[4]]
      if (length(expr) == 4L && identical(expr[[2]], as.name("N")) &&
          identical(expr[[3]], as.name("n")) && is.character(expr[[4]])) column <- expr[[4]]
    }
    if (!is.null(column)) {
      changed <<- changed + 1L
      if (scalar) { expr <- quote(N[n, "column"]); expr[[4]] <- column } else {
        expr <- quote(N[, "column"][n]); expr[[2]][[4]] <- column
      }
      return(expr)
    }
    for (i in seq_along(expr)[-1L]) expr[[i]] <- rewrite(expr[[i]])
    expr
  }
  body(fun) <- rewrite(body(fun))
  stopifnot(changed == 4L)
  object@impact <- fun
  object
}
if (group == "observed") {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  destination <- file.path(out, "prepare")
  status <- system2(rscript, c("--vanilla", "tools/benchmark-impact-numeric.R", "prepare",
    shQuote(args[1]), shQuote(destination)))
  stopifnot(status == 0L)
  flocking <- readRDS(file.path(destination, "fixtures.rds"))$flocking$arguments$object
  set.seed(1107)
  flocking_data <- hsim(flocking_interface(flocking, TRUE), size = 10001L)
}
cases <- if (group == "observed") c("counts_marks", "flocking_original", "flocking_scalar") else
  c("univariate", "marked", "multi_kernel", "fitted_poisson")
rows <- checks <- list()
variant_values <- list()
for (n in c(1000L, 10000L)) for (case in cases) {
  d <- if (case %in% c("univariate", "fitted_poisson")) 1L else 2L
  k <- if (case == "multi_kernel") 4L else d
  set.seed(7271L + d)
  dt <- c(0, rexp(n, d * .5))
  types <- c(0L, sample.int(d, n, replace = TRUE))
  marks <- c(0, sample(1:3, n, replace = TRUE))
  h <- new("hspec", mu = wrap(matrix(.3, d)),
    alpha = wrap(matrix(.15 / k, d, k)),
    beta = wrap(matrix(rep(seq(1, 2, length.out = k), each = d), d)),
    eta = if (case == "marked") wrap(matrix(.01, d, k)) else NULL,
    type_col_map = if (k == d) NULL else list(c(3L, 1L), c(4L, 2L)))
  initial <- .1
  if (case == "counts_marks") {
    h@eta <- matrix(.01, d, k)
    h@impact <- function(param = c(extra = .02), n, N, Nc, mark, ...) {
      matrix(param * (sum(N[n, ]) %% 3 + sum(Nc[n, ]) %% 4 + mark[n]) / 10, 2, 2)
    }
  }
  if (case %in% c("flocking_original", "flocking_scalar")) {
    h <- flocking_interface(flocking, case == "flocking_scalar")
    d <- 4L; k <- 4L
    dt <- flocking_data$inter_arrival[seq_len(n + 1L)]
    types <- flocking_data$type[seq_len(n + 1L)]
    marks <- flocking_data$mark[seq_len(n + 1L)]
    initial <- NULL
  }
  if (case == "fitted_poisson") {
    poisson <- new("hspec", mu = function(param = c(log_rate = log(.3))) matrix(exp(param)),
      alpha = function(param = numeric(0)) matrix(0), beta = function(param = numeric(0)) matrix(1))
    options(emhawkes.engine = "auto")
    fit <- suppressMessages(hfit(poisson, dt, lambda_component0 = 0, method = "NR"))
    h <- fitted_hspec(fit); initial <- fit_initial_state(fit)$lambda_component0
  }
  for (infer in c(FALSE, TRUE)) {
    invoke <- function(engine) {
      options(emhawkes.engine = engine)
      suppressMessages(logLik(h, dt, type = types, mark = marks,
        lambda_component0 = initial, infer = infer))
    }
    reference <- invoke("r")
    if (case == "flocking_original") variant_values[[paste(n, infer)]] <- reference
    if (case == "flocking_scalar") stopifnot(identical(reference, variant_values[[paste(n, infer)]]))
    for (engine in c("auto", "cpp")) {
      native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
      calls <- 0L
      assignInNamespace("hawkes_loglik_cpp", function(...) { calls <<- calls + 1L; native(...) }, "emhawkes")
      candidate <- invoke(engine)
      assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
      stopifnot(calls == 1L, isTRUE(all.equal(candidate, reference, tolerance = 1e-10)))
      checks[[length(checks) + 1L]] <- data.frame(case, events = n, infer, engine,
        native_calls = calls, exact = identical(candidate, reference),
        max_abs_difference = max(abs(unlist(candidate, use.names = FALSE) - unlist(reference, use.names = FALSE))))
    }
    invisible(invoke("r")); invisible(invoke("auto"))
    batch <- if (group == "observed") { if (n == 1000L) 10L else 3L } else
      if (n == 1000L) 50L else 10L
    for (round in 1:3) for (engine in if (round %% 2L) c("r", "auto") else c("auto", "r")) {
      gc()
      time <- system.time(for (iteration in seq_len(batch)) value <- invoke(engine))
      stopifnot(isTRUE(all.equal(value, reference, tolerance = 1e-10)))
      rows[[length(rows) + 1L]] <- data.frame(case, events = n, dimensions = d, kernels = k,
        infer, engine, round, batch, seconds = unname(time[["elapsed"]]) / batch,
        cpu_seconds = unname(time[["user.self"]] + time[["sys.self"]]) / batch, affinity)
    }
  }
}
timings <- do.call(rbind, rows)
summary <- aggregate(cbind(seconds, cpu_seconds) ~ case + events + dimensions + kernels + infer + engine,
  timings, median)
write.csv(timings, file.path(out, "timings.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(out, "accuracy.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "session-info.txt"))
runtime_md5 <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!startsWith(lines, "#'")]
  temporary <- tempfile()
  on.exit(unlink(temporary))
  writeBin(charToRaw(enc2utf8(paste(lines, collapse = "\n"))), temporary)
  unname(tools::md5sum(temporary))
}
write.csv(data.frame(file = c("R/hengine.R", "R/hllf.R"),
  runtime_md5 = vapply(c("R/hengine.R", "R/hllf.R"), runtime_md5, "")),
  file.path(out, "source-hashes.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
