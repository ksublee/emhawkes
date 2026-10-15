# Rscript --vanilla tools/benchmark-state-replay.R LIB OUT
# Public standalone logLik: fixed observations and independent callback defaults.
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(args[1])))
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
affinity <- Sys.getenv("EMHAWKES_BENCH_AFFINITY", "scheduler")
if (affinity != "scheduler") {
  stopifnot(.Platform$OS.type == "windows", grepl("^[1-9][0-9]*$", affinity))
  actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command",
    shQuote(sprintf("$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), affinity))), stdout = TRUE)
  stopifnot(length(actual) == 1L, trimws(actual) == affinity)
}
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(rscript, c("--vanilla", "tools/benchmark-bid-ask-numeric.R", "prepare",
  shQuote(args[1]), shQuote(file.path(out, "prepare"))))
stopifnot(status == 0L)
fixtures <- readRDS(file.path(out, "prepare", "fixtures.rds"))
set.seed(63001)
sim <- hsim(fixtures$bid_ask_current$arguments$object, size = 10001L,
  lambda_component0 = fixtures$bid_ask_current$arguments$lambda_component0)
rows <- checks <- memory <- list()
variants <- list()
for (n in c(1000L, 10000L)) for (case in c("observed_baseline", names(fixtures))) {
  if (case == "observed_baseline") {
    h <- new("hspec", mu = function(param = c(base = .3), n, Nc, ...) {
      matrix(param + .01 * (Nc[max(1, n - 1), ] %% 5), 2)
    }, alpha = matrix(.04, 2, 2), beta = matrix(c(1, 1.2, 1.5, 2), 2))
    set.seed(63002)
    dt <- c(0, rexp(n, 1)); types <- c(0, sample.int(2, n, TRUE))
    marks <- c(0, sample.int(3, n, TRUE)); initial <- .1
  } else {
    h <- fixtures[[case]]$arguments$object
    dt <- sim$inter_arrival[seq_len(n + 1L)]
    types <- sim$type[seq_len(n + 1L)]; marks <- sim$mark[seq_len(n + 1L)]
    initial <- fixtures[[case]]$arguments$lambda_component0
  }
  for (infer in c(FALSE, TRUE)) {
    invoke <- function(engine) {
      options(emhawkes.engine = engine)
      suppressMessages(logLik(h, dt, type = types, mark = marks,
        lambda_component0 = initial, infer = infer))
    }
    reference <- invoke("r")
    if (case == "bid_ask") variants[[paste(n, infer)]] <- reference
    if (case == "bid_ask_current") stopifnot(identical(reference, variants[[paste(n, infer)]]))
    for (engine in c("auto", "cpp")) {
      native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes"); calls <- 0L
      assignInNamespace("hawkes_loglik_cpp", function(...) { calls <<- calls + 1L; native(...) }, "emhawkes")
      candidate <- invoke(engine)
      assignInNamespace("hawkes_loglik_cpp", native, "emhawkes")
      stopifnot(calls == 1L, isTRUE(all.equal(candidate, reference, tolerance = 1e-10)))
      checks[[length(checks) + 1L]] <- data.frame(case, events=n, infer, engine,
        native_calls=calls, exact=identical(candidate, reference),
        max_abs_difference=max(abs(unlist(candidate, use.names=FALSE) - unlist(reference, use.names=FALSE))))
    }
    invisible(invoke("r")); invisible(invoke("auto"))
    batch <- if (n == 1000L) 5L else 2L
    for (round in 1:3) for (engine in if (round %% 2L) c("r", "auto") else c("auto", "r")) {
      gc()
      time <- system.time(for (iteration in seq_len(batch)) value <- invoke(engine))
      stopifnot(identical(value, reference))
      rows[[length(rows) + 1L]] <- data.frame(case, events=n, infer, engine, round, batch,
        seconds=unname(time[["elapsed"]])/batch,
        cpu_seconds=unname(time[["user.self"]]+time[["sys.self"]])/batch, affinity)
    }
    # One separate profiled call; parsing is excluded from both timing and R heap.
    for (engine in c("r", "auto")) {
      profile <- tempfile(); gc(reset=TRUE)
      Rprofmem(profile)
      value <- tryCatch(invoke(engine), finally=Rprofmem(NULL))
      heap <- gc()
      stopifnot(identical(value, reference))
      con <- file(profile, "r"); bytes <- allocations <- 0
      repeat {
        lines <- readLines(con, n=10000L, warn=FALSE)
        if (!length(lines)) break
        sizes <- suppressWarnings(as.numeric(sub(" .*", "", substr(lines, 1, 24))))
        bytes <- bytes + sum(sizes, na.rm=TRUE)
        allocations <- allocations + sum(!is.na(sizes))
      }
      close(con); unlink(profile)
      memory[[length(memory)+1L]] <- data.frame(case, events=n, infer, engine,
        allocated_bytes=bytes, allocations, max_R_heap_MiB=sum(heap[,6]))
    }
    cat(case, n, "infer", infer, "verified\n"); flush.console()
  }
}
timings <- do.call(rbind, rows)
summary <- aggregate(cbind(seconds, cpu_seconds) ~ case + events + infer + engine,
  timings, median)
write.csv(do.call(rbind, memory), file.path(out, "memory.csv"), row.names=FALSE)
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
write.csv(data.frame(file = c("R/hengine.R", "R/hllf.R", "src/hawkes_loglik.cpp"),
  runtime_md5 = vapply(c("R/hengine.R", "R/hllf.R", "src/hawkes_loglik.cpp"), runtime_md5, "")),
  file.path(out, "source-hashes.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
