# Rscript --vanilla tools/profile-simulation.R LIB OUT ID EVENTS MODE ROUND CASES
# MODE: time, memory, peak, profile, or probe. Optional args: VARIANT HSIM_SOURCE.
# VARIANT defaults to reference. CASES is simulation-profile-cases.R.
a <- commandArgs(TRUE)
lib <- a[1]; out <- a[2]; id <- a[3]; events <- as.integer(a[4]); mode <- a[5]; round <- as.integer(a[6])
.libPaths(c(lib, .libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))), normalizePath(lib)))
source(a[7]); x <- simulation_profile_case(id)
variant <- if (length(a) >= 8L) a[8] else "reference"
simulator <- simulation_profile_variant(variant, if (length(a) >= 9L) a[9] else "R/hsim.R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
key <- paste(id, events, variant, mode, round, sep = ".")
if (.Platform$OS.type == "windows" && nzchar(Sys.getenv("EMHAWKES_BENCH_AFFINITY"))) {
  mask <- Sys.getenv("EMHAWKES_BENCH_AFFINITY"); stopifnot(grepl("^[1-9][0-9]*$", mask))
  actual <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf(
    "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()", Sys.getpid(), mask))), stdout = TRUE)
  stopifnot(trimws(actual) == mask)
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
run <- function(n = events) simulator(x$object, size = n + 1L, lambda_component0 = x$initial)
for (i in 1:3) {set.seed(x$seed); invisible(run(min(events, 100L)))}
rss <- function() if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
  c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))), stdout = TRUE)) else NA_real_
gc(); heap_before <- gc(reset = TRUE); rss_before <- if (mode %in% c("memory", "peak")) rss() else NA_real_
set.seed(x$seed); seed_before <- .Random.seed
warnings <- character(); invoke <- function() withCallingHandlers(run(), warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
})
profile_file <- file.path(out, paste0(key, ".Rprof"))
allocation_file <- file.path(out, paste0(key, ".Rprofmem"))
counts <- memory <- profile <- NULL; elapsed <- cpu <- NA_real_
if (mode == "time") {
  timing <- system.time(value <- invoke()); elapsed <- unname(timing["elapsed"])
  cpu <- unname(timing["user.self"] + timing["sys.self"])
} else if (mode == "peak") {
  value <- invoke()
} else if (mode == "memory") {
  Rprofmem(allocation_file); value <- tryCatch(invoke(), finally = Rprofmem(NULL))
} else if (mode == "profile") {
  Rprof(profile_file, interval = .001, line.profiling = TRUE, gc.profiling = TRUE,
    numfiles = 1000L, bufsize = 1000000L)
  value <- tryCatch(invoke(), finally = Rprof(NULL))
} else if (mode == "probe") {
  .sim_probe <- new.env(parent = emptyenv())
  for (name in c("exp_calls", "exp_draws", "uniform_calls", "uniform_draws", "geom_calls", "callbacks", "arrivals")) .sim_probe[[name]] <- 0
  probe_run <- function() {
    for (name in c("rexp", "runif", "rgeom")) {
      tracer <- switch(name,
        rexp = quote({p <- get(".sim_probe", envir = .GlobalEnv); p$exp_calls <- p$exp_calls + 1; p$exp_draws <- p$exp_draws + n}),
        runif = quote({p <- get(".sim_probe", envir = .GlobalEnv); p$uniform_calls <- p$uniform_calls + 1; p$uniform_draws <- p$uniform_draws + n}),
        rgeom = quote({p <- get(".sim_probe", envir = .GlobalEnv); p$geom_calls <- p$geom_calls + 1}))
      invisible(capture.output(trace(name, where = asNamespace("stats"), tracer = tracer, print = FALSE)))
    }
    on.exit(for (name in c("rexp", "runif", "rgeom")) suppressMessages(untrace(name, where = asNamespace("stats"))), add = TRUE)
    invisible(capture.output(trace("call_component", where = asNamespace("emhawkes"),
      tracer = quote({p <- get(".sim_probe", envir = .GlobalEnv); p$callbacks <- p$callbacks + 1}), print = FALSE)))
    invisible(capture.output(trace("rarrival", where = asNamespace("emhawkes"),
      tracer = quote({p <- get(".sim_probe", envir = .GlobalEnv); p$arrivals <- p$arrivals + 1}), print = FALSE)))
    on.exit({suppressMessages(untrace("call_component", where = asNamespace("emhawkes"))); suppressMessages(untrace("rarrival", where = asNamespace("emhawkes")))}, add = TRUE)
    invoke()
  }
  value <- suppressMessages(probe_run()); counts <- as.list(.sim_probe)
} else stop("Unknown profiling mode.")
seed_after <- .Random.seed; heap_after <- gc(); rss_after <- if (mode %in% c("memory", "peak")) rss() else NA_real_
stopifnot(length(value$inter_arrival) == events + 1L, length(warnings) == 0L,
  all(is.finite(value$inter_arrival)), all(value$inter_arrival >= 0),
  all(value$type[-1] >= 1 & value$type[-1] <= x$meta$d))
hash <- function(object) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(object, f, compress = FALSE, version = 3)
  unname(tools::md5sum(f))
}
fingerprint <- hash(value[setdiff(names(value), "hspec")]); rng <- hash(seed_after)
if (mode == "memory") {
  con <- file(allocation_file, "r"); groups <- new.env(hash = TRUE, parent = emptyenv())
  total <- largest <- blocks <- 0
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE); if (!length(lines)) break
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines))); good <- is.finite(bytes)
    total <- total + sum(bytes[good]); blocks <- blocks + sum(good); if (any(good)) largest <- max(largest, bytes[good])
    tab <- table(lines[good]); for (key0 in names(tab)) {
      old <- groups[[key0]]; groups[[key0]] <- if (is.null(old)) unname(tab[[key0]]) else old + unname(tab[[key0]])
    }
  }
  close(con)
  keys <- ls(groups); frequencies <- vapply(keys, function(k) groups[[k]], 0)
  sizes <- as.numeric(sub(" .*", "", keys))
  allocations <- data.frame(bytes = sizes, count = frequencies, total_bytes = sizes * frequencies,
    stack = sub("^[0-9]+ *: *", "", keys), row.names = NULL)
  allocations <- allocations[order(-allocations$total_bytes), ]
  write.csv(allocations, file.path(out, paste0(paste(id, events, variant, sep = "."), ".allocations.csv")), row.names = FALSE)
  widths <- unique(c(1L, x$meta$d, x$meta$d * x$meta$k))
  shapes <- data.frame(width = widths, block_bytes = 8 * (events + 1) * widths + 48)
  shapes$count <- vapply(shapes$block_bytes, function(b) sum(frequencies[sizes == b]), 0)
  shapes$total_bytes <- shapes$block_bytes * shapes$count
  memory <- list(allocated_bytes = total, allocation_count = blocks, largest_block = largest, shapes = shapes)
}
if (mode == "profile") {
  profile <- summaryRprof(profile_file, lines = "both")
  write.csv(profile$by.self, file.path(out, paste0(id, ".", variant, ".self.csv")))
  write.csv(profile$by.total, file.path(out, paste0(id, ".", variant, ".total.csv")))
  if (!is.null(profile$by.line)) write.csv(profile$by.line, file.path(out, paste0(id, ".", variant, ".lines.csv")))
}
saveRDS(list(id = id, meta = x$meta, events = events, mode = mode, round = round, variant = variant,
  seconds = elapsed, cpu_seconds = cpu, memory = memory, profile = profile,
  heap_before_mib = sum(heap_before[, 2]), heap_peak_mib = sum(heap_after[, 6]),
  rss_before = rss_before, rss_peak = rss_after, returned_bytes = as.numeric(object.size(value)),
  fingerprint = fingerprint, rng = rng, seed = x$seed, counts = counts, warnings = warnings),
  file.path(out, paste0(key, ".rds")))
cat(id, events, variant, mode, round, "seconds", elapsed, "returned_MiB", as.numeric(object.size(value))/2^20, "\n")
