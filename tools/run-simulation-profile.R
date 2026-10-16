# Rscript --vanilla tools/run-simulation-profile.R plan|run LIB PROFILE_LIB OUT HSIM_SOURCE [stress]
# Run sequentially from the repository root. PROFILE_LIB is installed with
# R CMD INSTALL --with-keep.source. The optional stress run takes about 100 s.
a <- commandArgs(TRUE); action <- a[1]; lib <- a[2]; profile_lib <- a[3]; out <- a[4]; source <- a[5]
stopifnot(action %in% c("plan", "run"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cases <- c("univariate", "multivariate", "multi_equal", "multi_unequal", "large_equal", "large_unequal",
  "marked", "parameter_mu", "state_current", "state_dots", "state_observed", "signed")
plan <- list(); add <- function(id, events, mode, round = 1L, variant = "reference") {
  plan[[length(plan) + 1L]] <<- data.frame(id, events, mode, round, variant)
}
for (n in c(1000L, 10000L)) for (round in 1:3) for (id in cases)
  add(id, if (n == 10000L && id %in% c("large_equal", "large_unequal")) 3000L else n, "time", round)
for (id in c("univariate", "multi_equal", "multi_unequal", "large_equal", "state_current", "state_dots", "state_observed"))
  for (n in if (id == "large_equal") c(1000L, 3000L) else c(1000L, 3000L, 10000L)) add(id, n, "memory")
for (id in c("multi_equal", "large_equal", "state_current", "state_dots", "state_observed")) {
  n <- if (id == "large_equal") 3000L else 10000L
  for (round in 1:3) for (variant in if (round %% 2L) c("arrival_row", "selective_context", "both") else c("both", "selective_context", "arrival_row"))
    add(id, n, "time", round, variant)
  for (variant in c("arrival_row", "selective_context", "both")) add(id, n, "memory", 1L, variant)
}
for (id in c("univariate", "multi_equal", "multi_unequal", "large_equal", "state_current", "state_dots", "state_observed", "signed"))
  add(id, if (id == "large_equal") 3000L else 10000L, "profile")
for (id in cases) for (variant in c("reference", "both")) add(id, 1000L, "probe", 1L, variant)
for (id in c("multi_equal", "multi_unequal", "large_equal", "state_current"))
  add(id, if (id == "large_equal") 3000L else 10000L, "profile", 1L, "both")
for (round in 1:3) add("large_equal", 10000L, "time", round, "both")
for (id in c("large_equal", "state_current")) {
  for (round in 1:3) add(id, 100000L, "time", round, "both")
  add(id, 100000L, "peak", 1L, "both")
}
for (id in c("large_equal", "state_current")) add(id, 100000L, "profile", 1L, "both")
if (length(a) >= 6L && a[6] == "stress") add("large_equal", 10000L, "time")
plan <- do.call(rbind, plan); write.csv(plan, file.path(out, "design.csv"), row.names = FALSE)
if (action == "plan") {print(plan, row.names = FALSE); quit()}
# Pin installed code across resumes even when the package version is unchanged.
installed <- unlist(lapply(c(lib, profile_lib), function(path) {
  root <- file.path(path, "emhawkes")
  c(file.path(root, "DESCRIPTION"), list.files(file.path(root, "R"), pattern = "[.]rd[bx]$", full.names = TRUE),
    list.files(file.path(root, "libs"), pattern = "[.](dll|so|dylib)$", full.names = TRUE, recursive = TRUE))
}), use.names = FALSE)
stopifnot(length(installed) >= 6L, all(file.exists(installed)))
code_hashes <- data.frame(path = normalizePath(installed, winslash = "/"), md5 = unname(tools::md5sum(installed)))
record <- file.path(out, "installed-code.csv")
if (file.exists(record)) {
  stopifnot(identical(code_hashes, read.csv(record, stringsAsFactors = FALSE)))
} else write.csv(code_hashes, record, row.names = FALSE)
# Immutable workers avoid changing a running script when editing tools.
stage <- file.path(out, "worker-source"); dir.create(stage, showWarnings = FALSE)
inputs <- c("tools/profile-simulation.R", "tools/simulation-profile-cases.R", source)
for (input in inputs) {
  target <- file.path(stage, basename(input))
  if (file.exists(target)) stopifnot(unname(tools::md5sum(input)) == unname(tools::md5sum(target)))
  else stopifnot(file.copy(input, target))
}
rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (i in seq_len(nrow(plan))) {
  z <- plan[i, ]; file <- paste(z$id, z$events, z$variant, z$mode, z$round, "rds", sep = ".")
  if (file.exists(file.path(out, file))) next
  args <- c("--vanilla", file.path(stage, "profile-simulation.R"),
    if (z$mode == "profile") profile_lib else lib, out, z$id, z$events, z$mode, z$round,
    file.path(stage, "simulation-profile-cases.R"), z$variant, file.path(stage, basename(source)))
  status <- system2(rs, vapply(args, shQuote, "")); stopifnot(status == 0L)
}
