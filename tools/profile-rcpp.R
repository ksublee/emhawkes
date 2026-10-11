# Install the package version to measure, then run from the repository root:
# Rscript tools/profile-rcpp.R 10001
# An optional second argument selects the output directory.
library(emhawkes)

args <- commandArgs(trailingOnly = TRUE)
size <- if (length(args)) as.integer(args[1]) else 10001L
if (length(size) != 1L || is.na(size) || size < 2L) stop("Use at least two rows, including the initial row.")
output <- if (length(args) >= 2L) args[2] else tempfile("emhawkes-profile-")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)

dt <- c(0, rep(c(.1, .2, .3), length.out = size - 1L))
types <- c(0, rep(1:2, length.out = size - 1L))
marks <- c(0, rep(c(1, 2, 1), length.out = size - 1L))
constant <- new("hspec", mu = c(.2, .3), alpha = matrix(c(.2, .1, .1, .2), 2),
                beta = matrix(1.2, 2, 2))
multi <- new("hspec", mu = c(.2, .3), alpha = matrix(.08, 2, 4),
             beta = matrix(rep(c(.8, 1.2, 1.6, 2), each = 2), 2),
             type_col_map = list(1:2, 3:4))
state <- new("hspec",
  mu = function(param = c(base = .2, effect = .01), n, N) {
    counts <- if (n <= 1L) c(0, 0) else N[n - 1L, ]
    matrix(param[1] + param[2] * (counts %% 5), 2, 1)
  }, alpha = constant@alpha, beta = constant@beta,
  impact = function(param = c(effect = .02), mark, n, lambda_component_n) {
    matrix(param[1] * (mark[n] - 1), 2, 2)
  })
frpp <- set_residual(constant, param = c(a = .7, ell = 1.6),
                     dresidual = dtzexp, presidual = ptzexp,
                     qresidual = qtzexp, rresidual = rtzexp)
cases <- list(constant = constant, multi_kernel = multi, state = state, frpp = frpp)

measure <- function(object, name) {
  observation <- emhawkes:::normalize_observation(object, dt, types, marks)
  callbacks <- emhawkes:::model_callback_plans(object)
  evaluate <- function(infer = FALSE) {
    emhawkes:::loglik_hspec(object, observation, .1, infer, callbacks,
                           show_initial_message = FALSE)
  }
  evaluate()
  elapsed <- replicate(3, system.time(evaluate())[["elapsed"]])
  inferred <- replicate(3, system.time(evaluate(TRUE))[["elapsed"]])
  profile <- file.path(output, paste0(name, ".Rprof"))
  Rprof(profile, interval = .005, memory.profiling = TRUE)
  tryCatch(for (i in 1:12) evaluate(), finally = Rprof(NULL))
  summary <- summaryRprof(profile, memory = "both")
  capture.output(summary, file = file.path(output, paste0(name, "-summary.txt")))
  data.frame(case = name, rows = size, dimensions = object@dimens,
             kernels = ncol(emhawkes:::eval_param(object@beta)),
             median_loglik = median(elapsed), median_infer = median(inferred))
}

results <- do.call(rbind, lapply(names(cases), function(name) measure(cases[[name]], name)))
write.csv(results, file.path(output, "timings.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "session-info.txt"))
print(results, row.names = FALSE)
cat("Profiles and session information:", output, "\n")
