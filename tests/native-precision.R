library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")

# Exercise precision selection on the very first native call without a seed.
# The compiler may support long double while the running R intentionally does not.
had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
saved <- if (had_seed) .Random.seed else NULL
cancel <- new("hspec", mu = .3, alpha = matrix(0, 1, 3),
  beta = matrix(1, 1, 3), type_col_map = list(1:3))
initial <- matrix(c(1e16, 1, -1e16), 1)
one <- hsim(cancel, 1, lambda_component0 = initial)
reference <- r_sim(cancel, 1, lambda_component0 = initial)
stopifnot(identical(one, reference),
  identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), had_seed))
if (had_seed) stopifnot(identical(.Random.seed, saved))
expected <- .3 + rowSums(initial)
stopifnot(identical(unname(one$lambda[1, ]), unname(expected)))
old <- options(emhawkes.engine = "cpp")
replay <- logLik(cancel, 0, infer = TRUE, lambda_component0 = initial)
options(old)
stopifnot(identical(one$lambda, replay$lambda),
  identical(one$rambda, replay$rambda))

# This small mapped case was the first mismatch in the unmodified noLD build.
h <- new("hspec", mu = c(.3, .7), alpha = matrix(.015, 2, 4),
  beta = matrix(c(.8, 1.9), 2, 4), type_col_map = list(c(2L, 4L), c(3L, 4L)))
models <- list(h)
marked <- h
marked@eta <- matrix(.002, 2, 4)
marked@rmark <- function() runif(1, .5, 1.5)
marked@impact <- function(lambda_component_n) .01 * lambda_component_n
models[[2]] <- marked
for (model in models) for (seed in c(1, 62)) {
  set.seed(seed)
  reference <- r_sim(model, 31, lambda_component0 = .0075, N0 = 1:2, Nc0 = 3:4)
  rng <- .Random.seed; following <- rnorm(3)
  set.seed(seed)
  native <- hsim(model, 31, lambda_component0 = .0075, N0 = 1:2, Nc0 = 3:4)
  stopifnot(identical(reference, native), identical(rng, .Random.seed),
    identical(following, rnorm(3)))
  old <- options(emhawkes.engine = "cpp")
  replay <- logLik(model, native$inter_arrival, native$type, native$mark,
    infer = TRUE, lambda_component0 = .0075, N0 = 1:2, Nc0 = 3:4)
  options(old)
  for (field in c("lambda", "rambda", "lambda_component", "rambda_component")) {
    if (!identical(native[[field]], replay[[field]]))
      stop("Native precision mismatch: ", field, ", seed=", seed,
        ", marked=", !is.null(model@eta), ": ",
        paste(all.equal(native[[field]], replay[[field]], tolerance = 0), collapse = "; "))
  }
}
cat("R accumulation precision, initial-only RNG, mapped paths and callback histories passed.\n")
