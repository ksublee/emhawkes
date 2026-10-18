library(emhawkes)
reference_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")
# With no excitation jumps and beta*t below machine precision, the current
# component stays exactly one. Independent draws therefore give the path
# directly, without evaluating the simulator's arrival inversion formula.
for (beta in c(5e-324, 1e-320, 1e-308)) {
  h <- new("hspec", mu = 1e-20, alpha = 0, beta = beta)
  set.seed(4001)
  waiting <- vapply(seq_len(1000), function(i) min(rexp(1, 1e-20), -log(runif(1))), numeric(1))
  final_seed <- .Random.seed
  paths <- list()
  for (simulate in list(hsim, reference_sim)) {
    set.seed(4001)
    path <- simulate(h, size = 1001, lambda_component0 = 1)
    stopifnot(all(abs(path$inter_arrival[-1] - waiting) <=
        4 * .Machine$double.eps * pmax(1, waiting)),
      all(path$inter_arrival[-1] > 0), identical(.Random.seed, final_seed))
    paths[[length(paths) + 1L]] <- path
  }
  stopifnot(identical(paths[[1]], paths[[2]]))
}
# The native callback loop uses the same stable arrival calculation and keeps
# mark draws in the same position in R's RNG stream.
h@rmark <- function() 1 + runif(1)
set.seed(152)
native <- hsim(h, 30, lambda_component0 = 1)
seed <- .Random.seed
set.seed(152)
reference <- reference_sim(h, 30, lambda_component0 = 1)
stopifnot(identical(native, reference), identical(seed, .Random.seed))
