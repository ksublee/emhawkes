library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")

# Wide mapped kernels expose accidental cross-row reuse or changed component
# summation. Compare the complete public result and RNG with the R reference.
for (d in c(1L, 4L)) for (k in c(16L, 32L, 128L)) {
  mapping <- lapply(seq_len(d), function(i) rev(seq(i, k, by = d)))
  h <- new("hspec", mu = seq(.3, .7, length.out = d),
    alpha = matrix(.35 / k, d, k),
    beta = matrix(rep(seq(.8, 1.9, length.out = d), k), d, k),
    type_col_map = mapping)
  initial <- matrix(seq(.01, .03, length.out = d * k) / k, d, k)
  for (callback in c(FALSE, TRUE)) {
    model <- h
    if (callback) {
      model@impact <- function(lambda_component_n) .006 * lambda_component_n
      model@rmark <- function() runif(1, .5, 1.5)
      model@eta <- matrix(.001 / k, d, k)
    }
    for (size in c(1L, 2L, 251L)) for (seed in c(62L, 7361L)) {
      set.seed(seed)
      expected <- r_sim(model, size, lambda_component0 = initial, N0 = 2, Nc0 = 5)
      rng <- .Random.seed
      set.seed(seed)
      actual <- hsim(model, size, lambda_component0 = initial, N0 = 2, Nc0 = 5)
      stopifnot(identical(expected, actual), identical(rng, .Random.seed))
    }
    # Independently reconstruct every component before and after each event.
    for (engine in c("r", "cpp")) {
      old <- options(emhawkes.engine = engine)
      replay <- logLik(model, actual$inter_arrival, actual$type, actual$mark,
        lambda_component0 = initial, N0 = 2, Nc0 = 5, infer = TRUE)
      options(old)
      for (field in c("lambda", "rambda", "lambda_component", "rambda_component"))
        stopifnot(identical(actual[[field]], replay[[field]]))
      stopifnot(identical(rng, .Random.seed))
    }
  }
}
