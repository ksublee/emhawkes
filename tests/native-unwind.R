library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")

# R conditions raised inside a native callback must preserve the consumed RNG
# state, including when a subsequent native simulation enters a fresh scope.
exercise_unwind <- function(native, component, kind) {
  calls <- 0L
  fail <- function() {
    runif(1)
    calls <<- calls + 1L
    if (calls == 3L) {
      if (kind == "interrupt") stop(structure(list(message = "callback interrupted", call = NULL),
        class = c("interrupt", "condition")))
      stop("callback failed", call. = FALSE)
    }
  }
  h <- new("hspec", mu = .4, alpha = .03, beta = 1)
  if (component == "mark") h@rmark <- function() { fail(); 1 } else
    h@impact <- function(lambda_component_n) { fail(); .01 * lambda_component_n }
  simulation <- if (native) hsim else r_sim
  set.seed(7361)
  result <- tryCatch(simulation(h, 9, lambda_component0 = .02),
    error = function(e) e, interrupt = function(e) e)
  stopifnot(inherits(result, "condition"), calls == 3L)
  after_failure <- .Random.seed
  draws <- runif(3)
  recovery <- simulation(new("hspec", mu = .3, alpha = 0, beta = 1), 5,
                         lambda_component0 = 0)
  list(kind = if (inherits(result, "interrupt")) "interrupt" else "error",
    message = conditionMessage(result), calls = calls, after_failure = after_failure,
    draws = draws, recovery = recovery, final_rng = .Random.seed)
}
for (component in c("mark", "impact")) for (kind in c("error", "interrupt")) {
  expected <- exercise_unwind(FALSE, component, kind)
  actual <- exercise_unwind(TRUE, component, kind)
  stopifnot(identical(expected$kind, kind), identical(expected, actual))
}
