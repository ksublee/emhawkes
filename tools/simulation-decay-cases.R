source("tools/simulation-callback-cases.R")
.decay_previous_case <- simulation_profile_case
simulation_profile_case <- function(id) {
  if (id %in% c("univariate", "large_equal", "signed")) return(.callback_base_case(id))
  if (!id %in% c("wide32", "wide128", "impact32")) return(.decay_previous_case(id))
  x <- .callback_base_case("multi_equal")
  d <- 4L; k <- if (id == "wide128") 128L else 32L
  x$object <- new("hspec", mu = rep(.5, d), alpha = matrix(.35 / k, d, k),
    beta = matrix(rep(seq(1.1, 1.5, length.out = d), k), d),
    type_col_map = lapply(seq_len(d), function(i) seq(i, k, by = d)))
  x$initial <- matrix(.03 / k, d, k)
  if (id == "impact32") x$object@impact <- function(lambda_component_n) .006 * lambda_component_n
  x$meta$d <- d; x$meta$k <- k
  x
}
