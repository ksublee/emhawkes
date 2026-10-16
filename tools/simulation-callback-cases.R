source("tools/simulation-profile-cases.R")
.callback_base_case <- simulation_profile_case
simulation_profile_case <- function(id) {
  if (id %in% c("parameter_mu", "state_current", "marked", "state_dots", "state_observed",
                "multi_equal", "multi_unequal")) return(.callback_base_case(id))
  x <- .callback_base_case("multi_equal")
  if (id == "mark_only") x$object@rmark <- function() stats::rgeom(1, .7) + 1
  if (id == "eta_only") x$object@eta <- matrix(.015 / 8, 4, 8)
  if (id == "heavy_impact") x$object@impact <- function(lambda_component_n) {
    work <- crossprod(matrix(sin(seq_len(1024) + sum(lambda_component_n)), 32))
    .006 * lambda_component_n + sum(diag(work)) * 1e-9
  }
  stopifnot(id %in% c("mark_only", "eta_only", "heavy_impact"))
  x$meta$id <- id
  x
}
