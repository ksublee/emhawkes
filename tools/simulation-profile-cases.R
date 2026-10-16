# Fixed Markov fixtures for hsim profiling; no likelihood estimation is timed.
simulation_profile_case <- function(id) {
  design <- data.frame(id = c("univariate", "multivariate", "multi_equal", "multi_unequal",
    "large_equal", "large_unequal", "marked", "parameter_mu", "state_current",
    "state_dots", "state_observed", "signed"),
    d = c(1,4,4,4,8,8,4,4,4,4,4,2), k = c(1,4,8,8,16,16,8,8,8,8,8,4),
    stringsAsFactors = FALSE)
  z <- design[match(id, design$id), ]; stopifnot(nrow(z) == 1L, !is.na(z$id))
  d <- z$d; k <- z$k; rates <- seq(1.1, 1.5, length.out = d)
  beta <- matrix(rep(rates, k), d, k)
  unequal <- id %in% c("multi_unequal", "large_unequal")
  if (unequal) beta <- outer(rates, seq(.7, 2.2, length.out = k))
  alpha <- matrix(.35 / k, d, k); mu <- matrix(.5, d)
  mapping <- if (k == d) NULL else lapply(seq_len(d), function(i) seq(i, k, by = d))
  initial <- matrix(.03 / k, d, k)
  eta <- impact <- rmark <- NULL
  if (id == "parameter_mu") mu <- function(param = c(level = .5)) matrix(rep(param[1], d), d)
  if (id %in% c("marked", "state_observed")) rmark <- function() stats::rgeom(1, .7) + 1
  if (id == "marked") eta <- matrix(.015 / k, d, k)
  if (id == "state_current") impact <- function(lambda_component_n) .06 * lambda_component_n / (1 + lambda_component_n)
  if (id == "state_dots") impact <- function(lambda_component_n, ...) .06 * lambda_component_n / (1 + lambda_component_n)
  if (id == "state_observed") {
    mu <- function(n, N) matrix(.5 * (1 + .08 * (N[max(1L, n - 1L), ] %% 3)), d)
    impact <- function(n, mark, lambda_component_n) {
      .06 * (1 + .1 * (mark[n] - 1)) * lambda_component_n / (1 + lambda_component_n)
    }
  }
  if (id == "signed") {
    alpha <- matrix(c(.08, -.04, .05, -.02, .05, -.02, .08, -.04), 2, byrow = TRUE)
    mapping <- list(c(2L, 1L), c(4L, 3L)); initial <- alpha / 2
  }
  object <- new("hspec", mu = mu, alpha = alpha, beta = beta, eta = eta,
    type_col_map = mapping, impact = impact, rmark = rmark)
  list(object = object, initial = initial, meta = transform(z, arrival = if (unequal) "thinning" else "direct"),
    seed = if (id %in% c("multi_equal", "parameter_mu", "state_current", "state_dots")) 7361L else
      7361L + match(id, design$id))
}

# Diagnostic-only method copies. These are not installed or selected by hsim().
# The arrival probe changes only the private Hawkes helper's row index/input;
# every model in this study has rresidual = NULL. Public callbacks keep n.
simulation_profile_variant <- function(variant, source) {
  if (variant == "reference") return(emhawkes::hsim)
  stopifnot(variant %in% c("method_reference", "arrival_row", "selective_context", "both"))
  text <- paste(readLines(source, warn = FALSE), collapse = "\n")
  replace_one <- function(old, new) {
    stopifnot(lengths(regmatches(text, gregexpr(old, text, fixed = TRUE))) == 1L)
    text <<- sub(old, new, text, fixed = TRUE)
  }
  if (variant %in% c("arrival_row", "both")) replace_one(
    "res <- handle_arrival(n = n, rambda_component = rambda_component,",
    "res <- handle_arrival(n = 2L, rambda_component = rambda_component[n - 1L, , drop = FALSE],")
  if (variant %in% c("selective_context", "both")) {
    old <- 'context_at <- function(n) list(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                                   N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                                   lambda_component_n = lambda_component_n, mu = mu, alpha = alpha, beta = beta)'
    new <- 'context_at <- function(n, plan) {
      scope <- environment()
      lapply(alist(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
        N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
        lambda_component_n = lambda_component_n, mu = mu, alpha = alpha, beta = beta)[plan],
        eval, envir = scope)
    }'
    replace_one(old, new)
    for (index in c("1", "2", "n + 1")) replace_one(
      paste0("call_component(mu, callbacks$mu, context_at(", index, "))"),
      paste0("call_component(mu, callbacks$mu, context_at(", index, ", callbacks$mu))"))
    for (name in c("rmark", "impact")) replace_one(
      paste0("call_component(object@", name, ", callbacks$", name, ", context_at(n))"),
      paste0("call_component(object@", name, ", callbacks$", name, ", context_at(n, callbacks$", name, "))"))
  }
  expressions <- parse(text = text, keep.source = FALSE)
  methods <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("setMethod")), as.list(expressions))
  stopifnot(length(methods) == 1L)
  eval(methods[[1]][["definition"]], envir = asNamespace("emhawkes"))
}
