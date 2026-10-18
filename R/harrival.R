#' @include hspec.R
NULL

rarrival <- function(...){

  # Extract named arguments from ...
  args <- list(...)

  # Assign expected arguments to local variables
  rambda_component <- args$rambda_component
  mu <- args$mu
  beta <- args$beta
  dimens <- args$dimens
  type_col_map <- args$type_col_map
  n <- args$n
  current_component <- args$current_component
  if (is.null(current_component)) {
    current_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)
  }


  if (args$all_rows_beta_equal){

    #mu can be zero, so warning is turned off for a moment
    candidate_arrival <- suppressWarnings(stats::rexp(dimens, rate = mu))

    rowSums_rambda_component_n_miuns_1 <- rowSums(current_component)

    log_u <- log(stats::runif(dimens))
    product <- beta[, 1] * log_u
    scaled_log_u <- product / rowSums_rambda_component_n_miuns_1
    excited_arrival <- rep(Inf, dimens)
    finite_arrival <- rowSums_rambda_component_n_miuns_1 > 0 & scaled_log_u > -1
    excited_arrival[finite_arrival] <- -log1p(scaled_log_u[finite_arrival]) / beta[finite_arrival, 1]
    tiny <- which(rowSums_rambda_component_n_miuns_1 > 0 &
      (abs(product) < .Machine$double.xmin | abs(scaled_log_u) < .Machine$double.xmin))
    if (length(tiny)) {
      first_order <- -log_u[tiny] / rowSums_rambda_component_n_miuns_1[tiny]
      q <- beta[tiny, 1] * first_order
      correction <- rep(1, length(q))
      interior <- is.finite(q) & q > 0 & q < 1
      correction[interior] <- -log1p(-q[interior]) / q[interior]
      correction[!is.finite(q) | q >= 1] <- Inf
      excited_arrival[tiny] <- first_order * correction
    }
    candidate_arrival <- cbind(candidate_arrival, excited_arrival)

    inter_arrival <- min(candidate_arrival, na.rm = TRUE)

    # Rows identify event types; kernel column mappings do not change them.
    type <- if (is.finite(inter_arrival)) which(candidate_arrival == inter_arrival, arr.ind = TRUE)[1, 1] else NA_real_


  } else{

    # thinning algorithm

    rambda_component_n_miuns_1 <- current_component
    rowSums_rambda_component_n_miuns_1 <- rowSums(rambda_component_n_miuns_1)
    rambda_n_minus_1 <- rowSums_rambda_component_n_miuns_1 + mu

    candidate_arrival <- rep(0, dimens)

    # Loop through each component
    for(i in 1:dimens){

      s <- 0 # Initialize time s, which will represent the time to the next event.

      if(rambda_n_minus_1[i] == 0){

        candidate_arrival[i] <- Inf

      } else if (mu[i] == 0) {

        candidate_arrival[i] <- invert_exp_hazard(stats::rexp(1),
          rambda_component_n_miuns_1[i, ], 0, beta[i, ])

      } else {

        while(TRUE){

          # Sample the time until the next event from an exponential distribution
          # based on the current baseline + lambda intensity.
          s <- s + stats::rexp(1, rambda_n_minus_1[i])

          # Generate a uniform random number U.
          # U is used in the acceptance-rejection step of the thinning algorithm.
          U <- stats::runif(1, 0, rambda_n_minus_1[i])

          # Calculate the intensity at time s.
          # This is where the impact from past events (stored in current_lambda_component)
          # and the decay of this influence over time by beta are applied.
          f <- mu[i] + sum(rambda_component_n_miuns_1[i, ] * exp((- beta[i, ] * s)))
          if (is.na(U) || is.na(f)) {
            stop("Simulation encountered NaN during thinning; rescale time units to keep rates and waiting times representable.", call. = FALSE)
          }

          # If the random value U is less than (f + mu),
          # accept s as the next event time and exit the loop.
          if( U < f ) {
            candidate_arrival[i] <- s
            break
          }

        }
      }
    }


    inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    type <- if (is.finite(inter_arrival)) which.min(candidate_arrival) else NA_real_

  }

  c(inter_arrival = unname(inter_arrival), type = unname(type))
}

# Invert a nonnegative, decaying intensity. A finite total hazard leaves
# positive probability at an infinite waiting time.
invert_exp_hazard <- function(target, component, mu, beta) {
  if (length(target) != 1L || is.na(target) || target < 0) stop("Residual waiting times must be nonnegative.")
  if (is.infinite(target)) return(Inf)
  if (target == 0) return(0)
  if (mu == 0 && target >= sum(component / beta)) return(Inf)
  hazard <- function(time) mu * time + sum(integrate_exp_kernel(component, beta, time))
  upper <- target / (mu + sum(component))
  if (upper == 0) return(0)
  while (is.finite(upper) && hazard(upper) < target) upper <- upper * 2
  if (!is.finite(upper)) return(Inf)
  upper * stats::uniroot(function(x) hazard(upper * x) - target,
                        interval = c(0, 1), tol = 1e-12, check.conv = TRUE)$root
}


# for flexible model arrival

phi <- function(tau, rowSum_rambda_component_n_miuns_1, mu, beta){

  # rambda means after jump intensity
  # for arrival time n, row sum of rambda component at n-1 is needed

  mu * tau + integrate_exp_kernel(rowSum_rambda_component_n_miuns_1, beta, tau)

}

phi_wo_mu <- function(tau, rowSum_rambda_component_n_miuns_1,  beta){

  # rambda means after jump intensity
  # for arrival time n, row sum of rambda component at n-1 is needed

  integrate_exp_kernel(rowSum_rambda_component_n_miuns_1, beta, tau)

}


psi <- function(tau, rowSum_rambda_component_n_miuns_1, mu, beta){

  # rambda means after jump intensity
  # for arrival time n, row sum of rambda component at n-1 is needed

  mu +  rowSum_rambda_component_n_miuns_1  * exp(-beta * tau)

}


darrival <- function(...){

  # Extract named arguments from ...
  args <- list(...)

  # Assign expected arguments to local variables
  rambda_component <- args$rambda_component
  mu <- args$mu
  beta <- args$beta
  dimens <- args$dimens
  type_col_map <- args$type_col_map
  rresidual <- args$rresidual
  dresidual <- args$dresidual
  presidual <- args$presidual
  n <- args$n
  current_component <- args$current_component
  if (is.null(current_component)) {
    current_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)
  }
  if(is.null(args$tol)) tol <- 1e-09


  if (args$all_rows_beta_equal){

    #mu can be zero, so warning is turned off for a moment

    rowSums_rambda_component_n_miuns_1 <- rowSums(current_component)

    residuals_for_each_dimension <- rresidual(dimens)

    candidate_arrival <- numeric(length = dimens)

    for (i in 1:dimens){

      candidate_arrival[i] <- invert_exp_hazard(residuals_for_each_dimension[i],
        rowSums_rambda_component_n_miuns_1[i], mu[i], beta[i, 1])

    }

    inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE)

    type <- if (is.finite(inter_arrival)) minIndex[1] else NA_real_


    #if(length(type_col_map) > 0){
    #  for (j in 1:length(type_col_map)) {
    #    if (minIndex[1] %in% type_col_map[[j]]){
    #      type <- j
    #      break
    #    }
    #  }
    #}


  } else {

    stop("All elements of beta[i,] should be identical.")

  }

  c(inter_arrival = unname(inter_arrival), type = unname(type))
}



#' @title Expected Inter-Arrival Time
#' @description Computes the conditional expected time until the next event.
#' @param object An object of class \code{hspec}.
#' @param rambda_component Finite nonnegative scalar: the current total excitation
#' for the selected process dimension (sum its kernel components first).
#' @param type Process dimension index (default is 1).
#' @param mu Optional mu value (overrides object@mu if provided).
#' @param beta Optional beta value (overrides object@beta if provided).
#' @param tol Relative tolerance for numerical integration.
#' @param max_upper Upper integration limit.
#' @param subdivisions Number of subdivisions for numerical integration.
#' @return Expected value of next inter-arrival time.
#' @details
#' This calculation assumes a constant baseline and a common decay rate within
#' the selected row. With zero baseline, the result is `Inf` if there is positive
#' probability of no further event. A custom residual law then requires
#' `presidual` to determine this probability. A finite `max_upper` instead returns
#' the truncated first moment, excluding waiting times above that limit.
#'
#' @rdname expected_tau
#' @aliases expected_tau expected_tau,hspec-method
#' @export
setGeneric("expected_tau", function(object, rambda_component, type = 1,
                                    mu = NULL, beta = NULL,
                                    tol = .Machine$double.eps^0.25,
                                    max_upper = Inf,
                                    subdivisions = 400L)
  standardGeneric("expected_tau"))

setMethod("expected_tau", signature(object = "hspec"),
          function(object, rambda_component, type = 1,
                   mu = NULL, beta = NULL,
                   tol = .Machine$double.eps^0.25,
                   max_upper = Inf,
                   subdivisions = 400L) {

            if (!is.numeric(type) || is.complex(type) || length(type) != 1L ||
                !is.finite(type) || type != floor(type) || type < 1 || type > object@dimens) {
              stop("type must identify one process dimension.")
            }
            if (is.null(beta) && (is.function(object@beta) ||
                any(object@beta[type, ] != object@beta[type, 1]))) {
              stop("Provide a scalar beta or use equal decay rates within the selected row.")
            }
            if (is.null(mu) && is.function(object@mu)) stop("Provide a numeric scalar mu for a function-valued baseline.")
            # Use slot values if mu or beta not provided.
            mu <- if (is.null(mu)) object@mu[type] else mu
            beta <- if (is.null(beta)) object@beta[type, 1] else beta
            scalar <- function(x) is.numeric(x) && !is.complex(x) &&
              length(x) == 1L && is.finite(x)
            if (!scalar(mu) || mu < 0 || !scalar(beta) || beta <= 0 ||
                !scalar(rambda_component) || rambda_component < 0) {
              stop("mu and total excitation must be finite nonnegative scalars; beta must be finite and positive.")
            }
            if (!is.numeric(max_upper) || is.complex(max_upper) || length(max_upper) != 1L ||
                is.na(max_upper) || max_upper < 0) stop("max_upper must be nonnegative or Inf.")
            if (max_upper == 0) return(0)
            if (is.infinite(max_upper) && mu == 0) {
              if (is.null(object@dresidual)) return(Inf)
              if (is.null(object@presidual)) stop("Zero-baseline expectations require presidual for a custom residual law.")
              survival <- residual_cdf(object@presidual, rambda_component / beta,
                                       lower.tail = FALSE, log.p = TRUE)
              if (length(survival) != 1L || is.na(survival) || survival > 0)
                stop("presidual must return a valid scalar probability.")
              if (is.finite(survival)) return(Inf)
            }

            # E[T] = integral S(t) dt. Log-time coordinates resolve both the
            # short excitation scale and a long baseline tail without trying
            # to integrate a broad density on the original time scale.
            if (is.infinite(max_upper) &&
                (is.null(object@dresidual) || !is.null(object@presidual))) {
              if (is.null(object@dresidual) && rambda_component == 0) return(1 / mu)
              if (is.null(object@dresidual) &&
                  -rambda_component / beta - log(mu) > log(.Machine$double.xmax)) return(Inf)
              breaks <- sort(unique(c(-Inf, -log(mu), -log(beta),
                -log(max(mu, rambda_component)) - log1p(min(mu, rambda_component) / max(mu, rambda_component)), Inf)))
              survival_integrand <- function(log_time) {
                scaled_time <- exp(log(beta) + log_time)
                log_weight <- log(-expm1(-scaled_time)) - log(beta)
                log_weight[scaled_time == 0] <- log_time[scaled_time == 0]
                hazard <- exp(log(mu) + log_time) + exp(log(rambda_component) + log_weight)
                log_survival <- if (is.null(object@dresidual)) -hazard else
                  residual_cdf(object@presidual, hazard, lower.tail = FALSE, log.p = TRUE)
                exp(log_time + log_survival)
              }
              pieces <- vapply(seq_len(length(breaks) - 1L), function(i)
                stats::integrate(survival_integrand, breaks[i], breaks[i + 1L],
                  rel.tol = tol, abs.tol = 0, subdivisions = subdivisions)$value, numeric(1))
              return(sum(pieces))
            }

            # Integrate t*f(t) in log time. Raw-time quadrature can miss a
            # narrow density near zero, and its default absolute tolerance can
            # accept an answer with arbitrarily large relative error.
            if (mu == 0 && rambda_component == 0) return(0)
            upper <- log(max_upper)
            breaks <- sort(unique(c(-Inf, -log(mu), -log(beta),
              -log(max(mu, rambda_component)) -
                log1p(min(mu, rambda_component) / max(mu, rambda_component)))))
            breaks <- c(breaks[breaks < upper], upper)
            density_integrand <- function(log_time) {
              scaled_time <- exp(log(beta) + log_time)
              log_weight <- log(-expm1(-scaled_time)) - log(beta)
              log_weight[scaled_time == 0] <- log_time[scaled_time == 0]
              hazard <- exp(log(mu) + log_time) +
                exp(log(rambda_component) + log_weight)
              log_density <- if (is.null(object@dresidual)) stats::dexp(hazard, log = TRUE) else
                residual_density(object@dresidual, hazard, log = TRUE)
              log_component <- log(rambda_component) - scaled_time
              largest <- pmax(log(mu), log_component)
              log_intensity <- rep(-Inf, length(log_time))
              finite <- is.finite(largest)
              log_intensity[finite] <- largest[finite] +
                log1p(exp(pmin(log(mu), log_component[finite]) - largest[finite]))
              value <- exp(2 * log_time + log_density + log_intensity)
              value[!finite] <- 0
              value
            }
            pieces <- vapply(seq_len(length(breaks) - 1L), function(i)
              stats::integrate(density_integrand, breaks[i], breaks[i + 1L],
                rel.tol = tol, abs.tol = 0, subdivisions = subdivisions)$value, numeric(1))
            sum(pieces)
          })
