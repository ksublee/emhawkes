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


  if (args$all_rows_beta_equal){

    #mu can be zero, so warning is turned off for a moment
    candidate_arrival <- suppressWarnings(stats::rexp(dimens, rate = mu))

    rowSums_rambda_component_n_miuns_1 <- rowSums(matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE))

    rowSums_rambda_component_n_miuns_1[abs(rowSums_rambda_component_n_miuns_1) < 1e-10] <- 0

    matrixD <- 1 + beta[,1] * log(stats::runif(dimens)) / rowSums_rambda_component_n_miuns_1

    candidate_arrival <- cbind(candidate_arrival, -1 / beta[,1] * log(pmax(matrixD, 0)))

    inter_arrival <- min(candidate_arrival, na.rm = TRUE)

    minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE)

    type <- minIndex[1]

    if(length(type_col_map) > 0){
      for (j in 1:length(type_col_map)) {
        if (minIndex[1] %in% type_col_map[[j]]){
          type <- j
          break
        }
      }
    }


  } else{

    # thinning algorithm

    rambda_component_n_miuns_1 <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)
    rowSums_rambda_component_n_miuns_1 <- rowSums(rambda_component_n_miuns_1)
    rambda_n_minus_1 <- rowSums_rambda_component_n_miuns_1 + mu

    rambda_n_minus_1[abs(rambda_n_minus_1) < 1e-10] <- 0

    candidate_arrival <- rep(0, dimens)

    # Loop through each component
    for(i in 1:dimens){

      s <- 0 # Initialize time s, which will represent the time to the next event.

      if(rambda_n_minus_1[i] == 0){

        candidate_arrival[i] <- Inf

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
    type <- which(candidate_arrival == inter_arrival)

  }

  c(inter_arrival = inter_arrival, type = type)
}


# for flexible model arrival

phi <- function(tau, rowSum_rambda_component_n_miuns_1, mu, beta){

  # rambda means after jump intensity
  # for arrival time n, row sum of rambda component at n-1 is needed

  mu * tau + rowSum_rambda_component_n_miuns_1  * (-exp(-beta * tau) + 1) / beta

}

phi_wo_mu <- function(tau, rowSum_rambda_component_n_miuns_1,  beta){

  # rambda means after jump intensity
  # for arrival time n, row sum of rambda component at n-1 is needed

  rowSum_rambda_component_n_miuns_1  * (-exp(-beta * tau) + 1) / beta

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
  if(is.null(args$tol)) tol <- 1e-09


  if (args$all_rows_beta_equal){

    #mu can be zero, so warning is turned off for a moment

    rowSums_rambda_component_n_miuns_1 <- rowSums(matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE))

    residuals_for_each_dimension <- rresidual(dimens)

    candidate_arrival <- numeric(length = dimens)

    for (i in 1:dimens){

      candidate_arrival[i] <- uniroot(function(x) phi(x,
                                                      rowSums_rambda_component_n_miuns_1[i],
                                                      mu[i], beta[i, 1]) - residuals_for_each_dimension[i],
                                      interval = c(0, 10000), tol = 1e-09)$root

    }

    inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE)

    type <- minIndex[1]


    #if(length(type_col_map) > 0){
    #  for (j in 1:length(type_col_map)) {
    #    if (minIndex[1] %in% type_col_map[[j]]){
    #      type <- j
    #      break
    #    }
    #  }
    #}


  } else {

    Stop("All elements of beta[i,] should be identical.")

  }

  c(inter_arrival = inter_arrival, type = type)
}



#' @title Expected Inter-Arrival Time
#' @description Computes the conditional expected time until the next event.
#' @param object An object of class \code{hspec}.
#' @param rambda_component Rambda component.
#' @param type Process dimension index (default is 1).
#' @param mu Optional mu value (overrides object@mu if provided).
#' @param beta Optional beta value (overrides object@beta if provided).
#' @param tol Relative tolerance for numerical integration.
#' @param max_upper Upper integration limit.
#' @param subdivisions Number of subdivisions for numerical integration.
#' @return Expected value of next inter-arrival time.
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

            # Use slot values if mu or beta not provided
            mu <- if (is.null(mu)) object@mu[type] else mu
            beta <- if (is.null(beta)) object@beta[type, 1] else beta

            # Handle residual density
            if (is.null(object@dresidual)) {
              dresidual_fn <- function(x, param = NULL) dexp(x)
              default_param <- NULL
            } else {
              dresidual_fn <- object@dresidual
              default_param <- formals(dresidual_fn)$param
            }

            # Define integrand
            integrand <- function(s) {
              phi_val <- phi(s, rambda_component, mu, beta)
              dens_val <- dresidual_fn(phi_val, param = default_param)
              s * dens_val * psi(s, rambda_component, mu, beta)
            }

            # Numerical integration
            result <- integrate(integrand, lower = 0, upper = max_upper,
                                rel.tol = tol, subdivisions = subdivisions)

            return(result$value)
          })
