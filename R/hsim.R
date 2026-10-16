#' @include hspec.R hmoment.R utilities.R hinterface.R
NULL

#' Simulate an Exponential Hawkes Process
#'
#' @description
#' Simulate a univariate or multivariate path from an \code{\link{hspec-class}}
#' specification, including marked and multi-kernel models.
#'
#'
#' @param object \code{\link{hspec-class}}. S4 object that specifies the parameter values.
#' @param size Positive integer giving the number of output rows, including the
#' initial state. At most `size - 1` events are simulated.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' A finite scalar is expanded to every kernel component. A matrix must match
#' the shape of `alpha` and `beta`. Supply this argument when `type_col_map` is set.
#' Other size adjustments are
#' deprecated and produce a warning.
#' @param N0 Initial event counts: a finite scalar or one value per event type.
#' The default is zero.
#' @param Nc0 Initial mark totals: a finite scalar or one value per event type.
#' The default is zero.
#' @param verbose Logical. If \code{TRUE}, print progress messages during the simulation. Default is \code{FALSE}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An S3 \code{\link{hreal}} object containing the simulated path,
#' counting processes, and intensity components.
#' If no further event can occur, the result contains only the generated rows
#' and a warning reports the early termination. The initial row is included.
#'
#' @details
#' The first row is the initial state. Its `inter_arrival`, `type`, and `mark`
#' entries are zero; its `N` and `Nc` rows contain `N0` and `Nc0`.
#' Event types identify intensity rows. After a type `i` event, `type_col_map[[i]]`
#' selects the kernel columns to update; without a mapping, column `i` is updated.
#'
#' Simulation requires finite, nonnegative baselines. Signed kernel components
#' are allowed when their sum is nonnegative and decreasing between events.
#' For equal decay rates within a row, this requires a nonnegative row sum.
#' For unequal rates, the simulator uses a sufficient condition: group components
#' by decay rate, multiply each group sum by its rate, and check that cumulative
#' sums in increasing rate order are nonnegative.
#'
#' Simulation automatically uses an Rcpp event loop for plain numeric `mu`,
#' `alpha`, and `beta` and positive baselines in every dimension. This includes
#' multivariate and mapped multi-kernel models with equal or unequal decay rates
#' within a row. Equal-rate models use direct arrival sampling; unequal-rate
#' models retain the R algorithm's thinning and rejected draws. Models with
#' `eta`, component functions, event callbacks, custom residual distributions,
#' or zero baselines use R. This selection is independent of `emhawkes.engine`,
#' which controls likelihood calculations. Both paths use R's RNG in the same order.
#' Arrival calculations use the current component state.
#' Callbacks receive the arguments they request; callbacks with `...` continue
#' to receive the full supported context. Requesting full histories can increase
#' copying costs for long paths.
#'
#' To simulate an estimated model, use [fitted_hspec()] and
#' pass [fit_initial_state()] to `hsim()` with `do.call()`. These arguments describe
#' the start of the fitted observation window. Continuing after its last event
#' requires the state at that event.
#' @seealso [fitted_hspec()], [fit_initial_state()], [logLik()]
#'
#' @rdname hsim
#' @docType methods
#' @examples
#'
#' # example 1
#'
#' mu <- 1; alpha <- 1; beta <- 2
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' hsim(h, size=100)
#'
#'
#' # example 2
#' mu <- matrix(c(0.1, 0.1), nrow=2)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' print(res)
#'
#' @export
setGeneric("hsim", function(object, size = 100,
                            lambda_component0 = NULL,  N0 = NULL, Nc0 = NULL,
                            verbose = FALSE, ...)
  standardGeneric("hsim"))

#' @rdname hsim
setMethod("hsim", "hspec", function(object, size = 100,
    lambda_component0 = NULL, N0 = NULL, Nc0 = NULL, verbose = FALSE, ...) {
  simulate_hawkes(object, size, lambda_component0, N0, Nc0, verbose, ...)
})

# Internal reference selector: independent of the likelihood engine option.
simulate_hawkes_r <- function(object, size = 100, lambda_component0 = NULL,
    N0 = NULL, Nc0 = NULL, verbose = FALSE, ...) {
  simulate_hawkes(object, size, lambda_component0, N0, Nc0, verbose, ...,
                  .native = FALSE)
}

simulate_hawkes <- function(object, size = 100, lambda_component0 = NULL,
    N0 = NULL, Nc0 = NULL, verbose = FALSE, ..., .native = TRUE) {
    # Process additional arguments passed through '...'
    additional_argument <- list(...)
    if (length(size) != 1L || !is.finite(size) || size < 1 || size != floor(size)) {
      stop("size must be a positive integer, including the initial row.")
    }

    # Check if "lambda0" is among the additional arguments
    if ("lambda0" %in% names(additional_argument)) {

      # Issue a warning to inform the user that "lambda0" is deprecated
      warning("lambda0 is deprecated; instead use lambda_component0.")

      # Assign the value of "lambda0" to the new parameter name "lambda_component0"
      lambda_component0 <- additional_argument[["lambda0"]]
    }

    # pre-setting, especially for mu, alpha, beta, eta.
    # Useful when those are defined as functions.
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta
    callbacks <- model_callback_plans(object, simulation = TRUE)
    mu_call <- prepare_component_call(callbacks$mu)
    impact_call <- prepare_component_call(callbacks$impact)
    rmark_call <- prepare_component_call(callbacks$rmark)

    # Check for each type, beta are equal
    all_rows_beta_equal <- TRUE
    # Iterate over each row to check if all elements are identical
    for (i in seq_len(nrow(beta))) {
      # Check if all elements in the row are equal to the first element
      if (!all(beta[i, ] == beta[i, 1])) {
        all_rows_beta_equal <- FALSE
        break # Exit the loop early if any row is not identical
      }
    }

    dimens <- object@dimens
    rate_groups <- lapply(seq_len(dimens), function(i) {
      split(seq_len(ncol(beta)), factor(beta[i, ], levels = sort(unique(beta[i, ]))))
    })
    valid_components <- function(component) {
      if (any(!is.finite(component))) return(FALSE)
      if (all(component >= 0)) return(TRUE)
      component <- matrix(component, nrow = dimens)
      for (i in seq_len(dimens)) {
        if (all(component[i, ] >= 0)) next
        # Equal-rate resets can cancel to zero up to floating-point error.
        if (length(rate_groups[[i]]) == 1L) {
          weight <- component[i, ]
          cumulative <- sum(weight)
        } else {
          weight <- component[i, ] * beta[i, ]
          cumulative <- cumsum(vapply(rate_groups[[i]], function(j) sum(weight[j]), numeric(1)))
        }
        scale <- sum(abs(weight))
        # Abel summation: nonnegative cumulative derivative weights ensure a
        # decreasing excitation, whose limit is zero. It is thus nonnegative.
        if (!is.finite(scale) || any(!is.finite(cumulative)) ||
            any(cumulative < -32 * .Machine$double.eps * max(1, scale))) return(FALSE)
      }
      TRUE
    }

    N <- matrix(numeric(length = dimens * size), ncol = dimens)
    Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)
    colnames(Nc)  <- paste0("Nc", 1:dimens)
    colnames(N) <- paste0("N", 1:dimens)

    type <- inter_arrival <- numeric(length = size)
    mark <- c(0, rep(1, size-1))

    N0 <- normalize_count_initial(N0, dimens, "N0")
    Nc0 <- normalize_count_initial(Nc0, dimens, "Nc0")
    N[1,] <- N0
    Nc[1,] <- Nc0
    lambda <- lambda_component <- lambda_component_n <- NULL

    if (is.function(mu)){
      mu0 <- mu_call(mu, 1, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta)
    } else {
      mu0 <- mu
    }

    # default lambda_component0
    if(!is.null(object@type_col_map)){
      if(length(object@type_col_map) > 0 & is.null(lambda_component0)){
        stop("Initialization Error: 'lambda_component0' must be specified.
             This ensures the model starts with a valid configuration.
             Please review the input parameters and provide valid initial values for 'lambda_component0'.")
        }
    }


    if(!is.null(lambda_component0)){

      lambda_component0 <- normalize_lambda_initial(lambda_component0, dimens, ncol(beta), legacy = TRUE)

    } else {
      # default lambda_component0
      if(verbose == TRUE){
        message("Default initial values for 'lambda_component0' will be internally calculated and used for the simulation.")
      }

      lambda_component0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu, alpha = alpha, beta = beta)
    }

    lambda_component0 <- normalize_lambda_initial(lambda_component0, dimens, ncol(beta))
    lambda_component_n <- lambda_component0
    if (!valid_components(lambda_component0)) {
      stop("Initial kernel components must define nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).")
    }

    # Decide support before entering any RNG scope; never restart after draws.
    if (.native && is.null(native_simulation_reason(object, mu, alpha, beta,
                                                   all_rows_beta_equal))) {
      mapping <- if (is.null(object@type_col_map)) as.list(seq_len(dimens)) else object@type_col_map
      mapping <- lapply(mapping, function(columns) as.integer(columns) - 1L)
      out <- hawkes_simulate_cpp(as.numeric(mu), alpha, beta, lambda_component0,
        mapping, N, Nc, inter_arrival, type, mark)
      colnames(out$lambda) <- paste0("lambda", seq_len(dimens))
      colnames(out$rambda) <- paste0("rambda", seq_len(dimens))
      component_names <- as.vector(t(outer(paste0("lambda", seq_len(dimens)),
        seq_len(ncol(beta)), FUN = paste0)))
      colnames(out$lambda_component) <- component_names
      colnames(out$rambda_component) <- sub("^lambda", "rambda", component_names)
      if (out$generated < size) {
        rows <- seq_len(out$generated)
        for (name in c("inter_arrival", "type", "mark")) out[[name]] <- out[[name]][rows]
        for (name in c("N", "Nc", "lambda", "rambda", "lambda_component", "rambda_component"))
          out[[name]] <- out[[name]][rows, , drop = FALSE]
      }
      return(structure(list(hspec = object, inter_arrival = out$inter_arrival,
        arrival = cumsum(out$inter_arrival), type = out$type, mark = out$mark,
        N = out$N, Nc = out$Nc, lambda = out$lambda,
        lambda_component = out$lambda_component, rambda = out$rambda,
        rambda_component = out$rambda_component), class = "hreal"))
    }

    # Preallocation for lambdas and Ns and set initial values for lambdas
    #lambda_component <- matrix(sapply(lambda_component0, c, numeric(length = size - 1)),
    #                           ncol = dimens * ncol(beta))

    lambda_component <- matrix(0, nrow = size, ncol = ncol(beta) * dimens)
    # Set the initial values for the first row based on lambda_component0
    lambda_component[1, ] <- t(lambda_component0)



    # Initialize the lambda matrix with zeros
    lambda <- matrix(0, nrow = size, ncol = dimens)
    rowSums_lambda_component0 <- rowSums(matrix(lambda_component0, nrow=dimens))

    # Set the initial values for the first row based on the initial sums
    lambda[1, ] <- mu0 + rowSums_lambda_component0

    #lambda <- matrix(sapply(mu0 + rowSums_lambda_component0, c, numeric(length = size - 1)),
    #                 ncol = dimens)

    rambda <- lambda
    rambda_component <- lambda_component

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:ncol(beta), dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))

    colnames(rambda) <- paste0("rambda", 1:dimens)
    colnames(rambda_component) <- as.vector(t(outer(as.vector(outer("rambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))


    # Before entering the loop, compute mu after jump if exists,
    # This is equivalent to right continuous version of mu at n = 2.
    # This is intensity for interarrival between n=1 and n=2.
    # Mathematically, mu is a left continuous process, therefore, this is equivalent to mu(1 + dt)
    if (is.function(mu)){
      # Even if, we pass n=2 as argument, mu function probably use information up to n-1.
      rmu_n <- mu_call(mu, 2, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta)
    } else{
      rmu_n <- mu
    }

    # Helper function to handle arrival logic

    if (is.null(object@rresidual)) {

      handle_arrival <- rarrival

    } else {

      handle_arrival <- darrival

    }



    # Keep the state used for the next event separate from stored output.
    current_rambda_component <- lambda_component0
    generated <- 1L
    for (n in seq_len(size - 1L) + 1L) {

      # only piecewise constant mu is available, mu is left continuous.
      # rmu_n is previously defined right continuous version.
      mu_n <- rmu_n
      if (length(mu_n) != dimens || any(!is.finite(mu_n)) || any(mu_n < 0) ||
          !valid_components(current_rambda_component)) {
        stop("Simulation requires finite, nonnegative baselines and nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).")
      }


      ### Determine the next arrival, type and mark #######################################################
      res <- handle_arrival(n = n, current_component = current_rambda_component,
                            mu = mu_n, alpha = alpha, beta = beta, dimens = object@dimens,
                            type_col_map = object@type_col_map,
                            rresidual = object@rresidual, dresidual = object@dresidual,
                            presidual = object@presidual, all_rows_beta_equal = all_rows_beta_equal)

      if (is.infinite(res["inter_arrival"])) {
        warning("Simulation stopped early: no further event occurs on this path.")
        break
      }


      inter_arrival[n] <- res["inter_arrival"]
      type[n] <- res["type"]
      generated <- n


      N[n, ] <- N[n-1, ]
      N[n, type[n]] <- N[n-1, type[n]] + 1

      # generate a mark for Hawkes
      if( !is.null(object@rmark) ){
        # mark may depends on other variables
        mark[n] <- rmark_call(object@rmark, n, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta)
      }

      Nc[n, ] <- Nc[n-1, ]
      Nc[n, type[n]] <- Nc[n-1, type[n]] + mark[n]
      ######################################################################################################



      # lambda decayed due to time, impact due to mark is not added yet

      decayed <- exp(-beta * inter_arrival[n])

      decayed_lambda <- lambda_component_n <- current_rambda_component * decayed

      # update lambda
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu_n + rowSums(decayed_lambda)


      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)

      if( is.null(object@type_col_map)){
        types <- type[n]
      } else{
        types <- object@type_col_map[[type[n]]]
      }

      impact_alpha[ , types] <- alpha[ , types]

      new_lambda <- decayed_lambda + impact_alpha

      # additional impact by eta
      if(!is.null(eta)){

        impact_eta <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        new_lambda <- new_lambda + impact_eta
      }


      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(object@impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)

        impact_res <- impact_call(object@impact, n, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta)

        impact_mark[ , types] <- impact_res[ , types]

        new_lambda <- new_lambda + impact_mark

      }

      # update rambda
      # rambda_component = {"rambda11", "rambda12", ..., "rambda21", "rambda22", ...}
      rambda_component[n, ] <- t(new_lambda)
      # The history row used previously had no matrix dimnames.
      current_rambda_component <- matrix(new_lambda, nrow = dimens)



      if (is.function(mu)){
        rmu_n <- mu_call(mu, n + 1, mark, type, inter_arrival, N, Nc, lambda,
                     lambda_component, lambda_component_n, mu, alpha, beta)
      } else{
        rmu_n <- mu
      }
      rambda[n, ] <- rmu_n + rowSums(new_lambda)

    }

    rows <- seq_len(generated)
    inter_arrival <- inter_arrival[rows]; type <- type[rows]; mark <- mark[rows]
    N <- N[rows, , drop = FALSE]; Nc <- Nc[rows, , drop = FALSE]
    lambda <- lambda[rows, , drop = FALSE]
    rambda <- rambda[rows, , drop = FALSE]
    lambda_component <- lambda_component[rows, , drop = FALSE]
    rambda_component <- rambda_component[rows, , drop = FALSE]

    realization <- list(object, inter_arrival, cumsum(inter_arrival), type, mark,
                        N, Nc, lambda, lambda_component, rambda, rambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "type", "mark",
                            "N", "Nc", "lambda", "lambda_component", "rambda", "rambda_component")
    class(realization) <- c("hreal")

    realization
}

# Conservative support check. A function returning numeric parameters still uses R.
native_simulation_reason <- function(object, mu, alpha, beta, equal) {
  values <- list(object@mu, object@alpha, object@beta)
  if (any(vapply(values, is.function, logical(1)))) return("function_parameters")
  if (any(vapply(values, is.object, logical(1)))) return("classed_parameters")
  if (!is.null(object@eta)) return("eta")
  if (!is.null(object@impact) || !is.null(object@rmark)) return("event_callback")
  if (!is.null(object@rresidual) || !is.null(object@dresidual) ||
      !is.null(object@presidual) || !is.null(object@qresidual)) return("residual_distribution")
  if (any(!is.finite(mu)) || any(mu <= 0)) return("nonpositive_baseline")
  if (any(!is.finite(alpha)) || any(!is.finite(beta)) || any(beta <= 0)) return("invalid_parameters")
  if (is.object(object@type_col_map) || any(vapply(object@type_col_map, is.object, logical(1))))
    return("classed_mapping")
  NULL
}
