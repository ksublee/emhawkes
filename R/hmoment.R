#' @include hspec.R hinterface.R
NULL

setGeneric("get_lambda0", function(object, ...) standardGeneric("get_lambda0"))
# Compute component means after a mean policy has already been selected.
setMethod(
  f = "get_lambda0",
  signature(object = "hspec"),
  definition = function(object, inter_arrival = inter_arrival, type = type, mark = mark,
                        N = N, Nc = Nc,
                        mu = mu, alpha = alpha, beta = beta,
                        lambda = NULL, lambda_component = NULL, lambda_component_n = NULL,
                        ..., .legacy = FALSE){

    dimens <- object@dimens

    # parameter setting
    if (is.function(object@mu)){
      if(identical(names(formals(object@mu)), "param")){
        mu0 <- evalf(object@mu)
      } else {
        context <- list(n = 1, inter_arrival = inter_arrival, type = type, mark = mark,
                        N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                        lambda_component_n = lambda_component_n, mu = mu, alpha = alpha, beta = beta)
        mu0 <- call_component(mu, callback_plan(mu, names(context), "mu"), context)
      }
    } else{
      mu0 <- object@mu
    }
    if (is.function(object@alpha)){
      alpha <- evalf(object@alpha)
    } else{
      alpha <- object@alpha
    }
    if (is.function(object@beta)){
      beta <- evalf(object@beta)
    } else{
      beta <- object@beta
    }


    # Only replay of saved fits from before explicit initialization policies uses
    # the old half-mean and core-only approximation. New auto fits never use it.
    if (.legacy) {
      if (dimens == 1L) return(matrix((mu0[1] * beta[1] / (beta[1] - alpha[1]) - mu0[1]) / 2))
      stationary <- tryCatch(solve(diag(dimens) - alpha / beta) %*% mu0,
        error = function(e) NULL)
      if (is.null(stationary)) return(matrix(0, dimens, dimens))
      return(matrix(rep(stationary, dimens), dimens, byrow = TRUE) * alpha / beta)
    }
    result <- hawkes_mean_state(list(mu = mu0, alpha = alpha, beta = beta),
                                object@type_col_map)
    if (is.null(result)) stop("Mean initial intensity components require a stable nonnegative linear Hawkes model with finite component means; supply lambda_component0 explicitly or use 'zero'.", call. = FALSE)
    result$value
  }
)

# Mean values and the quantities needed for their analytic chain rule.
# A kernel column may receive events from several mapped event types.
hawkes_mean_state <- function(values, mapping = NULL) {
  mu <- values$mu; alpha <- values$alpha; beta <- values$beta
  if (any(!is.finite(c(mu, alpha, beta))) || any(mu < 0) ||
      any(alpha < 0) || any(beta <= 0)) return(NULL)
  dimensions <- nrow(beta); kernels <- ncol(beta)
  if (is.null(mapping)) mapping <- as.list(seq_len(dimensions))
  incidence <- matrix(0, kernels, dimensions)
  for (j in seq_len(dimensions)) incidence[mapping[[j]], j] <- 1
  branching <- alpha / beta
  integrated <- branching %*% incidence
  if (any(!is.finite(c(branching, integrated)))) return(NULL)
  radius <- if (dimensions == 1L) integrated[1] else
    max(Mod(eigen(integrated, only.values = TRUE)$values))
  if (!is.finite(radius) || radius >= 1) return(NULL)
  inverse <- tryCatch(solve(diag(dimensions) - integrated), error = function(e) NULL)
  if (is.null(inverse)) return(NULL)
  stationary <- inverse %*% mu
  rates <- as.vector(incidence %*% stationary)
  value <- matrix(rep(rates, dimensions), dimensions, byrow = TRUE) * branching
  if (any(!is.finite(value)) || any(value < 0)) return(NULL)
  list(value = value, inverse = inverse, stationary = stationary,
       rates = rates, branching = branching, incidence = incidence)
}

# Select once, before simulation draws or optimization. In internal objectives,
# NULL denotes parameter-dependent mean components; a matrix is a fixed state.
select_hawkes_initial <- function(object, initial, kernels, legacy = FALSE) {
  if (is.null(initial)) initial <- "auto"
  if (!is.character(initial)) return(list(requested = "explicit", used = "explicit",
    value = if (is.null(kernels)) initial else normalize_lambda_initial(initial, object@dimens, kernels, legacy), reason = NULL))
  if (length(initial) != 1L || is.na(initial) || !initial %in% c("auto", "mean", "zero"))
    stop("lambda_component0 must be 'auto', 'mean', 'zero', a finite scalar, or a kernel-shaped matrix.", call. = FALSE)
  requested <- initial
  values <- list(mu = object@mu, alpha = object@alpha, beta = object@beta)
  constant <- all(vapply(values, function(x) is.numeric(x) && !is.object(x), logical(1)))
  extra <- is.function(object@eta) || (!is.null(object@eta) && any(object@eta != 0)) ||
    !is.null(object@impact) || any(vapply(list(object@rresidual, object@dresidual,
      object@presidual, object@qresidual), function(x) !is.null(x), logical(1)))
  reason <- if (!constant || extra) "no supported unconditional mean formula for this model" else NULL
  mean <- if (initial != "zero" && is.null(reason)) hawkes_mean_state(values, object@type_col_map) else NULL
  if (initial != "zero" && is.null(reason) && is.null(mean))
    reason <- "finite stationary component means are unavailable at the starting parameters"
  if (initial == "mean" && !is.null(reason))
    stop("Mean initialization is unavailable: ", reason,
      ". Supply lambda_component0 explicitly or use 'zero'.", call. = FALSE)
  used <- if (initial != "zero" && !is.null(mean)) "mean" else "zero"
  list(requested = requested, used = used,
       value = if (used == "mean") NULL else if (is.null(kernels)) 0 else matrix(0, object@dimens, kernels),
       reason = if (used == "zero" && requested == "auto") reason else NULL)
}

message_hawkes_initial <- function(plan) {
  if (plan$used == "explicit") return(invisible(NULL))
  selection <- if (plan$requested == "auto") paste0("auto -> ", plan$used) else plan$used
  message("Initial intensity components: ", selection,
    if (plan$used == "mean") " (unconditional component means)." else
      paste0(if (!is.null(plan$reason)) paste0(" (", plan$reason, ")") else "", ". ",
        "Supplying lambda_component0 explicitly is recommended for this model."))
  invisible(NULL)
}


#' Compute Hawkes volatility
#'
#' This function computes Hawkes volatility for a constant bivariate Hawkes model
#' with one kernel per event type and a common decay rate within each row.
#'
#' @param object \code{\link{hspec-class}}
#' @param horizon Time horizon used to scale the stationary long-run variance rate.
#' @param inter_arrival Inter-arrival times of events which includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type  A vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark A vector of mark (jump) sizes. Start with zero. If omitted, all
#' observed events have unit marks. At least one event of each type is required.
#' @param dependence Dependence between mark and previous sigma-algebra.
#' @param lambda_component0 A matrix of initial intensity components, excluding
#' the baseline `mu`.
#' @param ... Further arguments passed to or from other methods.
#' @details
#' Numeric `mu`, `alpha`, and `beta` are required. A missing `eta` is treated as
#' zero. Custom impact functions, residual distributions, and noncanonical kernel
#' mappings are not supported by this moment formula.
#' The result is the square root of `horizon` times the stationary long-run
#' variance rate of the difference between the two marked counting processes.
#'
#' @rdname hvol
#' @export
setGeneric("hvol", function(object, horizon = 1,
                           inter_arrival = NULL,
                           type = NULL,
                           mark = NULL,
                           dependence = FALSE,
                           lambda_component0 = NULL, ...) standardGeneric("hvol"))

#' @rdname hvol
setMethod(
  f = "hvol",
  signature(object = "hspec"),
  definition = function(object, horizon = 1,
                        inter_arrival = NULL,
                        type = NULL,
                        mark = NULL,
                        dependence = FALSE,
                        lambda_component0 = NULL,
                        ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("'lambda0' is deprecated; use 'lambda_component0' instead.", call. = FALSE)

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    h <- object

    if (!is.numeric(horizon) || is.complex(horizon) || length(horizon) != 1L ||
        !is.finite(horizon) || horizon < 0) {
      stop("horizon must be a finite nonnegative scalar.")
    }
    if (!is.logical(dependence) || length(dependence) != 1L || is.na(dependence)) {
      stop("dependence must be TRUE or FALSE.")
    }
    if (h@dimens != 2) stop("hvol requires a bivariate model (dimens = 2).", call. = FALSE)
    if (any(vapply(list(h@mu, h@alpha, h@beta, h@eta), is.function, logical(1))))
      stop("hvol requires constant numeric mu, alpha, beta and eta; function-valued parameters are unsupported.", call. = FALSE)
    if (!identical(dim(h@alpha), c(2L, 2L)) || !identical(dim(h@beta), c(2L, 2L)))
      stop("hvol requires one kernel per type: alpha and beta must be 2 x 2 matrices.", call. = FALSE)
    if (any(h@beta[, 1] != h@beta[, 2]))
      stop("hvol requires equal decay rates within each row of beta.", call. = FALSE)
    if (!is.null(h@impact)) stop("hvol does not support a custom impact function.", call. = FALSE)
    if (any(!vapply(list(h@dresidual, h@presidual, h@qresidual, h@rresidual), is.null, logical(1))))
      stop("hvol does not support a custom residual distribution.", call. = FALSE)
    if (!is.null(h@type_col_map) && !identical(lapply(h@type_col_map, as.integer), list(1L, 2L)))
      stop("hvol requires one kernel per type with type_col_map = NULL or list(1, 2).", call. = FALSE)
    observation <- normalize_observation(h, inter_arrival, type, mark, needs = character())
    inter_arrival <- observation$inter_arrival
    type <- observation$type
    type[1] <- 0
    mark <- observation$mark
    if (!all(c(1, 2) %in% type)) stop("hvol requires at least one observed event of each type.")
    infer_result <- logLik(h, inter_arrival = inter_arrival, type = type, mark = mark,
                          lambda_component0 = lambda_component0, infer = TRUE)

    Mu <- matrix(h@mu, nrow=2)
    Alpha <- h@alpha
    Beta <- diag(diag(h@beta))
    Eta <- if (is.null(h@eta)) matrix(0, 2, 2) else h@eta

    up_move_idx <- type == 1
    down_move_idx <- type == 2


    #### compute Hawkes vol

    if (dependence){

      Z_1 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1])

      Z_2 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2])

      Z_sq_1 <- sum(mark[up_move_idx]^2 * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1])

      Z_sq_2 <- sum(mark[down_move_idx]^2 * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2])

      Z_ll_11 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,1] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1] * infer_result$lambda[up_move_idx,1])

      Z_ll_12 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,1] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,1] * infer_result$lambda[down_move_idx,2])

      Z_ll_21 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,2] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,2] * infer_result$lambda[up_move_idx,1])

      Z_ll_22 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,2] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2] * infer_result$lambda[down_move_idx,2])


      Z <- matrix(c(Z_1, Z_2,
                    Z_1, Z_2), nrow=2, byrow=TRUE)

      Z_sq <- matrix(c(Z_sq_1, Z_sq_2,
                       Z_sq_1, Z_sq_2), nrow=2, byrow=TRUE)

      Z_ll <- matrix(c(Z_ll_11, Z_ll_12,
                       Z_ll_21, Z_ll_22), nrow=2, byrow=TRUE)

      Z_nl <- Z


      E_lambda <- solve(Beta - Alpha + Eta - Eta * Z) %*% Beta %*% Mu

      # To find E_ll : E_lambda_lambda

      G <- (Alpha - Eta + Eta * Z) %*% diag(c(E_lambda)) %*% t(Alpha - Eta) +
        (Alpha - Eta) %*% diag(c(E_lambda)) %*% t(Eta * Z) +
        (Eta * Z_sq^(0.5)) %*% diag(c(E_lambda)) %*% t(Eta * Z_sq^(0.5))


      H <- kronecker(diag(2), Alpha - Beta) + kronecker(Alpha - Beta, diag(2)) +
        rbind(cbind(Eta %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ll)[,2] - 1))) +
        rbind(cbind(t(Eta)[,1] %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,1] %*% diag(t(Z_ll)[,2] - 1)),
              cbind(t(Eta)[,2] %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,2] %*% diag(t(Z_ll)[,2] - 1)))


      C <- Beta %*% Mu %*% t(E_lambda) + t(Beta %*% Mu  %*% t(E_lambda)) + G

      E_ll <- matrix(solve(H, -c(C)), nrow=2)

      HA <- kronecker(diag(2), Alpha - Beta) +
        rbind(cbind(Eta %*% diag(t(Z)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z)[,2] - 1)))

      CA <- (Beta %*% Mu) %*% t(E_lambda) %*% diag(diag(Z))

      A <- t(matrix(solve(HA, -c(CA)), nrow=2))
      # A <- diag(diag(Z)) %*% E_lambda %*% t(E_lambda)

      CB <- Z_ll * E_ll + ((Alpha - Eta) * Z + Eta * Z_sq) %*% diag(c(E_lambda)) - t(A)

      B <- t(matrix(solve(HA, -c(CB)), nrow=2))

      E_nl <- A * horizon + B

      P1 <- Z_nl * 0.5 * A + t(Z_nl * 0.5 * A)
      P2 <- Z_nl * B  + t(Z_nl * B)
      P3 <- Z_sq * diag(c(E_lambda))

      E_NN <- P1 * horizon^2 + P2 * horizon + P3 * horizon

      variance <- matrix(c(1, -1), nrow=1) %*% (P2 * horizon + P3 * horizon) %*% matrix(c(1, -1))

      (std <- sqrt(variance))

    } else {

      ##### independent

      Z_ind_1 <- mean(mark[up_move_idx])
      Z_ind_2 <- mean(mark[down_move_idx])

      Z_sq_ind_1 <- mean(mark[up_move_idx]^2)
      Z_sq_ind_2 <- mean(mark[down_move_idx]^2)

      Z_ind <- matrix(c(Z_ind_1, Z_ind_2,
                        Z_ind_1, Z_ind_2), nrow=2, byrow=TRUE)

      Z_sq_ind <-  matrix(c(Z_sq_ind_1, Z_sq_ind_2,
                            Z_sq_ind_1, Z_sq_ind_2), nrow=2, byrow=TRUE)

      E_lambda_ind <- solve(Beta - Alpha + Eta - Eta * Z_ind) %*% Beta %*% Mu

      G_ind <- (Alpha - Eta + Eta * Z_ind) %*% diag(c(E_lambda_ind)) %*% t(Alpha - Eta) +
        (Alpha - Eta) %*% diag(c(E_lambda_ind)) %*% t(Eta * Z_ind) +
        (Eta * Z_sq_ind^(0.5)) %*% diag(c(E_lambda_ind)) %*% t(Eta * Z_sq_ind^(0.5))

      H_ind <- kronecker(diag(2), Alpha - Beta) + kronecker(Alpha - Beta, diag(2)) +
        rbind(cbind(Eta %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ind)[,2] - 1))) +
        rbind(cbind(t(Eta)[,1] %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,1] %*% diag(t(Z_ind)[,2] - 1)),
              cbind(t(Eta)[,2] %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,2] %*% diag(t(Z_ind)[,2] - 1)))

      C_ind <- Beta %*% Mu %*% t(E_lambda_ind) + t(Beta %*% Mu  %*% t(E_lambda_ind)) + G_ind

      E_ll_ind <- matrix(solve(H_ind, -c(C_ind)), nrow=2)

      HA_ind <- kronecker(diag(2), Alpha - Beta) +
        rbind(cbind(Eta %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ind)[,2] - 1)))

      CA_ind <- (Beta %*% Mu) %*% t(E_lambda_ind) %*% diag(diag(Z_ind))

      A_ind <- t(matrix(solve(HA_ind, -c(CA_ind)), nrow=2))

      CB_ind <- Z_ind * E_ll_ind + ((Alpha - Eta) * Z_ind + Eta * Z_sq_ind) %*% diag(c(E_lambda_ind)) - t(A_ind)

      B_ind <- t(matrix(solve(HA_ind, -c(CB_ind)), nrow=2))

      E_nl_ind <- A_ind * horizon + B_ind

      P1_ind <- Z_ind * 0.5 * A_ind + t(Z_ind * 0.5 * A_ind)
      P2_ind <- Z_ind * B_ind  + t(Z_ind * B_ind)
      P3_ind <- diag(diag(Z_sq_ind)) %*% diag(c(E_lambda_ind))

      E_NN_ind <- P1_ind * horizon^2 + P2_ind * horizon + P3_ind * horizon

      variance_ind <- matrix(c(1, -1), nrow=1) %*% (P2_ind * horizon + P3_ind * horizon) %*% matrix(c(1, -1))

      (std_ind <- sqrt(variance_ind))

    }

  }
)

